"""
design slidelock probes # by Luyi Tian & Qiyu Gong
"""
import os
from collections import Counter
import pickle
from math import log
import mappy as mp
from Bio.SeqUtils import MeltingTemp as mt
from parse_gene_anno import GeneAnno
from interlap import reduce
import sys
import pandas as pd


probe_len = 52
min_CDS_probes = 510000 #100
GC_range = [0.10,15] #[0.45,0.65]
min_Tm = 10000 #50
max_Tm_diff = 1 #8
num_rep = 10000000 #5
min_ent = 0.1 #0.9

min_dist = 20

# min_dist = 30
# probe_len = 52
# min_Tm = 50
# max_Tm_diff = 8
# GC_range = [0.45, 0.65]
# prime_5end = int(probe_len/2)
# min_CDS_probes = 100
# max_intron_len = 4000  # only search first x/2 and last x/2 bp for introns larger than x bp


_c = {"A":"T","T":"A","G":"C","C":"G","N":"N"}

def _rc(seq):
    return "".join(_c[x] for x in seq[::-1])

def get_fa(fn):
    ch = ""
    seq = []
    for line in open(fn):
        if line[0] == ">":
            if ch != "":
                yield ch, "".join(seq)
            ch = line[1:].strip().split()[0]
            seq = []
        else:
            seq.append(line.strip().upper())
    yield ch, "".join(seq)


def align_seq_to_genome(seq,aln):
    return [(hit.ctg, hit.r_st, hit.r_en, hit.mapq) for hit in aln.map(seq)]



def build_index(fa_f,idx_f):
    mp.Aligner(fn_idx_in=fa_f, preset="sr",k=14,w=2,fn_idx_out=idx_f)


def get_GC_pct(seq):
    return (seq.count("G")+seq.count("C"))/float(len(seq))


def seq_entropy(seq):
    res = 0.
    for st in list(set(seq)):
        p = float(seq.count(st))/len(seq)
        res += -p*log(p)
    return res


def get_seq_complexity(seq, probe_len, prime_5end, num_rep=5, min_ent=0.9):
    prime_5end=int(probe_len/2)
    if seq.find("A"*num_rep)>0 or seq.find("T"*num_rep)>0 or seq.find("C"*num_rep)>0 or seq.find("G"*num_rep)>0:
        return False
    elif seq.find("N")>0:
        return False
    elif seq_entropy(seq[:prime_5end])<min_ent or seq_entropy(seq[-prime_5end:])<min_ent :
        return False
    else:
        return True


def split_list(li, trunk_size=100):
    new_li = []
    new_st = 0
    if len(li)<trunk_size:
        return [li]
    while True:
        if (new_st+trunk_size)<len(li):
            new_li.append(li[new_st:(new_st+trunk_size)])
            new_st += trunk_size
        else:
            new_li.append(li[new_st:] )
            break
    return new_li



def get_reference_paths(base_dir="/Users/gongqiyu/Documents/Probe/design/data/reference"):
    ref_dir = os.path.join(base_dir, "")
    genome = "GRCh38.primary_assembly.genome.fa"
    index = "GRCh38.primary_assembly.genome.fa"
    gtf = "gencode.v39.annotation.gtf"
    return {
        'genome_fa': os.path.join(ref_dir, genome),
        'index_f': os.path.join(ref_dir, index),
        'gtf_f': os.path.join(ref_dir, gtf)
    }


def get_matching_exons(selected_probe, data):
    selected_probe_exon = {}
    
    for gene, probes in selected_probe.items():
        gene_exon = data[data["gene_id"] == gene]
        selected_probe_exon[gene] = []

        for probe in probes:
            probe_position = probe[2]
            matching_exons = gene_exon[(gene_exon['start'] <= probe_position) & (gene_exon['end'] >= probe_position)]
            
            if not matching_exons.empty:
                exon_rank = matching_exons['exon_rank'].iloc[0]
                exon_id = matching_exons['exon_id'].iloc[0]
                selected_probe_exon[gene].append(probe + (exon_rank, exon_id))
    
    return selected_probe_exon


    
    
def _get_probe_from_intervals(an,aln,ix,ge,interval_l, tmp_probe, p_idx, region,result,status_cnt,probe_len, 
                              GC_range, min_Tm, max_Tm_diff, num_rep=5, min_ent=0.9):
    """
    ix: left coordinate of probe sequence
    ge: gene id
    interval_l: InterLap of current gene
    tmp_probe: probe sequence
    p_idx: probe NO. for current gene
    region: the region where probe located
    """

    prime_5end = int(probe_len/2)
    GC_pct = get_GC_pct(tmp_probe)
    seq_complex = get_seq_complexity(tmp_probe, probe_len, prime_5end, num_rep, min_ent)
    assert(len(tmp_probe)==probe_len)
    Tm5 = mt.Tm_NN(tmp_probe[:prime_5end], nn_table=mt.R_DNA_NN1,dnac1=100)
    Tm3 = mt.Tm_NN(tmp_probe[-prime_5end:], nn_table=mt.R_DNA_NN1,dnac1=100)

    if (GC_range[0]<GC_pct<GC_range[1]) and seq_complex and (Tm5>min_Tm) and (Tm3>min_Tm) and (abs(Tm5-Tm3)< max_Tm_diff):
        aln_res = align_seq_to_genome(tmp_probe,aln)
        exclude = True
        if len(aln_res)>0 and ((type(ix)==int and aln_res[0][1] == ix) or (type(ix)==tuple)):
            if len(aln_res)==1:
                status_cnt[ge]["passed_aln_unique"] += 1
                exclude = False
            else:
                exclude = False
                for res in aln_res[1:]:
                    if (res[1],res[2]) in an.chr_interval[res[0]]:
                        status_cnt[ge]["passed_aln_m2gene"] += 1
                        exclude = True
                        break
                if not exclude:
                    status_cnt[ge]["passed_aln_m2genome"] += 1   
        else:
            status_cnt[ge]["passed"] += 1
        if not exclude:
            probe_id = f"{ge}_{region}_{p_idx}"
            if type(interval_l)==int:
                num_overlap = interval_l
            else:
                num_overlap = len(list(interval_l.find((ix, ix+probe_len))))
            result[ge].append((probe_id, tmp_probe, ix, ge, region, num_overlap, GC_pct, seq_entropy(tmp_probe), (Tm5+Tm3)/2.0))
            # probe_id, probe sequence, probe_location, gene_id, region_type, number_of_overlap_transcripts, GC, entropy, Tm
            p_idx += 1
    else:
        status_cnt[ge]["not_passed"] += 1
    return p_idx

def get_exon_probe_list(an, aln,fa_dict, sel_gene=None, probe_len=52, min_CDS_probes=100, GC_range=(0.45,0.65), min_Tm= 50, max_Tm_diff=8, num_rep=5, min_ent=0.9):

    result = {}
    ixx=0
    status_cnt = {}
    for ch in an.chr_to_gene:
        processed_pos = {}
        for ge in an.chr_to_gene[ch]:
            if (sel_gene is not None) and (ge not in sel_gene):
                continue
            stnd = an.chr_to_gene[ch][ge][2]
            p_idx = 1
            result[ge] = []
            ixx += 1
            if ixx % 100==0:
                print(ixx,"genes processed")
            status_cnt[ge] = Counter()
            
            
            if len(an.gene_CDS_interval[ge])>0:
                region = "CDS"
                for ex in an.gene_CDS_interval[ge]:
                    if ex[1]-ex[0]>50:
                        for ix in range(ex[0],ex[1]-probe_len):
                            if ix in processed_pos:
                                continue
                            processed_pos[ix] = 0
                            tmp_probe = fa_dict[ch][ix:(ix+probe_len)]
                            
                            if len(tmp_probe) == 0:
                                continue
                            
                            if stnd=="+":
                                tmp_probe = _rc(tmp_probe) 

                            p_idx = _get_probe_from_intervals(an, aln, ix, ge, an.gene_CDS_interval[ge], tmp_probe, p_idx, region, result, status_cnt, probe_len, GC_range, min_Tm, max_Tm_diff, num_rep, min_ent)
                                      
            if status_cnt[ge]["passed_aln_unique"]+ status_cnt[ge]["passed_aln_m2genome"]<min_CDS_probes:
                region = "exon"
                for ex in an.gene_interval[ge]:
                    if ex[1]-ex[0]>40:
                        for ix in range(ex[0],ex[1]-probe_len):
                            if ix in processed_pos:
                                continue
                            processed_pos[ix] = 0
                            tmp_probe = fa_dict[ch][ix:(ix+probe_len)]
                            
                            if len(tmp_probe) == 0:
                                continue
                            
                            if stnd=="+":
                                tmp_probe = _rc(tmp_probe)
                                
                            p_idx = _get_probe_from_intervals(an, aln, ix, ge, an.gene_interval[ge], tmp_probe, p_idx, region,result,status_cnt, probe_len, GC_range, min_Tm, max_Tm_diff, num_rep, min_ent)
    return status_cnt, result


def exon2intron(ranges):
    itr = []
    if len(ranges)<2:
        return itr
    for ix in range(len(ranges)-1):
        assert(ranges[ix][1]<ranges[ix+1][0])
        itr.append((ranges[ix][1],ranges[ix+1][0]))
    return(itr)


def get_intron_probe_list(an, aln,fa_dict,sel_gene=None, probe_len=52, min_CDS_probes=100, GC_range=(0.45,0.65), min_Tm= 50, max_Tm_diff=8, max_intron_len = 4000):
    result = {}
    ixx=0
    status_cnt = {}
    for ch in an.chr_to_gene:
        processed_pos = {}
        for ge in an.chr_to_gene[ch]:
            if (sel_gene is not None) and (ge not in sel_gene):
                continue
            stnd = an.chr_to_gene[ch][ge][2]
            p_idx = 1
            result[ge] = []
            status_cnt[ge] = Counter()
            print(ge,stnd,len(an.gene_interval[ge]))
            if an.gene_interval[ge].max()>=len(fa_dict[ch]):
                print("exon exceed the chromosome sequence.",an.gene_interval[ge].max(),len(fa_dict[ch]))
                continue
            if len(an.gene_CDS_interval[ge])>1:
                ranges = [(it[0],it[1]) for it in an.gene_CDS_interval[ge]]
                ranges.sort(key=lambda x:x[0])  # sort by left coordinate
                ranges = reduce(ranges)
                if len(ranges)>1:
                    intron_range = exon2intron(ranges)
                    region = "intron"
                    for ex in intron_range:
                        if ex[1]-ex[0]>50:
                            if ex[1]-ex[0]>max_intron_len+1:
                                ix_range = list(range(ex[0],ex[0]+int(max_intron_len/2)-probe_len))
                                ix_range.extend(list(range(ex[1]-int(max_intron_len/2),ex[1]-probe_len)))
                            else:
                                ix_range = range(ex[0],ex[1]-probe_len)
                            for ix in ix_range:
                                if ix in processed_pos:
                                    continue
                                processed_pos[ix] = 0
                                tmp_probe = fa_dict[ch][ix:(ix+probe_len)]
                                if stnd=="+":
                                    tmp_probe = _rc(tmp_probe)
                                p_idx = _get_probe_from_intervals(an, aln, ix, ge, 1, tmp_probe, p_idx, region, result, status_cnt, probe_len, GC_range, min_Tm, max_Tm_diff) 
            if len(an.gene_interval[ge])>1 and status_cnt[ge]["passed_aln_unique"]+ status_cnt[ge]["passed_aln_m2genome"]<min_CDS_probes:
                ranges = [(it[0],it[1]) for it in an.gene_interval[ge]]
                ranges.sort(key=lambda x:x[0])  # sort by left coordinate
                ranges = reduce(ranges)
                if len(ranges)>1:
                    intron_range = exon2intron(ranges)
                    region = "intron"
                    for ex in intron_range:
                        if ex[1]-ex[0]>50:
                            if ex[1]-ex[0]>max_intron_len+1:
                                ix_range = list(range(ex[0],ex[0]+int(max_intron_len/2)-probe_len))
                                ix_range.extend(list(range(ex[1]-int(max_intron_len/2),ex[1]-probe_len)))
                            else:
                                ix_range = range(ex[0],ex[1]-probe_len)
                            for ix in ix_range:
                                if ix in processed_pos:
                                    continue
                                processed_pos[ix] = 0
                                tmp_probe = fa_dict[ch][ix:(ix+probe_len)]
                                if stnd=="+":
                                    tmp_probe = _rc(tmp_probe)
                                p_idx = _get_probe_from_intervals(an,aln,ix,ge,1, tmp_probe, p_idx, region,result,status_cnt, probe_len, GC_range, min_Tm, max_Tm_diff)
            print(ge,status_cnt[ge])
    return status_cnt, result


def get_junc_probe_list(an, aln,fa_dict,sel_gene=None, probe_len=52, GC_range=(0.45,0.65), min_Tm= 50, max_Tm_diff=8):
    prime_5end=int(probe_len/2)
    result = {}
    ixx=0
    status_cnt = {}
    for ch in an.chr_to_gene:
        for ge in an.chr_to_gene[ch]:
            if (sel_gene is not None) and (ge not in sel_gene):
                continue
            stnd = an.chr_to_gene[ch][ge][2]
            p_idx = 1
            result[ge] = []
            print(ge,an.gene_interval[ge])
            status_cnt[ge] = Counter()
            if len(an.gene_interval[ge])>1:
                junc_cnt = Counter()
                for tr in an.gene_to_transcript[ge]:
                    intron_range = exon2intron(an.transcript_to_exon[tr]._iset)
                    for r in intron_range:
                        junc_cnt[r] += 1
                region = "ee_junc"
                if len(junc_cnt)==0:
                    continue
                for junc in junc_cnt:
                    print(junc,junc_cnt[junc])
                    for dist2junc in range(-int(prime_5end/2),int(prime_5end/2)):
                        tmp_probe = fa_dict[ch][(junc[0]-prime_5end-dist2junc):junc[0]]+\
                            fa_dict[ch][junc[1]:(junc[1]+prime_5end-dist2junc)]
                        if stnd=="+":
                            tmp_probe = _rc(tmp_probe)
                        p_idx = _get_probe_from_intervals(an,aln,junc,ge,junc_cnt[junc], tmp_probe, p_idx, region,result,status_cnt, probe_len, GC_range, min_Tm, max_Tm_diff)
            print(status_cnt[ge])
    return status_cnt, result


def check_probes(gtf_f,probe_dict, out_file):
    an = GeneAnno(gtf_f)
    out_stat = open(out_file,"w")
    out_stat.write("probe_id,GC_pct,entropy,Tm,alignment_status,top_matching_gene\n")
    aln = mp.Aligner(fn_idx_in=index_f, preset="sr",best_n=4,min_chain_score=16)
    for pid in probe_dict:
        seq = probe_dict[pid]
        GC = get_GC_pct(seq)
        ent = seq_entropy(seq)
        Tm = mt.Tm_NN(seq)
        rec = [hit for hit in aln.map(seq)]
        aln_status = "unaligned"
        top_gene = "NA"
        if len(rec)>0:
            top_aln = list(an.chr_interval[rec[0].ctg].find((rec[0].r_st, rec[0].r_en)))
            if len(top_aln)>0:
                top_gene = top_aln[0][2]
            else:
                print(pid,rec[0])
                print("\t",list(an.chr_interval[rec[0].ctg].closest((rec[0].r_st, rec[0].r_en))))
            if len(rec)==1:
                aln_status = "unique_align"
            elif len(rec)>1:
                found = False
                for res in rec[1:]:
                    if (res.r_st, res.r_en) in an.chr_interval[res.ctg]:
                        aln_status = "multi_aln_gene"
                        found = True
                        break
                if not found:
                    aln_status = "multi_aln_genome"
        out_stat.write("{},{},{},{},{},{}\n".format(pid,GC,ent,Tm,aln_status,top_gene))
    out_stat.close()
            


# if __name__ == '__main__':
#     """
#     ref_dir = "/broad/thechenlab/Luyi/reference/mouse"
#     genome_fa = os.path.join(ref_dir,"GRCm39.primary_assembly.genome.fa")
#     index_f = os.path.join(ref_dir,"GRCm39.primary_assembly.genome.mmi")
#     gtf_f = os.path.join(ref_dir,"gencode.vM27.annotation.gtf")

#     ref_dir = "/broad/thechenlab/Luyi/reference/human"
#     genome_fa = os.path.join(ref_dir,"GRCh38.primary_assembly.genome.fa")
#     index_f = os.path.join(ref_dir,"GRCh38.primary_assembly.genome.fa")
#     gtf_f = os.path.join(ref_dir,"gencode.v39.annotation.gtf")
#     """
#     out_dir = os.path.join( "/broad/thechenlab/Chenlei/slidelock/probe_design",sys.argv[1])
#     sel_ge2 = [it.strip().split(",")[0] for it in open(os.path.join(out_dir,"markers.csv")).readlines()[1:]]
#     if "ENSMUSG" in sel_ge2[0]:
#         ref_dir = "/broad/thechenlab/Chenlei/slidelock/reference/mouse"
#         genome_fa = os.path.join(ref_dir,"GRCm39.primary_assembly.genome.fa")
#         index_f = os.path.join(ref_dir,"GRCm39.primary_assembly.genome.mmi")
#         gtf_f = os.path.join(ref_dir,"gencode.vM27.annotation.gtf")
#     else:
#         ref_dir = "/broad/thechenlab/Chenlei/slidelock/reference/human"
#         genome_fa = os.path.join(ref_dir,"GRCh38.primary_assembly.genome.fa")
#         index_f = os.path.join(ref_dir,"GRCh38.primary_assembly.genome.fa")
#         gtf_f = os.path.join(ref_dir,"gencode.v39.annotation.gtf")
#     print(out_dir)
#     fa_dict = {}
#     for c in get_fa(genome_fa):
#         fa_dict[c[0]] = c[1]
#         print(c[0],len(c[1]))
#     aln = mp.Aligner(fn_idx_in=index_f, preset="sr",best_n=4,min_chain_score=16)
#     an = GeneAnno(gtf_f)

#     #sel_ge = [it.strip().split(",")[0] for it in open(os.path.join(out_dir,"gene_list_large_withname_combine.tsv")).readlines()[1:]]
#     #sel_ge1 = [it.strip().split(",")[0] for it in open(os.path.join(out_dir,"gene_marker_withname.tsv")).readlines()[1:]]
    
#     #status_cnt1, probe_res1 = get_exon_probe_list(an, aln,fa_dict,sel_gene=sel_ge2)
#     status_cnt, probe_res_jun = get_junc_probe_list(an, aln,fa_dict,sel_gene=sel_ge2)
#     #status_cnt2, probe_res2 = get_intron_probe_list(an, aln,fa_dict,sel_gene=sel_ge2)

#     #pickle.dump(probe_res2, open(os.path.join(out_dir,"probe_intron.pkl"),"wb"))
#     #pickle.dump(probe_res1, open(os.path.join(out_dir,"probe_exon.pkl"),"wb"))
#     pickle.dump(probe_res_jun, open(os.path.join(out_dir,"probe_junction.pkl"),"wb"))
    
    

# def symbol_to_ens_biomart(gene_symbol, 
#                           local_path='data/reference/biomart_dataset_human.txt'):
#     if not os.path.exists(local_path):
#         download_biomart_dataset(local_path)
#     ens_ids = []
#     with open(local_path, 'r') as file:
#         for line in file:
#             parts = line.strip().split("\t")
#             if parts[0] == gene_symbol:
#                 ens_ids.append(parts[1])
#     return ens_ids

# symbol_to_ens_biomart('TRBC1')




def save_design_probes(data_folder, probe_request, selected_probe_exon, const_3, const_5_1, const_5_2, save_top=3, filter_same_exon=True):
    data_list = []
    for ge in selected_probe_exon:
        for rec in selected_probe_exon[ge]:
            probe_homo_seq = rec[1].lower()
            probe_3 = probe_homo_seq[0:25] + const_3
            probe_5 = const_5_1 + probe_homo_seq[27:] + const_5_2
            data = {
                "probe_id": rec[0],
                "probe_homology_sequence52bp": rec[1],
                "probe_start": rec[2],
                "gene_id": rec[3],
                "region_type": rec[4],
                "number_of_overlap_transcripts": rec[5],
                "GC": rec[6],
                "entropy": rec[7],
                "Tm": rec[8],
                "exon_rank": rec[9],
                "exon_id": rec[10],
                "probe_3": probe_3,
                "probe_5": probe_5
            }
            data_list.append(data)

    all_probes = pd.DataFrame(data_list)
    merged_df = all_probes.merge(probe_request[['gene_id', 'gene symbol']], on='gene_id', how='left')
    # merged_df['gene symbol'].fillna('', inplace=True)
    merged_df['gene symbol'] = merged_df['gene symbol'].fillna('')
    merged_df = merged_df.drop_duplicates()
    cols = merged_df.columns.tolist()
    merged_df = merged_df[[cols[-1]] + cols[:-1]]
    merged_df.to_csv(data_folder + 'Probe_exon_design_all_test_0.csv', index=False)
    
    # filter probes in the same exons
    def find_closest_to_55(group):
        group['distance'] = abs(group['Tm'] - 55)
        return group.nsmallest(6, 'distance')
    
    filtered_df = merged_df[merged_df['Tm'] >= 55]
    result_df = filtered_df.groupby('exon_id').filter(lambda x: len(x) >= 1).groupby('exon_id').apply(find_closest_to_55)
    result_df.to_csv(data_folder + 'Probe_exon_design_all_test_1.csv', index=False)

    if filter_same_exon:
        non_duplicates = filtered_df.groupby('exon_id').filter(lambda x: len(x) == 1)
        result_df = pd.concat([result_df, non_duplicates])
        result_df = result_df.drop(columns='distance').reset_index(drop=True)
    result_df.to_csv(data_folder + 'Probe_exon_design_all.csv', index=False)

    # save top probes
    if save_top > 0:
        selected = merged_df.copy()
        for gene in all_probes['gene_id'].unique():
            gene_list = all_probes[all_probes['gene_id'] == gene]
            gene_list = gene_list.sort_values('exon_rank', ascending=True)
            
            if len(gene_list) <= save_top:
                continue

            gene_list_keep = gene_list.groupby('exon_rank').head(1)
            
            if len(gene_list_keep) < save_top:
                gene_list_keep = gene_list[0:save_top]
            else:
                gene_list_keep = gene_list_keep.head(save_top)
                
            gene_list_del = gene_list[~gene_list['probe_id'].isin(gene_list_keep['probe_id'])]
            selected = selected[~selected['probe_id'].isin(gene_list_del['probe_id'])]

        df = selected
        df['gene_symbol_unique'] = df.groupby('gene symbol').cumcount() + 1
        df['gene_symbol_unique'] = df['gene symbol'] + '_' + df['gene_symbol_unique'].astype(str)
        df = df[['gene_symbol_unique'] + [col for col in df if col != 'gene_symbol_unique']]
        df = df.sort_values(['gene symbol', 'gene_id','exon_rank', 'gene_symbol_unique'])
        # df = df.drop_duplicates(subset=['probe_5'], keep='first')
        df = df.groupby('gene symbol').head(save_top).reset_index(drop=True)
        df.to_csv(data_folder + 'Probe_exon_design_top.csv', index=False)




def update_probe_request_gene_ids(probe_request, data, data_folder):

    probe_request.rename(columns={"Unnamed: 1": "gene_id"}, inplace=True)
    
    for index, row in probe_request.iterrows():
        if pd.isnull(row['gene_id']):
            matched_row = data[data['gene_name'] == row['gene symbol']]
            if not matched_row.empty:
                probe_request.at[index, 'gene_id'] = matched_row['gene_id'].iloc[0]
    
    na_rows = probe_request[probe_request['gene_id'].isna()]
    
    if na_rows.empty:
        probe_request.to_csv(data_folder + "Updated_probe_request.csv", index=False)
        return []
    else:
        return list(na_rows['gene symbol'])
    
    
    


def save_design_probes_intron(data_folder, probe_request, selected_probe_exon, const_3, const_5_1, const_5_2):

    data_list = []
    for ge in selected_probe_exon:
        for rec in selected_probe_exon[ge]:
            probe_homo_seq = rec[1].lower()
            probe_3 = probe_homo_seq[0:25] + const_3
            probe_5 = const_5_1 + probe_homo_seq[27:] + const_5_2
            data = {
                "probe_id": rec[0],
                "probe_homology_sequence52bp": rec[1],
                "probe_start": rec[2],
                "gene_id": rec[3],
                "region_type": rec[4],
                "number_of_overlap_transcripts": rec[5],
                "GC": rec[6],
                "entropy": rec[7],
                "Tm": rec[8],
                "probe_3": probe_3,
                "probe_5": probe_5
            }
            data_list.append(data)

    all_probes = pd.DataFrame(data_list)
    merged_df = all_probes.merge(probe_request[['gene_id', 'gene symbol']], on='gene_id', how='left')
    merged_df['gene symbol'].fillna('', inplace=True)
    merged_df = merged_df.drop_duplicates()
    cols = merged_df.columns.tolist()
    merged_df = merged_df[[cols[-1]] + cols[:-1]]
    filtered_df = merged_df[merged_df['Tm'] >= 55]
    filtered_df.to_csv(data_folder + 'Probe_intron_design_all.csv', index=False)
    return None
