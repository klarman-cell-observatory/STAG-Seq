"""
design slidelock probes
"""
import os
import sys
from collections import Counter
from math import log
import mappy as mp
from Bio.SeqUtils import MeltingTemp as mt
from parse_gene_anno import GeneAnno
from interlap import reduce


flanking_seq_new =   "TGGAATTCTCGGGTGCCACCGGTTTTTTTTTTTTTTTTTTTTNNNNNNNNNN"
flanking_seq_new1 = ["TGGAATTCTCGGGTGCCACCGGTTTTTTTTTTTTTTTTTTTANNNNNNNNNN",
                     "TGGAATTCTCGGGTGCCACCGGTTTTTTTTTTTTTTTTTTTCNNNNNNNNNN",
                     "TGGAATTCTCGGGTGCCACCGGTTTTTTTTTTTTTTTTTTTGNNNNNNNNNN"]
flanking_seq_intron = ["TGGAATTCTCGGGTGCCAAGCTTTTTTTTTTTTTTTTTTTTANNNNNNNNNN",
                     "TGGAATTCTCGGGTGCCAAGCTTTTTTTTTTTTTTTTTTTTCNNNNNNNNNN",
                     "TGGAATTCTCGGGTGCCAAGCTTTTTTTTTTTTTTTTTTTTGNNNNNNNNNN"]
flanking_seq_old = "NNNNNNNNNNTGGAATTCTCGGGTGCCACCGGTTTTTTTTTTTTTTTTTTTT"
sel_junc = ["TC","TA","CT","CA","TT"]
sel_junc_exclude = ["GC","GG","CG"]

def write_probe_seq(out_f,flanking_seq,probe_dict):
    pool_size=384
    with open(out_f,"w") as f:
        pool_ix = 1
        ix = 1
        for ge in probe_dict:
            for rec in probe_dict[ge]:
                probe_homo_seq = rec[1].lower()
                probe_seq = probe_homo_seq[16:]+flanking_seq+probe_homo_seq[:16]
                probe_seq = "/5Phos/"+probe_seq
                f.write("pool_"+str(pool_ix)+","+rec[0]+','+probe_seq+"\n")
                if ix % pool_size==0:
                    pool_ix += 1
                ix += 1

def write_probe_metadata(out_f,flanking_seq,selected_probe_dict):
    with open(out_f,"w") as f:
        f.write("probe_id,probe_homology_sequence,probe_start,gene_id,region_type,number_of_overlap_transcripts,GC,entropy,Tm,final_sequence\n")
        for ge in selected_probe_dict:
            # probe_id, probe sequence, probe_location, gene_id, region_type, number_of_overlap_transcripts, GC, entropy, Tm, probe_final_seq
            for rec in selected_probe_dict[ge]:
                probe_homo_seq = rec[1].lower()
                probe_seq = probe_homo_seq[16:]+flanking_seq+probe_homo_seq[:16]
                probe_seq = "/5Phos/"+probe_seq
                fmt = ",".join(str(it) for it in rec)
                f.write(fmt+","+probe_seq+"\n")


def write_probe_metadata1(out_f,flanking_seq,selected_probe_dict,known_probe=None):
    flanking_seq_new_rmT = "TGGAATTCTCGGGTGCCACCGGTTTTTTTTTTTTTTTTTTTVNNNNNNNNNN"
    ii = 0
    found_probe = []
    with open(out_f,"w") as f:
        f.write("probe_id,probe_homology_sequence16bp,probe_start,gene_id,region_type,number_of_overlap_transcripts,GC,entropy,Tm,probe_homology_sequence14bp,final_sequence14bp,final_sequence_NV,final_sequence_upper\n")
        for ge in selected_probe_dict:
            # probe_id, probe sequence, probe_location, gene_id, region_type, number_of_overlap_transcripts, GC, entropy, Tm, probe_final_seq
            for rec in selected_probe_dict[ge]:
                if (known_probe is None) or (rec[0] in known_probe):
                    found_probe.append(rec[0])
                    probe_homo_seq = rec[1].lower()
                    probe_seq = probe_homo_seq[16:]+flanking_seq[ii % 3]+probe_homo_seq[:16]
                    probe_seq = "/5Phos/"+probe_seq
                    probe_seq14bp = probe_homo_seq[16:30]+flanking_seq[ii % 3]+probe_homo_seq[2:16]
                    probe_seq14bp = "/5Phos/"+probe_seq14bp
                    probe_seqrmT = probe_homo_seq[16:30]+flanking_seq_new_rmT+probe_homo_seq[2:16]
                    probe_seqrmT = "/5Phos/"+probe_seqrmT
                    fmt = ",".join(str(it) for it in rec)
                    f.write(fmt+","+rec[1][2:30]+","+probe_seq14bp+","+probe_seqrmT+"\n")
                    ii += 1
        if (known_probe is not None):
            probe_not_found = [it for it in known_probe if it not in found_probe]
            print("probe not found:",probe_not_found)

def filter_probes(probe_dict,max_exon_probe,min_exon_probe,min_dist,prefer_region="CDS"):
    selected_probe_dict = {}
    noprobe_gene = []
    for ge in probe_dict:
        if len(probe_dict[ge])==0:
            noprobe_gene.append(ge)
            continue
        filtered_list = []
        selected_probe_dict[ge] = []
        for rec in probe_dict[ge]:
            if (rec[1][15:17] not in sel_junc_exclude) and rec[4]==prefer_region and rec[5]>1:  # prefered region, like CDS
                if len(filtered_list)>0 and rec[2]-filtered_list[-1][2]<min_dist:
                    if abs(filtered_list[-1][6]-0.55)>abs(rec[6]-0.55):
                        del filtered_list[-1]
                        filtered_list.append(rec)
                else:
                    filtered_list.append(rec)
        if len(filtered_list)==0:
            for rec in probe_dict[ge]:
                if (rec[1][15:17] not in sel_junc_exclude):
                    if len(filtered_list)>0 and rec[2]-filtered_list[-1][2]<min_dist:
                        if abs(filtered_list[-1][6]-0.55)>abs(rec[6]-0.55):
                            del filtered_list[-1]
                            filtered_list.append(rec)
                    else:
                        filtered_list.append(rec)
        if len(filtered_list)>max_exon_probe:
            gc_sort = [(rec,abs(rec[6]-0.55)) for rec in filtered_list if rec[1][15:17] in sel_junc]
            gc_sort2 = [(rec,abs(rec[6]-0.55)) for rec in filtered_list if rec[1][15:17] not in sel_junc]
            if len(gc_sort)>0:
                gc_sort.sort(key=lambda x:x[1])
            if len(gc_sort2)>0:
                gc_sort2.sort(key=lambda x:x[1])
                gc_sort.extend(gc_sort2)
            filtered_list = [it[0] for it in gc_sort][:max_exon_probe]
        elif len(filtered_list)<min_exon_probe:
            gc_sort = [(rec,abs(rec[6]-0.55)) for rec in probe_dict[ge] if rec[1][15:17] in sel_junc]
            gc_sort2 = [(rec,abs(rec[6]-0.55)) for rec in probe_dict[ge] if rec[1][15:17] not in sel_junc]
            if len(gc_sort)>0:
                gc_sort.sort(key=lambda x:x[1])
            if len(gc_sort2)>0:
                gc_sort2.sort(key=lambda x:x[1])
                gc_sort.extend(gc_sort2)
            for s_rec in gc_sort:
                if len(filtered_list)==0:
                    if s_rec[0] not in filtered_list:
                        filtered_list.append(s_rec[0])
                #elif min(abs(it[2]-s_rec[0][2]) for it in filtered_list)>min_dist:
                else:
                    if s_rec[0] not in filtered_list:
                        filtered_list.append(s_rec[0])
                    if len(filtered_list)>=max_exon_probe:
                        break
        if len(filtered_list)>0:
            selected_probe_dict[ge] = filtered_list
        else:
            print("NO probe after filtering.")
        print("\t",ge,len(filtered_list),len(probe_dict[ge]))
    print("Genes with no probe designed:", len(noprobe_gene)," / ",len(probe_dict))
    print("\t",noprobe_gene)
    return selected_probe_dict


# if __name__ == '__main__':
#     import pickle
#     """
#     out_dir = "/broad/thechenlab/Luyi/slidelock/SlidePipe/probe_design/human_4k"
#     probe_dict = pickle.load(open(os.path.join(out_dir, "probe_exon.pkl"), "rb") )
#     junction_dict = pickle.load(open(os.path.join(out_dir, "probe_junction.pkl"), "rb") )
#     #for ge in probe_dict:
#     #    if ge in junction_dict:
#     #        for rec in junction_dict[ge]:
#     #            if rec not in probe_dict[ge]:
#     #                probe_dict[ge].append(rec)
#     max_exon_probe = 4
#     min_exon_probe = 4
#     min_dist = 50
#     selected_probe_dict = filter_probes(probe_dict,max_exon_probe,min_exon_probe,min_dist)
#     write_probe_seq(os.path.join(out_dir,"exon_probes_finalsequence.txt"),flanking_seq_new,selected_probe_dict)
#     write_probe_metadata(os.path.join(out_dir,"exon_probes_homology_metadata.txt"),flanking_seq_new,selected_probe_dict)
#     """
    
#     out_dir = os.path.join("/broad/thechenlab/Chenlei/slidelock/probe_design", sys.argv[1])
#     probe_dict1 = pickle.load(open(os.path.join(out_dir, "probe_exon.pkl"), "rb") )
#     #junction_dict = pickle.load(open(os.path.join(out_dir, "probe_LN_res.pkl"), "rb") )
#     #for ge in probe_dict:
#     #    if ge in junction_dict:
#     #        for rec in junction_dict[ge]:
#     #            if rec not in probe_dict[ge]:
#     #                probe_dict[ge].append(rec)
#     max_exon_probe = int(sys.argv[2])
#     min_exon_probe = int(sys.argv[3])
#     min_dist = 30
#     selected_probe_dict1 = filter_probes(probe_dict1,max_exon_probe,min_exon_probe,min_dist)
    
#     #write_probe_seq(os.path.join(out_dir,"exon_probes_finalsequence.csv"),flanking_seq_new,selected_probe_dict)
#     #write_probe_metadata(os.path.join(out_dir,"exon_probes_homology_metadata.csv"),flanking_seq_new,selected_probe_dict)
#     write_probe_metadata1(os.path.join(out_dir,"exon_probes_homology_metadata_newdesign.csv"),flanking_seq_new1,selected_probe_dict1)
#     #known_probe = [it.strip().split(",")[0] for it in open(os.path.join(out_dir,"top_probe_selected.csv")).readlines()[1:]]
#     #write_probe_metadata1(os.path.join(out_dir,"exon_probes_homology_metadata_newdesign.csv"),flanking_seq_new1,probe_dict)
    
#     probe_dict2 = pickle.load(open(os.path.join(out_dir, "probe_intron.pkl"), "rb") )
#     selected_probe_dict2 = filter_probes(probe_dict2,max_exon_probe,min_exon_probe,min_dist)
#     write_probe_metadata1(os.path.join(out_dir,"intron_probes_homology_metadata_newdesign.csv"),flanking_seq_intron,selected_probe_dict2)
#     """
#     out_dir = "/broad/thechenlab/Luyi/slidelock/SlidePipe/probe_design/mouse_intron"
#     probe_dict = pickle.load(open(os.path.join(out_dir, "probe_intron.pkl"), "rb") )
#     max_exon_probe = 32
#     min_exon_probe = 32ish -l h_vmem=20G
#     min_dist = 50
#     selected_probe_dict = filter_probes(probe_dict,max_exon_probe,min_exon_probe,min_dist)
#     write_probe_seq(os.path.join(out_dir,"intron_probes_finalsequence.txt"),flanking_seq_new,selected_probe_dict)
#     write_probe_metadata(os.path.join(out_dir,"intron_probes_homology_metadata.txt"),flanking_seq_new,selected_probe_dict)

#     probe_dict = pickle.load(open(os.path.join(out_dir, "probe_snRNA.pkl"), "rb") )
#     max_exon_probe = 1
#     min_exon_probe = 1
#     min_dist = 50
#     selected_probe_dict = filter_probes(probe_dict,max_exon_probe,min_exon_probe,min_dist)
#     write_probe_seq(os.path.join(out_dir,"snRNA_probes_finalsequence.txt"),flanking_seq_new,selected_probe_dict)
#     write_probe_metadata(os.path.join(out_dir,"snRNA_probes_homology_metadata.txt"),flanking_seq_new,selected_probe_dict)
#     """
