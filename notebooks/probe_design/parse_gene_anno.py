#parse gene annotation # # by Luyi Tian
from collections import namedtuple, defaultdict
import copy
from interlap import InterLap

# Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)
Pos = namedtuple('Pos', ["chr",'start', 'end', "strand","parent_id"])
def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""#
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        ret[key] = value
    return ret


### parse gtf file

# Initialized GeneInfo named tuple. Note: namedtuple is immutable
# gtf and gff have the same column specification so we can use `GFFRecord` for gtf files
def parseGTFAttributes(attributeString):
    """Parse the GTF attribute column and return a dict"""#
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        if len(attribute)>0:
            items = attribute.split("\"")
            if len(items)<2:
                continue
            key = items[0].strip()
            value = items[1].strip()
            ret[key] = value
    return ret


def parseGFF3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.
    can also be used to parse gtf file
    Supports transparent gzip decompression.
    """
    #Parse with transparent decompression
    openFunc = gzip.open if filename.endswith(".gz") else open
    attrFunc = parseGTFAttributes if filename.endswith("gtf.gz") or filename.endswith(".gtf") else parseGFFAttributes
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            if len(parts) != len(gffInfoFields):
                print(line)
                print(parts)
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else parts[0],
                "source": None if parts[1] == "." else parts[1],
                "type": None if parts[2] == "." else parts[2],
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else parts[6],
                "phase": None if parts[7] == "." else parts[7],
                "attributes": attrFunc(parts[8])
            }
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield GFFRecord(**normalizedInfo)
######


def guess_annotation_source(gff_f):
    idx = 0
    for line in open(gff_f):
        idx += 1
        if "GENCODE" in line:
            print("parse GENCODE annotation.")
            return "GENCODE"
        elif "1\tEnsembl" in line:
            print("parse Ensembl annotation.")
            return "Ensembl"
        if idx > 1000:  # should be within the first 1000 lines
            return "Ensembl"
    return "Ensembl"  # default is Ensembl


def clean_geneid(gene_id):
    return gene_id.split(".")[0]

def _parse_gff_tree(gff_f):
    chr_to_gene = {}
    transcript_dict = {}
    gene_to_transcript = {}
    transcript_to_exon = {}
    a_s = guess_annotation_source(gff_f)
    if a_s == "Ensembl":
        for rec in parseGFF3(gff_f):
            #if rec.seqid != "1":  # for test
            #    break
            if "gene_id" in rec.attributes:
                chr_to_gene.setdefault(rec.seqid,[]).append(clean_geneid(rec.attributes["gene_id"]))
            if "Parent" in rec.attributes and (rec.attributes["Parent"].split(':')[0] == "gene"):  # transcript
                gene_id = clean_geneid(rec.attributes["Parent"].split(':')[1])
                gene_to_transcript.setdefault(gene_id, []).append(rec.attributes["transcript_id"])
                transcript_dict[rec.attributes["transcript_id"]] = Pos(rec.seqid, rec.start-1, rec.end, rec.strand, gene_id)  # `-1` convert 1 based to 0 based
            elif rec.type == "exon":
                if rec.attributes["Parent"].split(':')[0] != "transcript":
                    print(rec)
                    print("format error.")
                    raise Exception
                transcript_to_exon.setdefault(rec.attributes["Parent"].split(':')[1], []).append([rec.start-1, rec.end])  # `-1` convert 1 based to 0 based
    elif a_s == "GENCODE":
        for rec in parseGFF3(gff_f):
            if rec.type == "gene":
                chr_to_gene.setdefault(rec.seqid,[]).append(clean_geneid(rec.attributes["gene_id"]))
            elif rec.type == "transcript":
                gene_id = clean_geneid(rec.attributes["Parent"])
                if gene_id not in chr_to_gene[rec.seqid]:
                    chr_to_gene.setdefault(rec.seqid,[]).append(gene_id)
                gene_to_transcript.setdefault(gene_id, []).append(rec.attributes["transcript_id"])
                transcript_dict[rec.attributes["transcript_id"]] = Pos(rec.seqid, rec.start-1, rec.end, rec.strand, gene_id)  # `-1` convert 1 based to 0 based
            elif rec.type == "exon":
                transcript_to_exon.setdefault(rec.attributes["Parent"], []).append([rec.start-1, rec.end])  # `-1` convert 1 based to 0 based
    for tr in transcript_to_exon:
        transcript_to_exon[tr].sort(key=lambda x: x[0])  # the GENCODE annotation might be un-ordered.
        if len(transcript_to_exon[tr])>1 and transcript_to_exon[tr][0][0] == transcript_to_exon[tr][1][0]:  # for genes in XY, there might be duplicates.
            new_ex = [transcript_to_exon[tr][0]]
            for ex in transcript_to_exon[tr]:
                if ex[0] != new_ex[-1][0] and ex[1] != new_ex[-1][0]:
                    new_ex.append(ex)
            transcript_to_exon[tr] = new_ex
    for ge in gene_to_transcript:
        gene_to_transcript[ge] = list(set(gene_to_transcript[ge]))
    return chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon


def _parse_gtf_tree(gtf_f):
    chr_to_gene = {}
    transcript_dict = {}
    gene_to_transcript = {}
    transcript_to_exon = defaultdict(InterLap)
    transcript_to_CDS = defaultdict(InterLap)
    chr_interval = defaultdict(InterLap)
    for rec in parseGFF3(gtf_f):
        if rec.type == "gene":
            chr_to_gene.setdefault(rec.seqid,{})[clean_geneid(rec.attributes["gene_id"])] = (rec.start-1, rec.end, rec.strand)
            chr_interval[rec.seqid].add((rec.start-1, rec.end, clean_geneid(rec.attributes["gene_id"])))
        elif rec.type == "transcript":
            if "gene_id" not in rec.attributes:
                Warning("transcript dont have `gene_id` attributes: {}".format(rec.attributes))
                continue
            gene_id = clean_geneid(rec.attributes["gene_id"])
            chr_interval[rec.seqid].add((rec.start-1, rec.end, gene_id))
            if gene_id not in chr_to_gene[rec.seqid]:
                chr_to_gene.setdefault(rec.seqid,{})[gene_id] = (rec.start-1, rec.end, rec.strand)
            gene_to_transcript.setdefault(gene_id, []).append(rec.attributes["transcript_id"])
            transcript_dict[rec.attributes["transcript_id"]] = Pos(rec.seqid, rec.start-1, rec.end, rec.strand, gene_id)  # `-1` convert 1 based to 0 based
        elif rec.type == "exon":
            if "gene_id" not in rec.attributes:
                Warning("exon dont have `gene_id` attributes: {}".format(rec.attributes))
                continue
            if rec.seqid not in chr_to_gene or clean_geneid(rec.attributes["gene_id"]) not in chr_to_gene[rec.seqid]:
                chr_to_gene.setdefault(rec.seqid,{})[clean_geneid(rec.attributes["gene_id"])] = (rec.start-1, rec.end, rec.strand)
            if "transcript_id" in rec.attributes:
                if rec.attributes["transcript_id"] not in transcript_dict:
                    gene_to_transcript.setdefault(clean_geneid(rec.attributes["gene_id"]),[]).append(rec.attributes["transcript_id"])
                    transcript_dict[rec.attributes["transcript_id"]] = Pos(rec.seqid, -1, 1, rec.strand, clean_geneid(rec.attributes["gene_id"]))
                transcript_to_exon[rec.attributes["transcript_id"]].add((rec.start-1, rec.end))  # `-1` convert 1 based to 0 based
            else:
                Warning("exon dont have `transcript_id` attributes: {}".format(rec.attributes))
        elif rec.type == "CDS":
            if "gene_id" not in rec.attributes:
                Warning("CDS dont have `gene_id` attributes: {}".format(rec.attributes))
                continue
            if rec.seqid not in chr_to_gene or rec.attributes["gene_id"] not in chr_to_gene[rec.seqid]:
                chr_to_gene.setdefault(rec.seqid,{})[gene_id] = (rec.start-1, rec.end, rec.strand)
            if "transcript_id" in rec.attributes:
                if rec.attributes["transcript_id"] not in transcript_dict:
                    gene_to_transcript.setdefault(rec.attributes["gene_id"],[]).append(rec.attributes["transcript_id"])
                    transcript_dict[rec.attributes["transcript_id"]] = Pos(rec.seqid, -1, 1, rec.strand, clean_geneid(rec.attributes["gene_id"]))
                transcript_to_CDS[rec.attributes["transcript_id"]].add((rec.start-1, rec.end))  # `-1` convert 1 based to 0 based
            else:
                Warning("exon dont have `transcript_id` attributes: {}".format(rec.attributes))
    for tr in transcript_to_exon:
       transcript_dict[tr] = Pos(transcript_dict[tr].chr, transcript_to_exon[tr].min(), transcript_to_exon[tr].max(), transcript_dict[tr].strand, transcript_dict[tr].parent_id)
    for ge in gene_to_transcript:
        gene_to_transcript[ge] = list(set(gene_to_transcript[ge]))
    return chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon, transcript_to_CDS, chr_interval

# naive warpper
def parse_gff_tree(gff_f):
    if gff_f.endswith("gtf.gz") or gff_f.endswith(".gtf"):
        return _parse_gtf_tree(gff_f)
    else:
        return _parse_gff_tree(gff_f)


def get_gene_flat(gene_to_transcript, transcript_to_exon):
    gene_dict = defaultdict(InterLap)
    for g in gene_to_transcript:
        for tr in gene_to_transcript[g]:
            for ex in transcript_to_exon[tr]:
                gene_dict[g].add(ex)
    return gene_dict


def get_gene_conserve(gene_to_transcript, transcript_to_exon):
    pass


class GeneAnno(object):
    """
    store gene annotation from the gtf/gff files
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon, transcript_to_CDS, chr_interval
    """
    def __init__(self,gff_f):
        res = parse_gff_tree(gff_f)
        self.chr_to_gene = res[0]
        self.transcript_dict = res[1]
        self.gene_to_transcript = res[2]
        self.transcript_to_exon = res[3]
        self.transcript_to_CDS = res[4]
        self.chr_interval = res[5]
        self.gene_interval = get_gene_flat(self.gene_to_transcript, self.transcript_to_exon)
        self.gene_CDS_interval = get_gene_flat(self.gene_to_transcript, self.transcript_to_CDS)


