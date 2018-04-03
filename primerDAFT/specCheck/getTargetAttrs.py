from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from primerDAFT.specCheck.getMaskedSeq import getMaskedSeq

def getTargetAttrs(target,side,idx,data,side_cols,target_cols,pysam_fasta):
    primer_seq = data["result"]["pairs"][idx]["PRIMER_"+side.upper()]["SEQUENCE"]
    attr_names=["sseqid","sstart","send","strand"]
    json_names=["chr","start","end","strand"]


    cols = [x+"_"+side for x in side_cols] + target_cols
    attrs = [getattr(target,x) for x in cols]
    attr_dict=dict(zip(side_cols + target_cols,attrs))

    out_dict = dict(zip(json_names,[attr_dict.get(x) for x in attr_names]))
    out_dict["start"] = int(out_dict["start"])
    out_dict["end"] = int(out_dict["end"])


    if attr_dict["strand"] is "+":
        region_str=str(out_dict["chr"])+":"+str(out_dict["start"])+"-"+str(out_dict["end"])
        match_seq = pysam_fasta.fetch(region=region_str)
        match_seq = Seq(match_seq,IUPAC.unambiguous_dna)
    else:
        region_str=str(out_dict["chr"])+":"+str(out_dict["end"])+"-"+str(out_dict["start"])
        match_seq = pysam_fasta.fetch(region=region_str)
        match_seq = Seq(match_seq,IUPAC.unambiguous_dna).reverse_complement()

    masked_seq = getMaskedSeq(primer_seq,match_seq)
    out_dict["masked_seq"] = str(masked_seq)
    out_dict["match_seq"] = str(match_seq)
    return(out_dict)
