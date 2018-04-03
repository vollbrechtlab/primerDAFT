#pysam
import pysam

#BioPython
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Pandas
import pandas as pd
from pandas.compat import StringIO

#primerDAFT
from primerDAFT.specCheck.getMaskedSeq import getMaskedSeq

def getOfftargetAttrs(off_target,side,idx,data,side_cols,target_cols,pysam_fasta):
    primer_seq = data["result"]["pairs"][idx]["PRIMER_"+side.upper()]["SEQUENCE"]
    attr_names=["sseqid","sstart","send","strand"]
    json_names=["chr","start","end","strand"]
    cols = [x+"_"+side for x in side_cols] + target_cols
    attrs = [getattr(off_target,x) for x in cols]
    attr_dict=dict(zip(side_cols+target_cols,attrs))
    out_dict=None

    if attr_dict["qlen"]>attr_dict["length"]:
        if attr_dict["strand"] is "+":
            attr_dict["sstart"] = attr_dict["sstart"]- (attr_dict["qstart"]-1)
            attr_dict["send"] =   attr_dict["send"]  + (attr_dict["qlen"]-attr_dict["qend"])
        else:
            attr_dict["sstart"] = attr_dict["sstart"]+ (attr_dict["qstart"]-1)
            attr_dict["send"] =   attr_dict["send"]  - (attr_dict["qlen"]-attr_dict["qend"])
        out_dict = dict(zip(json_names,[attr_dict.get(x) for x in attr_names]))
    else:
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
