#basepackage
import primerDAFT
from primerDAFT.specCheck.getTargetAttrs import getTargetAttrs
from primerDAFT.specCheck.getOfftargetAttrs import getOfftargetAttrs
import traceback

#Base Python
import json, os, re, sys, configparser

#printing
import pprint as pp

#pysam
import pysam

#Pandas
import pandas as pd
from pandas.compat import StringIO

#Bio Python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

class SpecError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def specCheck(task, result, configFile):
    '''
    Dealing with input files both task and task results
    '''

    config = configparser.ConfigParser()
    config.read(configFile)

    # with open(configFile) as config_file:
    #     config = json.load(config_file)

    task_data = task
    data = result

    tmp_id = task_data['taskId'] + "_tmp"


    '''
    Defining the basic variables
    '''

    #blast parameters
    blast_cols = config["BLAST"]["blast_cols"].split(",")
    # ["qseqid","qseq","qlen","qstart","qend","sseqid","sseq","sstart","send","length","mismatch","evalue"]
    out_fmt="'6 "+" ".join(blast_cols)+"'"
    num_primer_pairs = len(data["result"]["pairs"])
    genome_fasta = config["GENOMES"][task_data["spec_check_data"]["GENOME"]]
    if os.path.isfile(genome_fasta):
        pysam_fasta = pysam.FastaFile(genome_fasta)
    elif os.path.isfile(genome_fasta+".gz"):
        pysam_fasta = pysam.FastaFile(genome_fasta+".gz")
    else:
        raise SpecError("Genome not Found")
    blastPath = os.path.dirname(os.path.dirname(__file__))+"/bin/blastn"


    if "blastn_path" in config and config["blastn_path"] is not None:
        blastPath = config["blastn_path"]
    print("using " + blastPath)

    '''
    Runnnig the steps
    '''

    #Saving primers as fasta file
    primer_fa=config["GENERAL"]["cache_directory"]+tmp_id+".primer.fa"
    seqs = []
    counter=0
    for pair in data["result"]["pairs"]:
        #pp.pprint(pair)
        left_seq = Seq(pair['PRIMER_LEFT']['SEQUENCE'],IUPAC.unambiguous_dna)
        right_seq = Seq(pair['PRIMER_RIGHT']['SEQUENCE'],IUPAC.unambiguous_dna)
        left_sr = SeqRecord(left_seq,id="left_"+str(counter),description="")
        right_sr = SeqRecord(right_seq,id="right_"+str(counter),description="")
        seqs.append(left_sr)
        seqs.append(right_sr)
        counter = counter+1
    SeqIO.write(seqs,primer_fa,"fasta")


    #Run biopython BLASTN
    blast_cline = NcbiblastnCommandline(cmd=blastPath,query=primer_fa,db=genome_fasta, task="blastn-short",outfmt=out_fmt,evalue=config["BLAST"]["evalue"],num_threads=config["BLAST"]["num_threads"])
    stdout, stderr = blast_cline()

    pd_data = StringIO(stdout)
    df = pd.read_csv(pd_data,sep="\t",header=None,names=blast_cols)
    primer_bl_out = re.sub("fa","bl.out",primer_fa)
    df.to_csv(primer_bl_out,index=None)
    df['strand'] = df.apply(lambda x: "+" if x['send'] > x['sstart'] else "-",axis=1)

    #get_gen_start = lambda row: row['send'] if row['sstart'] < row['send'] else row['sstart']
    #df['pstart'] = df.apply(get_gen_start,axis=1)


    #Clean the blast out to remove targets with

    df["side"],df["pair"] = zip(*df["qseqid"].map(lambda x:x.split("_")))
    df["all_mismatch"] = df.mismatch+(df.qlen-df.length)
    df_filt = df[df.all_mismatch<task_data["spec_check_data"]["TOTAL_MISMATCH_IGNORE"]]



    sel_cols=["side","pair","qseqid","sseqid","qseq","qstart","qend","sstart","send","strand","qlen","length","all_mismatch"]
    idx_col = ["pair","sseqid"]
    '''
        Get potential targets
    '''
    try:
        pot_target_data = get_pot_targets(df_filt,data,task_data,sel_cols,idx_col,pysam_fasta)
    except Exception as e:
        #print(e)
        #traceback.print_exc()
        raise ValueError
    else:
        data = pot_target_data

    try:
        off_target_data = get_off_targets(df_filt,data,task_data,sel_cols,idx_col,pysam_fasta)
    except Exception as e:
        #print(e)
        #traceback.print_exc()
        raise ValueError
    else:
        data = off_target_data

    data = rank_primers(data)

    return data


def get_off_targets(df_filt,data,task_data,sel_cols,idx_col,pysam_fasta):
    max_target_size = int(task_data["spec_check_data"]["MAX_TARGET_SIZE"])
    a = df_filt[((df_filt.all_mismatch > 0) & (df_filt.all_mismatch <= task_data["spec_check_data"]["TOTAL_SPECIFICITY_MISMATCH"])) & (df_filt.side=="left")][sel_cols]
    b = df_filt[((df_filt.all_mismatch > 0) & (df_filt.all_mismatch <= task_data["spec_check_data"]["TOTAL_SPECIFICITY_MISMATCH"])) & (df_filt.side=="right")][sel_cols]
    c = pd.merge(a,b,on=idx_col,how="inner",suffixes=["_left","_right"])
    c["size"] = c["sstart_left"]-c["sstart_right"]
    # off_targets = c[(abs(c.sstart_left-c.sstart_right)<task_data["spec_check_data"]["MAX_TARGET_SIZE"]) & (c.strand_left != c.strand_right)]
    off_targets = c[((c["size"]<=max_target_size) & (c["size"] > 0) & (c["strand_left"] == "-")) | ((c["size"]>=-1*max_target_size) & (c["size"] < 0) & (c["strand_left"] == "+")) & (c["strand_left"] != c["strand_right"])]

    for x in list(pd.unique(off_targets["pair"])):
        x = int(x)
        curr_targets = off_targets[off_targets.pair==str(x)].reset_index(None)
        all_off_targets=[]
        for off_target in curr_targets.itertuples():
            off_target_dict = {"left":{},"right":{},"target_seq":"","prod_size":0}
            off_target_dict["left"] = getOfftargetAttrs(off_target,"left",x,data,side_cols,target_cols,pysam_fasta)
            off_target_dict["right"] = getOfftargetAttrs(off_target,"right",x,data,side_cols,target_cols,pysam_fasta)

            if off_target_dict["left"]["strand"] is "+":
                gen_start=off_target_dict["left"]["start"]
                gen_end=off_target_dict["right"]["start"]
            else:
                gen_start=off_target_dict["right"]["start"]
                gen_end=off_target_dict["left"]["start"]

            off_target_dict["prod_size"] = gen_end - gen_start + 1
            region_str=str(off_target_dict["left"]["chr"])+":"+str(gen_start)+"-"+str(gen_end)

            if gen_start < gen_end:
                off_target_dict["target_seq"] = pysam_fasta.fetch(region=region_str)
                all_off_targets.append(off_target_dict)

        data["result"]["pairs"][x]["all_targets"]["off_targets"] = all_off_targets
    return data

def get_pot_targets(df_filt,data,task_data,sel_cols,idx_col,pysam_fasta):
    max_target_size = int(task_data["spec_check_data"]["MAX_TARGET_SIZE"])
    a = df_filt[(df_filt.all_mismatch==0) & (df_filt.length == df_filt.qlen) & (df_filt.side=="left")][sel_cols]
    b = df_filt[(df_filt.all_mismatch==0) & (df_filt.length == df_filt.qlen) & (df_filt.side=="right")][sel_cols]
    c = pd.merge(a,b,on=idx_col,how="inner",suffixes=["_left","_right"])
    c["size"] = c["sstart_left"]-c["sstart_right"]
    pot_targets = c[((c["size"]<=max_target_size) & (c["size"] > 0) & (c["strand_left"] == "-")) | ((c["size"]>=-1*max_target_size) & (c["size"] < 0) & (c["strand_left"] == "+")) & (c["strand_left"] != c["strand_right"])]


    if len(pot_targets) == 0:
        print("No potential targets in the genome")
        raise ValueError

    target_cols = ["sseqid"]
    side_cols = ["sstart","send","strand","qstart","qend","length","qlen"]

    for x in list(pd.unique(pot_targets["pair"])):
        x = int(x)
        curr_targets = pot_targets[pot_targets.pair==str(x)].reset_index(None)
        all_pot_targets=[]

        for target in curr_targets.itertuples():

            target_dict = {"left":{},"right":{},"target_seq":"","prod_size":0}

            target_dict["left"] = getTargetAttrs(target,"left",x,data,side_cols,target_cols,pysam_fasta)
            target_dict["right"] = getTargetAttrs(target,"right",x,data,side_cols,target_cols,pysam_fasta)

            if target_dict["left"]["strand"] == "+":
                gen_start=target_dict["left"]["start"]
                gen_end=target_dict["right"]["start"]
            else:
                gen_start=target_dict["right"]["start"]
                gen_end=target_dict["left"]["start"]
            target_dict["prod_size"] = gen_end - gen_start + 1

            target_dict["target_seq"] = pysam_fasta.fetch(region=str(target_dict["left"]["chr"])+":"+str(gen_start)+"-"+str(gen_end))
            all_pot_targets.append(target_dict)
        data["result"]["pairs"][x]["all_targets"] = {"pot_targets":all_pot_targets}

    return data

def rank_primers(data):
    idx = pd.Series(range(len(data["result"]["pairs"])))
    targets = pd.Series([len(primer["all_targets"]["pot_targets"]) for primer in data["result"]["pairs"]])
    # targets = pd.Series(np.random.permutation([1,2,4,3,5]))
    off_targets = pd.Series([ len(primer["all_targets"]["off_targets"]) if "off_targets" in primer["all_targets"] else 0 for primer in data["result"]["pairs"]])
    # off_targets = pd.Series([4,3,2,1,1])
    df = pd.DataFrame({
        "idx":idx,
        "targets":targets,
        "off_targets":off_targets
    })
    df_sort = df.sort_values(by=['targets'])
    sort_pairs = []
    sort_pairs = []
    for x in df_sort["idx"]:
        tmp_data = data["result"]["pairs"][x]
        tmp_data["off_targets"] = int(df_sort[df_sort.idx == x].off_targets)
        tmp_data["targets"] = int(df_sort[df_sort.idx == x].targets)
        sort_pairs.append(tmp_data)
    data["result"]["pairs"] = sort_pairs
    return(data)
