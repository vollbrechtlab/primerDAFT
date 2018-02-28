################### SPEC CHECK ######################

'''
    #Importing modules
'''
#Base Python
import json
import re

#printing
import pprint as pp

#Bio Python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

#pysam
import pysam

#Pandas
import pandas as pd
from pandas.compat import StringIO


'''
Defining functions to be used in the file
'''

# def get_masked_seq(primer_seq,genome_seq):
#     results = map(lambda x,y:"." if x==y else y,primer_seq,genome_seq)
#     out_seq = "".join(list(results))
#     return(out_seq)


# def get_target_attrs(target,side,idx,data,side_cols,target_cols,pysam_fasta):
#     primer_seq = data["result"]["pairs"][idx]["PRIMER_"+side.upper()]["SEQUENCE"]
#     attr_names=["sseqid","sstart","send","strand"]
#     json_names=["chr","start","end","strand"]
#
#     cols = [x+"_"+side for x in side_cols] + target_cols
#     attrs = [getattr(target,x) for x in cols]
#     attr_dict=dict(zip(side_cols + target_cols,attrs))
#
#     out_dict = dict(zip(json_names,[attr_dict.get(x) for x in attr_names]))
#     out_dict["start"] = int(out_dict["start"])
#     out_dict["end"] = int(out_dict["end"])
#
#     if attr_dict["strand"] is "+":
#         match_seq = pysam_fasta.fetch(region=out_dict["chr"]+":"+str(out_dict["start"])+"-"+str(out_dict["end"]))
#         match_seq = Seq(match_seq,IUPAC.unambiguous_dna)
#     else:
#         match_seq = pysam_fasta.fetch(region=out_dict["chr"]+":"+str(out_dict["end"])+"-"+str(out_dict["start"]))
#         match_seq = Seq(match_seq,IUPAC.unambiguous_dna).reverse_complement()
#
#     masked_seq = get_masked_seq(primer_seq,match_seq)
#     out_dict["masked_seq"] = str(masked_seq)
#     out_dict["match_seq"] = str(match_seq)
#     return(out_dict)


# def get_offtarget_attrs(off_target,side,idx,data,side_cols,target_cols,pysam_fasta):
#     primer_seq = data["result"]["pairs"][idx]["PRIMER_"+side.upper()]["SEQUENCE"]
#     attr_names=["sseqid","sstart","send","strand"]
#     json_names=["chr","start","end","strand"]
#     cols = [x+"_"+side for x in side_cols] + target_cols
#     attrs = [getattr(off_target,x) for x in cols]
#     attr_dict=dict(zip(side_cols+target_cols,attrs))
#     out_dict=None
#
#     if attr_dict["qlen"]>attr_dict["length"]:
#         if attr_dict["strand"] is "+":
#             attr_dict["sstart"] = attr_dict["sstart"]- (attr_dict["qstart"]-1)
#             attr_dict["send"] =   attr_dict["send"]  + (attr_dict["qlen"]-attr_dict["qend"])
#         else:
#             attr_dict["sstart"] = attr_dict["sstart"]+ (attr_dict["qstart"]-1)
#             attr_dict["send"] =   attr_dict["send"]  - (attr_dict["qlen"]-attr_dict["qend"])
#         out_dict = dict(zip(json_names,[attr_dict.get(x) for x in attr_names]))
#     else:
#         out_dict = dict(zip(json_names,[attr_dict.get(x) for x in attr_names]))
#
#     out_dict["start"] = int(out_dict["start"])
#     out_dict["end"] = int(out_dict["end"])
#
#     if attr_dict["strand"] is "+":
#         match_seq = pysam_fasta.fetch(region=out_dict["chr"]+":"+str(out_dict["start"])+"-"+str(out_dict["end"]))
#         match_seq = Seq(match_seq,IUPAC.unambiguous_dna)
#     else:
#         match_seq = pysam_fasta.fetch(region=out_dict["chr"]+":"+str(out_dict["end"])+"-"+str(out_dict["start"]))
#         match_seq = Seq(match_seq,IUPAC.unambiguous_dna).reverse_complement()
#
#     masked_seq = get_masked_seq(primer_seq,match_seq)
#     out_dict["masked_seq"] = str(masked_seq)
#     out_dict["match_seq"] = str(match_seq)
#     return(out_dict)


# def specCheck(task, taskResult):
#     '''
#     Dealing with input files both task and task results
#     '''
#
#     with open('config.json') as config_file:
#         config = json.load(config_file)
#
#     task_data = task
#     data = taskResult
#
#     tmp_id = task_data['task_id'] + "_taskResult.json"
#
#
#     '''
#     Defining the basic variables
#     '''
#
#     #blast parameters
#     blast_cols = ["qseqid","qseq","qlen","qstart","qend","sseqid","sseq","sstart","send","length","mismatch","evalue"]
#     out_fmt="'6 "+" ".join(blast_cols)+"'"
#     num_primer_pairs = len(data["result"]["pairs"])
#     genome_fasta = config[task_data["spec_check"]["GENOME"]]
#     pysam_fasta = pysam.FastaFile(genome_fasta)
#
#     blastPath = os.path.dirname(__file__) + "/blastn"
#     if "blastn_path" in config and config["blastn_path"] is not None:
#         blastPath = config["blastn_path"]
#
#     print("using " + blastPath)
#
#
#     '''
#     input_seqs = []
#     input_seq = Seq(task_data['input_data']["SEQUENCE_TEMPLATE"],IUPAC.unambiguous_dna)
#     input_sr = SeqRecord(input_seq,id="input",description="")
#     input_seqs.append(input_sr)
#     SeqIO.write(input_seqs,"input_seq.fa","fasta")
#     blast_cline = NcbiblastnCommandline(query="input_seq.fa",db=genome_fasta, task="megablast",outfmt=out_fmt,evalue=1)
#     stdout, stderr = blast_cline()
#
#
#     # In[427]:
#
#
#     pd_data = StringIO(stdout)
#     df_input = pd.read_csv(pd_data,sep="\t",header=None,names=blast_cols)
#     df_input.to_csv("input_blast.out",index=None)
#     '''
#
#     # In[430]:
#
#     '''
#     Runnnig the steps
#     '''
#
#     #Saving primers as fasta file
#     primer_fa=tmp_id+".primer.fa"
#     seqs = []
#     counter=0
#     for pair in data["result"]["pairs"]:
#         #pp.pprint(pair)
#         left_seq = Seq(pair['PRIMER_LEFT']['SEQUENCE'],IUPAC.unambiguous_dna)
#         right_seq = Seq(pair['PRIMER_RIGHT']['SEQUENCE'],IUPAC.unambiguous_dna)
#         left_sr = SeqRecord(left_seq,id="left_"+str(counter),description="")
#         right_sr = SeqRecord(right_seq,id="right_"+str(counter),description="")
#         seqs.append(left_sr)
#         seqs.append(right_sr)
#         counter = counter+1
#     SeqIO.write(seqs,primer_fa,"fasta")
#
#
#
#     #Run biopython BLASTN
#     blast_cline = NcbiblastnCommandline(cmd=blastPath,query=primer_fa,db=genome_fasta, task="blastn-short",outfmt=out_fmt,evalue=150,num_threads=4)
#     stdout, stderr = blast_cline()
#
#     pd_data = StringIO(stdout)
#     df = pd.read_csv(pd_data,sep="\t",header=None,names=blast_cols)
#     primer_bl_out = re.sub("fa","bl.out",primer_fa)
#     df.to_csv(primer_bl_out,index=None)
#     df['strand'] = df.apply(lambda x: "+" if x['send'] > x['sstart'] else "-",axis=1)
#
#
#     #get_gen_start = lambda row: row['send'] if row['sstart'] < row['send'] else row['sstart']
#     #df['pstart'] = df.apply(get_gen_start,axis=1)
#
#
#     #Clean the blast out to remove targets with
#
#     df["side"],df["pair"] = zip(*df["qseqid"].map(lambda x:x.split("_")))
#     df["all_mismatch"] = df.mismatch+(df.qlen-df.length)
#     df_filt = df[df.all_mismatch<task_data["spec_check"]["TOTAL_MISMATCH_IGNORE"]]
#
#
#     sel_cols=["side","pair","qseqid","sseqid","qseq","qstart","qend","sstart","send","strand","qlen","length","all_mismatch"]
#     idx_col = ["pair","sseqid"]
#     a = df_filt[(df_filt.all_mismatch<2) & (df_filt.length == df_filt.qlen) & (df_filt.side=="left")][sel_cols]
#     b = df_filt[(df_filt.all_mismatch<2) & (df_filt.length == df_filt.qlen) & (df_filt.side=="right")][sel_cols]
#     c = pd.merge(a,b,on=idx_col,how="outer",suffixes=["_left","_right"])
#     pot_targets = c[(abs(c["sstart_left"]-c["sstart_right"])<task_data["spec_check"]["MAX_TARGET_SIZE"]) & (c.strand_left != c.strand_right)]
#
#
#
#     target_cols = ["sseqid"]
#     side_cols = ["sstart","send","strand","qstart","qend","length","qlen"]
#     for x in list(pd.unique(pot_targets["pair"])):
#         x = int(x)
#         curr_targets = pot_targets[pot_targets.pair==str(x)].reset_index(None)
#         all_pot_targets=[]
#
#         for target in curr_targets.itertuples():
#             target_dict = {"left":{},"right":{},"target_seq":"","prod_size":0}
#             target_dict["left"] = get_target_attrs(target,"left",x,data,side_cols,target_cols,pysam_fasta)
#             target_dict["right"] = get_target_attrs(target,"right",x,data,side_cols,target_cols,pysam_fasta)
#             if target_dict["left"]["strand"] == "+":
#                 gen_start=target_dict["left"]["start"]
#                 gen_end=target_dict["right"]["start"]
#             else:
#                 gen_start=target_dict["right"]["start"]
#                 gen_end=target_dict["left"]["start"]
#             target_dict["prod_size"] = gen_end - gen_start + 1
#             target_dict["target_seq"] = pysam_fasta.fetch(region=target_dict["left"]["chr"]+":"+str(gen_start)+"-"+str(gen_end))
#             all_pot_targets.append(target_dict)
#         data["result"]["pairs"][x]["all_targets"] = {"pot_targets":all_pot_targets}
#
#
#     a = df_filt[((df_filt.mismatch > 0) | (df_filt.qlen > df_filt.length)) & (df_filt.side=="left")][sel_cols]
#     b = df_filt[((df_filt.mismatch > 0) | (df_filt.qlen > df_filt.length)) & (df_filt.side=="right")][sel_cols]
#     c = pd.merge(a,b,on=idx_col,how="outer",suffixes=["_left","_right"])
#     off_targets = c[(abs(c.sstart_left-c.sstart_right)<task_data["spec_check"]["MAX_TARGET_SIZE"]) & (c.strand_left != c.strand_right)]
#
#
#     for x in list(pd.unique(off_targets["pair"])):
#         x = int(x)
#         curr_targets = off_targets[off_targets.pair==str(x)].reset_index(None)
#         all_off_targets=[]
#         for off_target in curr_targets.itertuples():
#             off_target_dict = {"left":{},"right":{},"target_seq":"","prod_size":0}
#             off_target_dict["left"] = get_offtarget_attrs(off_target,"left",x,data,side_cols,target_cols,pysam_fasta)
#             off_target_dict["right"] = get_offtarget_attrs(off_target,"right",x,data,side_cols,target_cols,pysam_fasta)
#
#             if off_target_dict["left"]["strand"] is "+":
#                 gen_start=off_target_dict["left"]["start"]
#                 gen_end=off_target_dict["right"]["start"]
#             else:
#                 gen_start=off_target_dict["right"]["start"]
#                 gen_end=off_target_dict["left"]["start"]
#
#             off_target_dict["prod_size"] = gen_end - gen_start + 1
#             region=off_target_dict["left"]["chr"]+":"+str(gen_start)+"-"+str(gen_end)
#
#             if gen_start < gen_end:
#                 off_target_dict["target_seq"] = pysam_fasta.fetch(region=region)
#                 all_off_targets.append(off_target_dict)
#
#         data["result"]["pairs"][x]["all_targets"]["off_targets"] = all_off_targets
#
#     return data
