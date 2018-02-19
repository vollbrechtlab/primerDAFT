#!/usr/bin/env

# this file makes the package

import os
from primer3.bindings import designPrimers
import json

def transformInput(data):
    """ separate input to seq_args and global_args
    Args:
        data: input data
    Returns:
        separated input data
    """
    p3py_data = {}
    p3py_data['seq_args'] = {}
    p3py_data['global_args'] = {}
    for key in data.keys():
        if('SEQUENCE_' in key.upper()):
            p3py_data['seq_args'][key.upper()] = data[key]
        elif('PRIMER_' in key.upper()):
            p3py_data['global_args'][key.upper()] = data[key]

    return p3py_data


def createBetterResult(result):
    """ Create a primer3 result in a better format
    Args:
        result: primer3 result
    Returns:
        better result
    """

    betterResult = {}
    betterResult['pairs'] = []

    for i in range(max([result['PRIMER_PAIR_NUM_RETURNED'],result['PRIMER_LEFT_NUM_RETURNED'],result['PRIMER_RIGHT_NUM_RETURNED'],result['PRIMER_INTERNAL_NUM_RETURNED']])):
        betterResult['pairs'].append({});

    for i in range(result['PRIMER_PAIR_NUM_RETURNED']):
        betterResult['pairs'][i]['COMPL_ANY_TH'] = result['PRIMER_PAIR_{}_COMPL_ANY_TH'.format(i)]
        betterResult['pairs'][i]['COMPL_END_TH'] = result['PRIMER_PAIR_{}_COMPL_END_TH'.format(i)]
        betterResult['pairs'][i]['PENALLTY'] = result['PRIMER_PAIR_{}_PENALTY'.format(i)]
        betterResult['pairs'][i]['PRODUCT_SIZE'] = result['PRIMER_PAIR_{}_PRODUCT_SIZE'.format(i)]

    for i in range(result['PRIMER_LEFT_NUM_RETURNED']):
        betterResult['pairs'][i]['PRIMER_LEFT'] = {}
        betterResult['pairs'][i]['PRIMER_LEFT']['START'] = result['PRIMER_LEFT_{}'.format(i)][0]
        betterResult['pairs'][i]['PRIMER_LEFT']['LENGTH'] = result['PRIMER_LEFT_{}'.format(i)][1]
        betterResult['pairs'][i]['PRIMER_LEFT']['END_STABILITY'] = result['PRIMER_LEFT_{}_END_STABILITY'.format(i)]
        betterResult['pairs'][i]['PRIMER_LEFT']['GC_PERCENT'] = result['PRIMER_LEFT_{}_GC_PERCENT'.format(i)]
        betterResult['pairs'][i]['PRIMER_LEFT']['HAIRPIN_TH'] = result['PRIMER_LEFT_{}_HAIRPIN_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_LEFT']['PENALTY'] = result['PRIMER_LEFT_{}_PENALTY'.format(i)]
        betterResult['pairs'][i]['PRIMER_LEFT']['SELF_ANY_TH'] = result['PRIMER_LEFT_{}_SELF_ANY_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_LEFT']['SELF_END_TH'] = result['PRIMER_LEFT_{}_SELF_END_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_LEFT']['SEQUENCE'] = result['PRIMER_LEFT_{}_SEQUENCE'.format(i)]
        betterResult['pairs'][i]['PRIMER_LEFT']['TM'] = result['PRIMER_LEFT_{}_TM'.format(i)]

    for i in range(result['PRIMER_RIGHT_NUM_RETURNED']):
        betterResult['pairs'][i]['PRIMER_RIGHT'] = {}
        betterResult['pairs'][i]['PRIMER_RIGHT']['START'] = result['PRIMER_RIGHT_{}'.format(i)][0]
        betterResult['pairs'][i]['PRIMER_RIGHT']['LENGTH'] = result['PRIMER_RIGHT_{}'.format(i)][1]
        betterResult['pairs'][i]['PRIMER_RIGHT']['END_STABILITY'] = result['PRIMER_RIGHT_{}_END_STABILITY'.format(i)]
        betterResult['pairs'][i]['PRIMER_RIGHT']['GC_PERCENT'] = result['PRIMER_RIGHT_{}_GC_PERCENT'.format(i)]
        betterResult['pairs'][i]['PRIMER_RIGHT']['HAIRPIN_TH'] = result['PRIMER_RIGHT_{}_HAIRPIN_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_RIGHT']['PENALTY'] = result['PRIMER_RIGHT_{}_PENALTY'.format(i)]
        betterResult['pairs'][i]['PRIMER_RIGHT']['SELF_ANY_TH'] = result['PRIMER_RIGHT_{}_SELF_ANY_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_RIGHT']['SELF_END_TH'] = result['PRIMER_RIGHT_{}_SELF_END_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_RIGHT']['SEQUENCE'] = result['PRIMER_RIGHT_{}_SEQUENCE'.format(i)]
        betterResult['pairs'][i]['PRIMER_RIGHT']['TM'] = result['PRIMER_RIGHT_{}_TM'.format(i)]

    for i in range(result['PRIMER_INTERNAL_NUM_RETURNED']):
        betterResult['pairs'][i]['PRIMER_INTERNAL'] = {}
        betterResult['pairs'][i]['PRIMER_INTERNAL']['START'] = result['PRIMER_INTERNAL_{}'.format(i)][0]
        betterResult['pairs'][i]['PRIMER_INTERNAL']['LENGTH'] = result['PRIMER_INTERNAL_{}'.format(i)][1]
        betterResult['pairs'][i]['PRIMER_INTERNAL']['GC_PERCENT'] = result['PRIMER_INTERNAL_{}_GC_PERCENT'.format(i)]
        betterResult['pairs'][i]['PRIMER_INTERNAL']['HAIRPIN_TH'] = result['PRIMER_INTERNAL_{}_HAIRPIN_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_INTERNAL']['PENALTY'] = result['PRIMER_INTERNAL_{}_PENALTY'.format(i)]
        betterResult['pairs'][i]['PRIMER_INTERNAL']['SELF_ANY_TH'] = result['PRIMER_INTERNAL_{}_SELF_ANY_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_INTERNAL']['SELF_END_TH'] = result['PRIMER_INTERNAL_{}_SELF_END_TH'.format(i)]
        betterResult['pairs'][i]['PRIMER_INTERNAL']['SEQUENCE'] = result['PRIMER_INTERNAL_{}_SEQUENCE'.format(i)]
        betterResult['pairs'][i]['PRIMER_INTERNAL']['TM'] = result['PRIMER_INTERNAL_{}_TM'.format(i)]

    return betterResult


def findPrimers(inputData, resultFormat="better"):
    """ return primer3 result with given format
    Args:
        inputData: input data
        resultFormat: result format (raw/better)
    Returns:
        result
    """
    p3pyInputData = transformInput(inputData)

    result = {}
    try:
        result = designPrimers(p3pyInputData['seq_args'], p3pyInputData['global_args'])
    except: # input data is broken
        raise Exception('input data is broken')

    if resultFormat == "better":
        return createBetterResult(result);

    return result

def findPrimersFromTask(task):
    """ return primer3 result from task. Checks exception
    Args:
        inputData: input data
    Returns:
        task result
    """

    taskResult = {}
    try: # try to get result
        if 'format' in task:
            taskResult['result'] = findPrimers(task['input_data'], task['format'])
        else:
            taskResult['result'] = findPrimers(task['input_data'])

    except Exception as e:
        taskResult['status'] = 'error'
        taskResult['error_statement'] = str(e)

    else: # no problem
        taskResult['status'] = 'ok'

    return taskResult


def findPrimersFile(taskPath, taskResultPath):
    """ find primers from the task file and store the result to a task result file
    Args:
        taskPath: location of the task file
        taskResultPath: location of the task result file to store
    """

    #taskPath = 'cache/'+taskId+'_task.json'
    #taskResultPath = 'cache/'+taskId+'_taskResult.json'
    #taskResultPath = os.path.basename(taskPath)
    #taskResultPath = "maybe_result.json"
    #taskResultPath = os.path.basename(taskPath) + ".result.json"

    taskFile = None
    try:
        taskFile = open(taskPath, 'r')
    except IOError as e:
        raise Exception("input file does not exist")

    task = None
    try: # try to load input json file
        task = json.load(taskFile)
    except:
        raise Exception('input data is broken')

    taskResult = {}
    try: # try to get result
        if 'format' in task:
            taskResult['result'] = findPrimers(task['input_data'], task['format'])
        else:
            taskResult['result'] = findPrimers(task['input_data'])

    except Exception as e:
        taskResult['status'] = 'error'
        taskResult['error_statement'] = str(e)

    else: # no problem
        taskResult['status'] = 'ok'

    # save the result to a file
    with open(taskResultPath, 'w') as newTaskResultFile:
        json.dump(taskResult, newTaskResultFile, sort_keys = True, indent = 4, ensure_ascii = False)


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

def get_masked_seq(primer_seq,genome_seq):
    results = map(lambda x,y:"." if x==y else y,primer_seq,genome_seq)
    out_seq = "".join(list(results))
    return(out_seq)


def get_target_attrs(target,side,idx,data,side_cols,target_cols,pysam_fasta):
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
        match_seq = pysam_fasta.fetch(region=out_dict["chr"]+":"+str(out_dict["start"])+"-"+str(out_dict["end"]))
        match_seq = Seq(match_seq,IUPAC.unambiguous_dna)
    else:
        match_seq = pysam_fasta.fetch(region=out_dict["chr"]+":"+str(out_dict["end"])+"-"+str(out_dict["start"]))
        match_seq = Seq(match_seq,IUPAC.unambiguous_dna).reverse_complement()

    masked_seq = get_masked_seq(primer_seq,match_seq)
    out_dict["masked_seq"] = str(masked_seq)
    out_dict["match_seq"] = str(match_seq)
    return(out_dict)


def get_offtarget_attrs(off_target,side,idx,data,side_cols,target_cols,pysam_fasta):
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
        match_seq = pysam_fasta.fetch(region=out_dict["chr"]+":"+str(out_dict["start"])+"-"+str(out_dict["end"]))
        match_seq = Seq(match_seq,IUPAC.unambiguous_dna)
    else:
        match_seq = pysam_fasta.fetch(region=out_dict["chr"]+":"+str(out_dict["end"])+"-"+str(out_dict["start"]))
        match_seq = Seq(match_seq,IUPAC.unambiguous_dna).reverse_complement()

    masked_seq = get_masked_seq(primer_seq,match_seq)
    out_dict["masked_seq"] = str(masked_seq)
    out_dict["match_seq"] = str(match_seq)
    return(out_dict)


def specCheck(task, taskResult):
    '''
    Dealing with input files both task and task results
    '''

    with open('genome_config.json') as genome_config_file:
        genome_config = json.load(genome_config_file)

    task_data = task
    data = taskResult

    tmp_id = task_data['task_id'] + "_taskResult.json"


    '''
    Defining the basic variables
    '''

    #blast parameters
    blast_cols = ["qseqid","qseq","qlen","qstart","qend","sseqid","sseq","sstart","send","length","mismatch","evalue"]
    out_fmt="'6 "+" ".join(blast_cols)+"'"
    num_primer_pairs = len(data["result"]["pairs"])
    genome_fasta = genome_config[task_data["spec_check"]["GENOME"]]
    pysam_fasta = pysam.FastaFile(genome_fasta)


    '''
    input_seqs = []
    input_seq = Seq(task_data['input_data']["SEQUENCE_TEMPLATE"],IUPAC.unambiguous_dna)
    input_sr = SeqRecord(input_seq,id="input",description="")
    input_seqs.append(input_sr)
    SeqIO.write(input_seqs,"input_seq.fa","fasta")
    blast_cline = NcbiblastnCommandline(query="input_seq.fa",db=genome_fasta, task="megablast",outfmt=out_fmt,evalue=1)
    stdout, stderr = blast_cline()


    # In[427]:


    pd_data = StringIO(stdout)
    df_input = pd.read_csv(pd_data,sep="\t",header=None,names=blast_cols)
    df_input.to_csv("input_blast.out",index=None)
    '''

    # In[430]:

    '''
    Runnnig the steps
    '''

    #Saving primers as fasta file
    primer_fa=tmp_id+".primer.fa"
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
    blast_cline = NcbiblastnCommandline(query=primer_fa,db=genome_fasta, task="blastn-short",outfmt=out_fmt,evalue=150,num_threads=4)
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
    df_filt = df[df.all_mismatch<task_data["spec_check"]["TOTAL_MISMATCH_IGNORE"]]


    sel_cols=["side","pair","qseqid","sseqid","qseq","qstart","qend","sstart","send","strand","qlen","length","all_mismatch"]
    idx_col = ["pair","sseqid"]
    a = df_filt[(df_filt.all_mismatch<2) & (df_filt.length == df_filt.qlen) & (df_filt.side=="left")][sel_cols]
    b = df_filt[(df_filt.all_mismatch<2) & (df_filt.length == df_filt.qlen) & (df_filt.side=="right")][sel_cols]
    c = pd.merge(a,b,on=idx_col,how="outer",suffixes=["_left","_right"])
    pot_targets = c[(abs(c["sstart_left"]-c["sstart_right"])<task_data["spec_check"]["MAX_TARGET_SIZE"]) & (c.strand_left != c.strand_right)]



    target_cols = ["sseqid"]
    side_cols = ["sstart","send","strand","qstart","qend","length","qlen"]
    for x in list(pd.unique(pot_targets["pair"])):
        x = int(x)
        curr_targets = pot_targets[pot_targets.pair==str(x)].reset_index(None)
        all_pot_targets=[]

        for target in curr_targets.itertuples():
            target_dict = {"left":{},"right":{},"target_seq":"","prod_size":0}
            target_dict["left"] = get_target_attrs(target,"left",x,data,side_cols,target_cols,pysam_fasta)
            target_dict["right"] = get_target_attrs(target,"right",x,data,side_cols,target_cols,pysam_fasta)
            if target_dict["left"]["strand"] == "+":
                gen_start=target_dict["left"]["start"]
                gen_end=target_dict["right"]["start"]
            else:
                gen_start=target_dict["right"]["start"]
                gen_end=target_dict["left"]["start"]
            target_dict["prod_size"] = gen_end - gen_start + 1
            target_dict["target_seq"] = pysam_fasta.fetch(region=target_dict["left"]["chr"]+":"+str(gen_start)+"-"+str(gen_end))
            all_pot_targets.append(target_dict)
        data["result"]["pairs"][x]["all_targets"] = {"pot_targets":all_pot_targets}


    a = df_filt[((df_filt.mismatch > 0) | (df_filt.qlen > df_filt.length)) & (df_filt.side=="left")][sel_cols]
    b = df_filt[((df_filt.mismatch > 0) | (df_filt.qlen > df_filt.length)) & (df_filt.side=="right")][sel_cols]
    c = pd.merge(a,b,on=idx_col,how="outer",suffixes=["_left","_right"])
    off_targets = c[(abs(c.sstart_left-c.sstart_right)<task_data["spec_check"]["MAX_TARGET_SIZE"]) & (c.strand_left != c.strand_right)]


    for x in list(pd.unique(off_targets["pair"])):
        x = int(x)
        curr_targets = off_targets[off_targets.pair==str(x)].reset_index(None)
        all_off_targets=[]
        for off_target in curr_targets.itertuples():
            off_target_dict = {"left":{},"right":{},"target_seq":"","prod_size":0}
            off_target_dict["left"] = get_offtarget_attrs(off_target,"left",x,data,side_cols,target_cols,pysam_fasta)
            off_target_dict["right"] = get_offtarget_attrs(off_target,"right",x,data,side_cols,target_cols,pysam_fasta)

            if off_target_dict["left"]["strand"] is "+":
                gen_start=off_target_dict["left"]["start"]
                gen_end=off_target_dict["right"]["start"]
            else:
                gen_start=off_target_dict["right"]["start"]
                gen_end=off_target_dict["left"]["start"]

            off_target_dict["prod_size"] = gen_end - gen_start + 1
            region=off_target_dict["left"]["chr"]+":"+str(gen_start)+"-"+str(gen_end)

            if gen_start < gen_end:
                off_target_dict["target_seq"] = pysam_fasta.fetch(region=region)
                all_off_targets.append(off_target_dict)

        data["result"]["pairs"][x]["all_targets"]["off_targets"] = all_off_targets

    return data



def run(task):
    """ run findPrimers and specCheck
    Args:
        task (dic): task data
        saveTmp (bool): true if user wants to save temporary files
    Returns:
        dic: result dictionary 
    """

    taskResult = findPrimersFromTask(task)
    taskResultPath = task['task_id'] + "_taskResult.json"
     # save the primer result to a file
    with open(taskResultPath, 'w') as newTaskResultFile:
        json.dump(taskResult, newTaskResultFile, sort_keys = True, indent = 4, ensure_ascii = False)
    
    specResult = specCheck(task, taskResult)

    return specResult
