#!/usr/bin/env

__all__ = ['designPrimers', 'specCheck']

import os, json, sys

from primerDAFT.designPrimers.findPrimers import findPrimers
from primerDAFT.specCheck.specCheck import specCheck

def run(task, configFile):
    """ run findPrimers and specCheck

    Parameters
    ----------
        task (dic): task data
        configFile (string): config file path

    Returns:
        dict: result dictionary

    """

    result = {}
    result['taskId'] = task['taskId']

    try: # try to get primer3 result
        if 'format' in task:
            result['result'] = findPrimers(task['primer3_data'], task['format'])
        else:
            result['result'] = findPrimers(task['primer3_data'])
    except Exception as e:
        result['status'] = 'error'
        result['error_statement'] = 'primer3_data is broken'
        return result

    result['status'] = 'primer3 ok'

    if ('spec_check_data' in task) and (task['spec_check_data'] != None) and (task['spec_check_data']['RUN_SPEC_CHECK'] == 1):
        try:
            specCheck(task, result, configFile)
        except Exception as e:
            #print(e)
            result['status'] += ' speck check error'
            result['error_statement'] = 'no match for this genome'
        else:
            result['status'] = 'all ok'

    return result


def createCSV(result):
    if 'ok' in result['status']:
        resultCSV = 'primer pair,product size,penalty,any th,end th'

        leftExists = False
        rightExists = False
        internalExists = False
        if 'PRIMER_LEFT' in result['result']['pairs'][0]:
            leftExists = True
            resultCSV += ',left_start,left_length,left_tm,left_gc,left_seq,left_hairpin_th,left_penalty,left_self_any_th,left_self_end_th,left_end_stability'
        if 'PRIMER_RIGHT' in result['result']['pairs'][0]:
            rightExists = True
            resultCSV += ',right_start,right_length,right_tm,right_gc,right_seq,right_hairpin_th,right_penalty,right_self_any_th,right_self_end_th,right_end_stability'
        if 'PRIMER_INTERNAL' in result['result']['pairs'][0]:
            internalExists = True
            resultCSV += ',internal_start,internal_length,internal_tm,internal_gc,internal_seq,internal_hairpin_th,internal_penalty,internal_self_any_th,internal_self_end_th'
        
        if 'all ok' in result['status']:
            resultCSV += ',targets,off targets'
        resultCSV += '\n'

        i = 1
        for pair in result['result']['pairs']:
            resultCSV += "{},{},{},{},{}".format(i, pair['PRODUCT_SIZE'], pair['PENALTY'], pair['COMPL_ANY_TH'], pair['COMPL_END_TH'])

            if leftExists:
                resultCSV += ",{},{},{},{},{},{},{},{},{},{}".format(pair['PRIMER_LEFT']['START'], pair['PRIMER_LEFT']['LENGTH'] , pair['PRIMER_LEFT']['TM'], pair['PRIMER_LEFT']['GC_PERCENT'], pair['PRIMER_LEFT']['SEQUENCE'], pair['PRIMER_LEFT']['HAIRPIN_TH'], pair['PRIMER_LEFT']['PENALTY'], pair['PRIMER_LEFT']['SELF_ANY_TH'], pair['PRIMER_LEFT']['SELF_END_TH'], pair['PRIMER_LEFT']['END_STABILITY'])

            if rightExists:
                resultCSV += ",{},{},{},{},{},{},{},{},{},{}".format(pair['PRIMER_LEFT']['START'], pair['PRIMER_LEFT']['LENGTH'] , pair['PRIMER_LEFT']['TM'], pair['PRIMER_LEFT']['GC_PERCENT'], pair['PRIMER_LEFT']['SEQUENCE'], pair['PRIMER_LEFT']['HAIRPIN_TH'], pair['PRIMER_LEFT']['PENALTY'], pair['PRIMER_LEFT']['SELF_ANY_TH'], pair['PRIMER_LEFT']['SELF_END_TH'], pair['PRIMER_LEFT']['END_STABILITY'])

            if internalExists:
                resultCSV += ",{},{},{},{},{},{},{},{},{}".format(pair['PRIMER_LEFT']['START'], pair['PRIMER_LEFT']['LENGTH'] , pair['PRIMER_LEFT']['TM'], pair['PRIMER_LEFT']['GC_PERCENT'], pair['PRIMER_LEFT']['SEQUENCE'], pair['PRIMER_LEFT']['HAIRPIN_TH'], pair['PRIMER_LEFT']['PENALTY'], pair['PRIMER_LEFT']['SELF_ANY_TH'], pair['PRIMER_LEFT']['SELF_END_TH'])
       
            if 'all ok' in result['status']:
                resultCSV += ",{},{}".format(pair['targets'],pair['off_targets'])
            resultCSV += '\n'
            i += 1

        return resultCSV


    return ''
