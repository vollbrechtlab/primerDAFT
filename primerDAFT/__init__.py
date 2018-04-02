#!/usr/bin/env
# this file makes the package
"""
Created on 29 July 2012
@author: Lisa Simpson
"""

__all__ = ['designPrimers', 'specCheck']

import os
from primerDAFT.designPrimers.findPrimers import findPrimers
from primerDAFT.specCheck.specCheck import specCheck
import json,sys

def run(task, configFile):
    """ run findPrimers and specCheck

    Parameters
    ----------
        task (dic): task data
        saveTmp (bool): true if user wants to save temporary files

    Returns:
        dict: result dictionary

    """

    result = {}

    try: # try to get primer3 result
        if 'format' in task:
            result['result'] = findPrimers(task['primer3_data'], task['format'])
        else:
            result['result'] = findPrimers(task['primer3_data'])

    except Exception as e:
        result['status'] = 'error'
        result['error_statement'] = 'primer3_data is broken'

    else: # no problem for primer3
        result['status'] = 'primer3 ok'

        if task["spec_check"]["RUN_SPEC_CHECK"] == 1:
            try:
                specCheck_result = specCheck(task, result, configFile)
            except Exception as e:
                print(e)
                result['status'] += ' speck check error'
                result['error_statement'] = 'no match for this genome'
            else:
                result['status'] = 'all ok'
                result["specCheck_result"] = specCheck_result['result']

    return result
