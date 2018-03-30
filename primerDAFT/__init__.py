#!/usr/bin/env
# this file makes the package
"""
Created on 29 July 2012
@author: Lisa Simpson
"""

import os
from primerDAFT.designPrimers.findPrimers import findPrimers
from primerDAFT.specCheck.specCheck import specCheck
import json,sys
from pprint import pprint

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

    try: # try to get result
        if 'format' in task:
            result['result'] = findPrimers(task['primer3_data'], task['format'])
        else:
            result['result'] = findPrimers(task['primer3_data'])
    
    except Exception as e:
        result['status'] = 'error'
        result['error_statement'] = 'primer3_data is broken'
    
    else: # no problem
        result['status'] = 'ok'

    if result['status'] == 'ok':
        try:
            result['result'] = specCheck(task, result, configFile)
        except Exception as e:
            print(e)
            result['status'] = 'error'
            result['error_statement'] = 'spec_check is broken'

    return result
