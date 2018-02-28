#!/usr/bin/env
# this file makes the package
"""
Created on 29 July 2012
@author: Lisa Simpson
"""

import os
from primer_blast_dx.designPrimers.findPrimersFromTask import findPrimersFromTask
from primer_blast_dx.specCheck.specCheck import specCheck
import json,sys
import pprint as pp

def run(task,configFile):
    """ run findPrimers and specCheck

    Parameters
    ----------
        task (dic): task data
        saveTmp (bool): true if user wants to save temporary files

    Returns:
        dict: result dictionary

    """

    taskResult = findPrimersFromTask(task)

    taskResultPath = task['task_id'] + "_taskResult.json"
     # save the primer result to a file
    with open(taskResultPath, 'w') as newTaskResultFile:
        json.dump(taskResult, newTaskResultFile, sort_keys = True, indent = 4, ensure_ascii = False)

    specResult = specCheck(task, taskResult,configFile)
    sys.exit()
    return specResult
