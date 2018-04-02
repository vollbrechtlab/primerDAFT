#!/usr/bin/env

import primerDAFT
from pprint import pprint
import json

def test1_specCheck():
    result = primerDAFT.specCheck("test_task.json", "test_result.json")
    pprint(result)

def test_findPrimers():
    taskFilePath = "test_task.json"
    with open(taskFilePath) as taskFile:
        taskData = json.load(taskFile)
    result = primerDAFT.findPrimers(taskData['input_data'])
    pprint(result)

def test_run():
    taskFilePath = "test_data/test_task_maize.json"
    # taskFilePath = "test_data/test_task.json"
    with open(taskFilePath) as taskFile:
        taskData = json.load(taskFile)
    with open('cache/'+taskData['taskId']+'_task.json', 'w') as outfile:
        json.dump(taskData, outfile)
    configFile = "primer-dx.conf"
    result = primerDAFT.run(taskData, configFile)
    # pprint(result)
    with open('cache/'+taskData['taskId']+'_result.json', 'w') as outfile:
        json.dump(result, outfile,indent=4)

test_run()
