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
    taskFilePath = "test_task.json"
    with open(taskFilePath) as taskFile:
        taskData = json.load(taskFile)
    configFile = "primer-dx.conf"
    result = primerDAFT.run(taskData, configFile)
    #pprint(result)

test_run()
