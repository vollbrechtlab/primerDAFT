#!/usr/bin/env

import primerDAFT
from pprint import pprint
import json, os, csv

task_f="erica_task.json"
result_f="erica_result.json"

def test_run():
    with open(task_f, 'r') as f:
        taskData = json.load(f)
        #taskData['spec_check_data'] = None

    result = primerDAFT.run(taskData, "pdaft.conf")
    with open(taskData['taskId']+'_result.json', 'w') as f:
        json.dump(result, f, indent=4)

def test_createCSV():
    with open(result_f, 'r') as f:
        result = json.load(f)

    resultCSV = primerDAFT.createCSV(result)
    with open(result['taskId']+'_result.csv', 'w') as f:
        f.write(resultCSV)

test_run()
test_createCSV()
