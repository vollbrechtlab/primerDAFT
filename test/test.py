#!/usr/bin/env

import primerDAFT
from pprint import pprint
import json, os, csv

def test_run():
    with open("C8nSJqbgGqrv5GUT_task.json", 'r') as f:
        taskData = json.load(f)
        #taskData['spec_check_data'] = None

    result = primerDAFT.run(taskData, "pdaft.conf")
    with open(taskData['taskId']+'_result.json', 'w') as f:
        json.dump(result, f, indent=4)

def test_createCSV():
    with open("C8nSJqbgGqrv5GUT_result.json", 'r') as f:
        result = json.load(f)
    resultCSV = primerDAFT.createCSV(result)
    with open(result['taskId']+'_result.csv', 'w') as f:
        f.write(resultCSV)

test_run()
#test_createCSV()
