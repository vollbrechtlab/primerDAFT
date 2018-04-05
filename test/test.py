#!/usr/bin/env

import primerDAFT
from pprint import pprint
import json, os

# create cache folder if it doesnt exist yet
if not os.path.exists("cache"):
    os.makedirs("cache")

def test_run():
    taskFilePath = "test_data/test_task_maize.json"
    with open(taskFilePath) as taskFile:
        taskData = json.load(taskFile)
    with open('cache/'+taskData['taskId']+'_task.json', 'w') as outfile:
        json.dump(taskData, outfile, indent=4)

    result = primerDAFT.run(taskData, "pdaft.conf")
    with open('cache/'+taskData['taskId']+'_result.json', 'w') as outfile:
        json.dump(result, outfile, indent=4)

test_run()
