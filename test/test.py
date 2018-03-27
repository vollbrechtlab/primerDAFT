#!/usr/bin/env

import primerDAFT
from pprint import pprint
import json

#primerDAFT.findPrimersFile("0Y7Cnlt4E37SP2W9_task.json", "someresult.json")

def test1_specCheck():
    result = primerDAFT.specCheck("0Y7Cnlt4E37SP2W9_task.json", "someresult.json")
    print(result)

def test2_specCheck():
    result = primerDAFT.specCheck("0Y7Cnlt4E37SP2W9_task.json", "0Y7Cnlt4E37SP2W9_taskResult.json")
    print(result)


def test_findPrimers():

    taskFilePath = "0Y7Cnlt4E37SP2W9_task.json"
    with open(taskFilePath) as taskFile:
        taskData = json.load(taskFile)

    result = primerDAFT.findPrimers(taskData['input_data'])
    print(result)

def test_run():
    taskFilePath = "0Y7Cnlt4E37SP2W9_task.json"
    with open(taskFilePath) as taskFile:
        taskData = json.load(taskFile)
    configFile = "primer-dx.conf"
    result = primerDAFT.run(taskData, configFile)
    pprint(result)

test_run()

#test1_specCheck()


"""
out_json=open(tmp_id+".out.json","w")
#pp.pprint(data)
json.dump(data,out_json,indent="\t",sort_keys=True)
out_json.close()
"""
