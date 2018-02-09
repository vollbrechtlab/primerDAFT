#!/usr/bin/env

import primer_blast_dx
from pprint import pprint
import json

#primer_blast_dx.findPrimersFile("0Y7Cnlt4E37SP2W9_task.json", "someresult.json")

def test_specCheck():
    result = primer_blast_dx.specCheck("0Y7Cnlt4E37SP2W9_task.json", "someresult.json")
    print(result)

def test_findPrimers():

    taskFilePath = "0Y7Cnlt4E37SP2W9_task.json"
    with open(taskFilePath) as taskFile:
        taskData = json.load(taskFile)

    result = primer_blast_dx.findPrimers(taskData['input_data'])
    print(result)

def test_run():
    taskFilePath = "0Y7Cnlt4E37SP2W9_task.json"
    with open(taskFilePath) as taskFile:
        taskData = json.load(taskFile)

    result = primer_blast_dx.run(taskData)
    print(result)

test_run()


"""
out_json=open(tmp_id+".out.json","w")
#pp.pprint(data)
json.dump(data,out_json,indent="\t",sort_keys=True)
out_json.close()
"""
