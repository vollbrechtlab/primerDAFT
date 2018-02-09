#!/usr/bin/env

import primer_blast_dx

#primer_blast_dx.findPrimersFile("0Y7Cnlt4E37SP2W9_task.json", "someresult.json")

result = primer_blast_dx.specCheck("0Y7Cnlt4E37SP2W9_task.json", "someresult.json")
print(result)

"""
out_json=open(tmp_id+".out.json","w")
#pp.pprint(data)
json.dump(data,out_json,indent="\t",sort_keys=True)
out_json.close()
"""