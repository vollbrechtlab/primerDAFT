#!/usr/bin/env

"""
	Author: Takao Shibamoto
	Date: 9/30/2017
	Description: Excute primer3, format output
"""

from primer3.bindings import designPrimers
import json

def printJson(data):
	print(json.dumps(data, sort_keys=True, indent=4))

def transformInput(data):
	""" separate input to seq_args and global_args
	Args: 
		param1: input data
	Returns:
		separated input data
	"""
	p3py_data = {}
	p3py_data['seq_args'] = {}
	p3py_data['global_args'] = {}
	for key in data.keys():
		if('SEQUENCE_' in key.upper()):
			p3py_data['seq_args'][key.upper()] = data[key]
		elif('PRIMER_' in key.upper()):
			p3py_data['global_args'][key.upper()] = data[key]

	return p3py_data


def createBetterResult(result):
	""" Create a primer3 result in a better format
	Args: 
		param1: primer3 result
	Returns:
		better result
	"""

	betterResult = {}
	betterResult['pairs'] = []

	print('pair:', result['PRIMER_PAIR_NUM_RETURNED'])
	print('left:', result['PRIMER_LEFT_NUM_RETURNED'])
	print('right:', result['PRIMER_RIGHT_NUM_RETURNED'])
	print('internal:', result['PRIMER_INTERNAL_NUM_RETURNED'])

	for i in range(max([result['PRIMER_PAIR_NUM_RETURNED'],result['PRIMER_LEFT_NUM_RETURNED'],result['PRIMER_RIGHT_NUM_RETURNED'],result['PRIMER_INTERNAL_NUM_RETURNED']])):
		betterResult['pairs'].append({});

	for i in range(result['PRIMER_PAIR_NUM_RETURNED']):
		betterResult['pairs'][i]['COMPL_ANY_TH'] = result['PRIMER_PAIR_{}_COMPL_ANY_TH'.format(i)]
		betterResult['pairs'][i]['COMPL_END_TH'] = result['PRIMER_PAIR_{}_COMPL_END_TH'.format(i)]
		betterResult['pairs'][i]['PENALLTY'] = result['PRIMER_PAIR_{}_PENALTY'.format(i)]
		betterResult['pairs'][i]['PRODUCT_SIZE'] = result['PRIMER_PAIR_{}_PRODUCT_SIZE'.format(i)]

	for i in range(result['PRIMER_LEFT_NUM_RETURNED']):
		betterResult['pairs'][i]['PRIMER_LEFT'] = {}
		betterResult['pairs'][i]['PRIMER_LEFT']['START'] = result['PRIMER_LEFT_{}'.format(i)][0]
		betterResult['pairs'][i]['PRIMER_LEFT']['LENGTH'] = result['PRIMER_LEFT_{}'.format(i)][1]
		betterResult['pairs'][i]['PRIMER_LEFT']['END_STABILITY'] = result['PRIMER_LEFT_{}_END_STABILITY'.format(i)]
		betterResult['pairs'][i]['PRIMER_LEFT']['GC_PERCENT'] = result['PRIMER_LEFT_{}_GC_PERCENT'.format(i)]
		betterResult['pairs'][i]['PRIMER_LEFT']['HAIRPIN_TH'] = result['PRIMER_LEFT_{}_HAIRPIN_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_LEFT']['PENALTY'] = result['PRIMER_LEFT_{}_PENALTY'.format(i)]
		betterResult['pairs'][i]['PRIMER_LEFT']['SELF_ANY_TH'] = result['PRIMER_LEFT_{}_SELF_ANY_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_LEFT']['SELF_END_TH'] = result['PRIMER_LEFT_{}_SELF_END_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_LEFT']['SEQUENCE'] = result['PRIMER_LEFT_{}_SEQUENCE'.format(i)]
		betterResult['pairs'][i]['PRIMER_LEFT']['TM'] = result['PRIMER_LEFT_{}_TM'.format(i)]

	for i in range(result['PRIMER_RIGHT_NUM_RETURNED']):
		betterResult['pairs'][i]['PRIMER_RIGHT'] = {}
		betterResult['pairs'][i]['PRIMER_RIGHT']['START'] = result['PRIMER_RIGHT_{}'.format(i)][0]
		betterResult['pairs'][i]['PRIMER_RIGHT']['LENGTH'] = result['PRIMER_RIGHT_{}'.format(i)][1]
		betterResult['pairs'][i]['PRIMER_RIGHT']['END_STABILITY'] = result['PRIMER_RIGHT_{}_END_STABILITY'.format(i)]
		betterResult['pairs'][i]['PRIMER_RIGHT']['GC_PERCENT'] = result['PRIMER_RIGHT_{}_GC_PERCENT'.format(i)]
		betterResult['pairs'][i]['PRIMER_RIGHT']['HAIRPIN_TH'] = result['PRIMER_RIGHT_{}_HAIRPIN_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_RIGHT']['PENALTY'] = result['PRIMER_RIGHT_{}_PENALTY'.format(i)]
		betterResult['pairs'][i]['PRIMER_RIGHT']['SELF_ANY_TH'] = result['PRIMER_RIGHT_{}_SELF_ANY_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_RIGHT']['SELF_END_TH'] = result['PRIMER_RIGHT_{}_SELF_END_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_RIGHT']['SEQUENCE'] = result['PRIMER_RIGHT_{}_SEQUENCE'.format(i)]
		betterResult['pairs'][i]['PRIMER_RIGHT']['TM'] = result['PRIMER_RIGHT_{}_TM'.format(i)]

	for i in range(result['PRIMER_INTERNAL_NUM_RETURNED']):
		betterResult['pairs'][i]['PRIMER_INTERNAL'] = {}
		betterResult['pairs'][i]['PRIMER_INTERNAL']['START'] = result['PRIMER_INTERNAL_{}'.format(i)][0]
		betterResult['pairs'][i]['PRIMER_INTERNAL']['LENGTH'] = result['PRIMER_INTERNAL_{}'.format(i)][1]
		betterResult['pairs'][i]['PRIMER_INTERNAL']['GC_PERCENT'] = result['PRIMER_INTERNAL_{}_GC_PERCENT'.format(i)]
		betterResult['pairs'][i]['PRIMER_INTERNAL']['HAIRPIN_TH'] = result['PRIMER_INTERNAL_{}_HAIRPIN_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_INTERNAL']['PENALTY'] = result['PRIMER_INTERNAL_{}_PENALTY'.format(i)]
		betterResult['pairs'][i]['PRIMER_INTERNAL']['SELF_ANY_TH'] = result['PRIMER_INTERNAL_{}_SELF_ANY_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_INTERNAL']['SELF_END_TH'] = result['PRIMER_INTERNAL_{}_SELF_END_TH'.format(i)]
		betterResult['pairs'][i]['PRIMER_INTERNAL']['SEQUENCE'] = result['PRIMER_INTERNAL_{}_SEQUENCE'.format(i)]
		betterResult['pairs'][i]['PRIMER_INTERNAL']['TM'] = result['PRIMER_INTERNAL_{}_TM'.format(i)]
	
	return betterResult


def findPrimers(inputData, resultFormat="raw"):
	""" return primer3 result with given format
	Args: 
		param1: input data
		param2: result format (raw/better)
	Returns:
		result
	"""
	p3pyInputData = transformInput(inputData)

	result = {}
	try:
		result = designPrimers(p3pyInputData['seq_args'], p3pyInputData['global_args'])
	except: # input data is broken
		raise Exception('input data is broken')

	if resultFormat == "better":
		return createBetterResult(result);

	return result


def findPrimersFile(taskId):
	""" Create a primer3 result in a better format
	Args: 
		param1: task ID
	"""

	taskPath = 'cache/'+taskId+'_task.json'
	taskResultPath = 'cache/'+taskId+'_taskResult.json'

	taskFile = None
	try: 
		taskFile = open(taskPath, 'r')
	except IOError as e:
		raise Exception("input file does not exist")

	task = None
	try: # try to load input json file
		task = json.load(taskFile)
	except:
		raise Exception('input data is broken')

	taskResult = {}
	try: # try to get result
		if 'format' in task:
			taskResult['result'] = findPrimers(task['input_data'], task['format'])
		else:
			taskResult['result'] = findPrimers(task['input_data'])
	
	except Exception as e:
		taskResult['status'] = 'error'
		taskResult['error_statement'] = str(e)
	
	else: # no problem
		taskResult['status'] = 'ok'

	# save the result to a file
	with open(taskResultPath, 'w') as newTaskResultFile:
		json.dump(taskResult, newTaskResultFile, sort_keys = True, indent = 4, ensure_ascii = False)

	
def saveTask(newTask, taskId):
	""" Save a task to a file
	"""

	newTaskPath = 'cache/'+taskId+'_task.json'

	# save input data as json file
	with open(newTaskPath, 'w') as newTaskFile:
		json.dump(newTask, newTaskFile)

	print("New task (ID:{}) added".format(taskId))


def openTaskResultFile(taskId):
	return open('cache/'+taskId+'_taskResult.json', 'r') 