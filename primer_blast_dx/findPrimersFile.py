def findPrimersFile(taskPath, taskResultPath):
    """ find primers from the task file and store the result to a task result file
    Args:
        taskPath: location of the task file
        taskResultPath: location of the task result file to store
    """

    #taskPath = 'cache/'+taskId+'_task.json'
    #taskResultPath = 'cache/'+taskId+'_taskResult.json'
    #taskResultPath = os.path.basename(taskPath)
    #taskResultPath = "maybe_result.json"
    #taskResultPath = os.path.basename(taskPath) + ".result.json"

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
