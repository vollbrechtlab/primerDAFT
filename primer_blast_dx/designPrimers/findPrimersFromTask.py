from primer_blast_dx.designPrimers.findPrimers import findPrimers

def findPrimersFromTask(task):
    """ return primer3 result from task. Checks exception
    
    Args:
        inputData: input data
    Returns:
        task result
    """

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

    return taskResult
