from primer3.bindings import designPrimers

def findPrimers(inputData, resultFormat="better"):
    """ return primer3 result with given format
    Args:
        inputData: input data
        resultFormat: result format (raw/better)
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
