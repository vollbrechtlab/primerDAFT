from primer3.bindings import designPrimers
from primerDAFT.designPrimers.transformInput import transformInput
from primerDAFT.designPrimers.createBetterResult import createBetterResult

def findPrimers(inputData, resultFormat="better"):
    """ return primer3 result with given format

    Parameters
    ----------
    inputData : input data
        asd
    resultFormat : result format (raw/better),optional
        sadsa

    Returns
    ----------
    result : sad
    """
    inputData = transformInput(inputData)

    result = {}
    try:
        result = designPrimers(inputData['seq_args'], inputData['global_args'])
    except: # input data is broken
        raise Exception('input data is broken')

    if resultFormat == "better":
        return createBetterResult(result);

    return result
