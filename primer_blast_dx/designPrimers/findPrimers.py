from primer3.bindings import designPrimers
from primer_blast_dx.designPrimers.transformInput import transformInput
from primer_blast_dx.designPrimers.createBetterResult import createBetterResult

def findPrimers(inputData, resultFormat="better"):
    """ return primer3 result with given format

    Parameters
    ----------
    inputData: input data
        asd
    resultFormat: result format (raw/better),optional
        sadsa

    Returns
    ----------
    result:
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
