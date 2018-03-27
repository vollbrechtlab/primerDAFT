def transformInput(data):
    """ separate input to seq_args and global_args
    
    Args:
        data: input data
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
