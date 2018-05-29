def createBetterResult(result):
    """ Create a primer3 result in a better format

    Parameters
    ----------
    result : dict
        primer3 result.

    Returns
    -------
    betterResult : dict
        primer3 better results

    """

    betterResult = {}
    betterResult['pairs'] = []

    for i in range(max([result['PRIMER_PAIR_NUM_RETURNED'],result['PRIMER_LEFT_NUM_RETURNED'],result['PRIMER_RIGHT_NUM_RETURNED'],result['PRIMER_INTERNAL_NUM_RETURNED']])):
        betterResult['pairs'].append({});

    for i in range(result['PRIMER_PAIR_NUM_RETURNED']):
        betterResult['pairs'][i]['COMPL_ANY_TH'] = result['PRIMER_PAIR_{}_COMPL_ANY_TH'.format(i)]
        betterResult['pairs'][i]['COMPL_END_TH'] = result['PRIMER_PAIR_{}_COMPL_END_TH'.format(i)]
        betterResult['pairs'][i]['PENALTY'] = result['PRIMER_PAIR_{}_PENALTY'.format(i)]
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
