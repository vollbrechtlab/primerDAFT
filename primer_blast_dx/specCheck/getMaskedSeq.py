def getMaskedSeq(primer_seq,genome_seq):
    results = map(lambda x,y:"." if x==y else y,primer_seq,genome_seq)
    out_seq = "".join(list(results))
    return(out_seq)
