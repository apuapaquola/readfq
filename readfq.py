def readfq(fp):
    """generator function that parses fastq or fasta files, yielding
       (name, sequence, qual) for fastq files, or (name, sequence, None) 
       for fasta files, for each record
    """

    last = None # this is a buffer keeping the last unprocessed line

    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                l = l.rstrip()
                if l[0] in '>@': # fasta/q header line
                    last = l     # save this line
                    break
                    
        if not last: break

        name = last[1:]
        seqs = []
        last = None
        for l in fp: # read the sequence
            l = l.rstrip()
            if l[0] in '@+>':
                last = l
                break
            seqs.append(l)

        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq = ''.join(seqs)
            leng = 0
            seqs = []
            for l in fp: # read the quality
                l = l.rstrip()
                seqs.append(l)
                leng += len(l)
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

if __name__ == "__main__":
    import sys
    n, slen, qlen = 0, 0, 0
    for name, seq, qual in readfq(sys.stdin):
        n += 1
        slen += len(seq)
        qlen += qual and len(qual) or 0
    print(n, '\t', slen, '\t', qlen)
