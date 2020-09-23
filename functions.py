from globalvars import srange


def namer(count, proteinseq, optionalname=None):
    ll = len(str(len(proteinseq)))
    if optionalname is None:
        primername = 'p_' + '0' * (ll - len(str(count + 1))) \
                     + str(count + 1) + '_' \
                     + proteinseq[count].upper()
    else:
        primername = optionalname
    return primername


def tmcalculator(seq):
    numgc = seq.count('c') + seq.count('g')
    tm = 64.9 + 41 * ((numgc - 16.4) / len(seq))
    tm = round(tm, 2)
    return tm


def checktm(tm_int, tm_target):
    if tm_int < (tm_target - srange):
        return 1
    else:
        if tm_int < (tm_target + srange):
            return 2
        else:
            return 3


def gccontent(seq):
    return 100 * (seq.count('g') + seq.count('c')) / float(len(seq))


def complement(seq):
    compseq = seq.replace('a', '1').replace('g', '2') \
        .replace('t', '3').replace('c', '4') \
        .replace('1', 't').replace('2', 'c') \
        .replace('3', 'a').replace('4', 'g')
    return compseq


def reversecomplement(seq, NNSpresent):
    tmp = []
    if NNSpresent:
        wingseq = seq.split('NNS')
        for wing in wingseq:
            tmp.append(complement(wing))
        rcseq = tmp[1][::-1] + 'SNN' + tmp[0][::-1]
    else:
        cseq = complement(seq)
        rcseq = cseq[::-1]
    return rcseq


def checktrim10nt5p(seq, p_start):
    for m in range(10):
        if seq[0] == 'c' or seq[0] == 'g':
            break
        seq = seq[1:]
        p_start += 1
    return seq, p_start


def checktrim10nt3p(seq, p_end):
    for m in range(10):
        if seq[-1] == 'c' or seq[-1] == 'g':
            break
        seq = seq[:-1]
        p_end -= 1
    return seq, p_end


def trimleft(nucseq, seq, p_start, p_end, tm_target):
    # This function starts trimming the sequence from
    # left then right consequtively...
    k = 0
    while k < 10 and checktm(tmcalculator(seq), tm_target) != 2:
        if k % 2 == 0:
            seq = seq[1:]
            p_start += 1
        else:
            seq = seq[:-1]
            p_end -= 1
        seq, p_start = checktrim10nt5p(seq, p_start)
        seq, p_end = checktrim10nt3p(seq, p_end)
        k += 1
    return seq


def trimright(nucseq, seq, p_start, p_end, tm_target):
    # This function starts trimming the sequence from right
    # then left consequtively...
    k = 0
    while k < 10 and checktm(tmcalculator(seq), tm_target) != 2:
        if k % 2 == 0:
            seq = seq[:-1]
            p_end -= 1
        else:
            seq = seq[1:]
            p_start += 1
        seq, p_start = checktrim10nt5p(seq, p_start)
        seq, p_end = checktrim10nt3p(seq, p_end)
        # print 'TRIMRIGHT-Trial',k,'SEQUENCE',seq
        k += 1
    return seq


def trimboth(nucseq, seq, p_start, p_end, tm_target):
    k = 0
    while k < 10 and checktm(tmcalculator(seq), tm_target) != 2:
        seq = seq[1:-1]
        p_start += 1
        p_end -= 1
        seq, p_start = checktrim10nt5p(seq, p_start)
        seq, p_end = checktrim10nt3p(seq, p_end)
        k += 1
    return seq


def checkadd10nt5p(nucseq, seq, p_start):
    for m in range(10):
        if seq[0] == 'c' or seq[0] == 'g':
            break
        seq = nucseq[p_start - 1] + seq
        p_start -= 1
    return seq, p_start


def checkadd10nt3p(nucseq, seq, p_end):
    for m in range(10):
        if seq[-1] == 'c' or seq[-1] == 'g':
            break
        seq = seq + nucseq[p_end + 1]
        p_end += 1
    return seq, p_end


def addleft(nucseq, seq, p_start, p_end, tm_target, kup=10):
    k = 0
    while k < kup and checktm(tmcalculator(seq), tm_target) != 2:
        if len(seq) >= 39:
            break
        if k % 2 == 0:
            seq = nucseq[p_start - 1] + seq
            p_start -= 1
        else:
            seq = seq + nucseq[p_end + 1]
            p_end += 1
        seq, p_start = checkadd10nt5p(nucseq, seq, p_start)
        seq, p_end = checkadd10nt3p(nucseq, seq, p_end)
        k += 1
        if len(seq) >= 39:
            break
    return seq


def addright(nucseq, seq, p_start, p_end, tm_target, kup=10):
    k = 0
    while k < kup and checktm(tmcalculator(seq), tm_target) != 2:
        if len(seq) >= 39:
            break
        if k % 2 == 0:
            seq = seq + nucseq[p_end + 1]
            p_end += 1

        else:
            seq = nucseq[p_start - 1] + seq
            p_start -= 1
        seq, p_start = checkadd10nt5p(nucseq, seq, p_start)
        seq, p_end = checkadd10nt3p(nucseq, seq, p_end)
        k += 1
        if len(seq) >= 39:
            break
    return seq


def addboth(nucseq, seq, p_start, p_end, tm_target, kup=10):
    k = 0
    while k < kup and checktm(tmcalculator(seq), tm_target) != 2:
        if len(seq) >= 39:
            break
        seq = nucseq[p_start - 1] + seq + nucseq[p_end + 1]
        p_start -= 1
        p_end += 1
        seq, p_start = checkadd10nt5p(nucseq, seq, p_start)
        seq, p_end = checkadd10nt3p(nucseq, seq, p_end)
        k += 1
        if len(seq) >= 39:
            break
    return seq


################################################################
def logic(nucseq, seq, p_start, p_end, tm_target, k=10):
    if tmcalculator(seq) < (tm_target - srange):
        seqleft = addleft(nucseq, seq, p_start, p_end, tm_target, k)
        seqright = addright(nucseq, seq, p_start, p_end, tm_target, k)
        seqboth = addboth(nucseq, seq, p_start, p_end, tm_target, k)
        return seqboth, seqleft, seqright
    else:
        if tmcalculator(seq) < (tm_target + srange):
            seqleft = None
            seqright = None
            return seq, seqleft, seqright
        else:
            seqleft = trimleft(nucseq, seq, p_start, p_end, tm_target)
            seqright = trimright(nucseq, seq, p_start, p_end, tm_target)
            seqboth = trimboth(nucseq, seq, p_start, p_end, tm_target)
            return seqboth, seqleft, seqright


################################################################

def justtrim(nucsequence, seq, p_start, p_end, tm_target):
    for m in range(10):
        if seq[0] == 'c' or seq[0] == 'g':
            break
        seq = seq[1:]
        p_start += 1
    for m in range(10):
        if seq[-1] == 'c' or seq[-1] == 'g':
            break
        seq = seq[:-1]
        p_end -= 1
    seq, seqleft, seqright = logic(nucsequence, seq, p_start, p_end, tm_target)
    return seq, seqleft, seqright


def score(seq, tm_target, flag=0):
    tm = tmcalculator(seq)
    tm_lowlim = tm_target - srange
    tm_uplim = tm_target + srange
    gc_content = ((seq.count('g') + seq.count('c')) / float(len(seq) - 3))
    gc_score = (100 * gc_content * (1. - gc_content))
    tm_score = -((tm - tm_lowlim) * (tm - tm_uplim)) / ((srange ** 2) / 25.)
    len_score =-(len(seq) ** 2 - 70 * len(seq) + 1200)
    # This scoring is optimized for primers with length in between 30-40 nucleotides.Score values: [-Inf to 25]
    if flag == 0:
        wing = (len(seq.split('NNS')[0]) - float(len(seq.split('NNS')[1])))
        wing_score = (25 - wing ** 2)
        totalscore = gc_score + tm_score + len_score + wing_score
        return gc_score, tm_score, len_score, wing_score, totalscore
    else:
        # Since common primers do not have NNS, they also do not have a wing score.
        # For that reason, total score is normalized to 100.
        totalscore = (float(gc_score) + float(tm_score) + float(len_score)) * 4. / 3.
        return gc_score, tm_score, len_score, totalscore
