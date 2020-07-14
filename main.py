import os, time, sys
import pandas as pd
from globalvars import nucsequence, protsequence, \
    tm_target, genename, nucseqcount
from functions import *

print('\n\nPosition example: aaacccATGacacactgtgtct --> POSITION=7')
nucstart = int(input('Provide the position of A in ...ATG...: ')) - 1

df = pd.DataFrame([], columns=['Residue #', 'Primer Name', 'Primer Seq(Fwd)',
                               'Primer Sequence (Rev-Comp)', 'Tm', 'GC Content(%)',
                               'Primer Length', 'GC Score', 'Tm Score', 'Length Score',
                               'Wing Score', 'Total Score'])

primers = {}
primers_left = {}
primers_right = {}

aprimers = {}
aprimers_left = {}
aprimers_right = {}

# Before the CDS starts at least 60 nucleotide needs to be provided.
# This is important for common primer design

# COMMON PRIMER DESIGN
forwcom, forwcom_left, forwcom_right = justtrim(nucsequence, nucsequence[10:50], 10, 50,
                                                tm_target)  # Common forward primer for all reactions
if forwcom_left == None:
    tmfcom = tmcalculator(forwcom)
    primers_right['forward'] = primers_left['forward'] = primers['forward'] = [forwcom, None, tmfcom]
else:
    tmfcom_left = tmcalculator(forwcom_left)
    tmfcom_right = tmcalculator(forwcom_right)
    tmfcom = tmcalculator(forwcom)
    primers['forward'] = [forwcom, None, tmfcom]
    primers_left['forward'] = [forwcom_left, None, tmfcom_left]
    primers_right['forward'] = [forwcom_right, None, tmfcom_right]

revcom, revcom_left, revcom_right = justtrim(nucsequence, nucsequence[-50:-10], len(nucsequence) - 50,
                                             len(nucsequence) - 10,
                                             tm_target)  # Common reverse primer for all reactions

if revcom_left == None:
    revcom = reversecomplement(revcom, 0)
    tmrcom = tmcalculator(revcom)
    primers_right['reverse'] = primers_left['reverse'] = primers['reverse'] = [revcom, None, tmrcom]

else:
    revcom_left = reversecomplement(revcom_left, 0)
    revcom_right = reversecomplement(revcom_right, 0)
    revcom = reversecomplement(revcom, 0)
    tmrcom_left = tmcalculator(revcom_left)
    tmrcom_right = tmcalculator(revcom_right)
    tmrcom = tmcalculator(revcom)

    primers['reverse'] = [revcom, None, tmrcom]
    primers_left['reverse'] = [revcom_left, None, tmrcom_left]
    primers_right['reverse'] = [revcom_right, None, tmrcom_right]

if forwcom_left == None:
    a, b, c, d = score(forwcom, tm_target, 1)
    df.loc[len(df)] = [0, 'Forward Common', str(forwcom), '-', str(tmfcom),
                       gccontent(forwcom), len(forwcom), a, b, c, '-', d]

    a, b, c, d = score(revcom, tm_target, 1)
    df.loc[len(df)] = [0, 'Reverse Common', str(revcom), '-', str(tmrcom),
                       gccontent(revcom), len(revcom), a, b, c, '-', d]
else:

    a, b, c, d = score(forwcom, tm_target, 1)
    df.loc[len(df)] = [0, 'Forward Common (Both)', str(forwcom), '-',
                       str(tmfcom), gccontent(forwcom), len(forwcom),
                       a, b, c, '-', d]

    a, b, c, d = score(forwcom_left, tm_target, 1)
    df.loc[len(df)] = [0, 'Forward Common (Left)', str(forwcom_left),
                       '-', str(tmfcom_left), gccontent(forwcom_left),
                       len(forwcom_left), a, b, c, '-', d]

    a, b, c, d = score(forwcom_right, tm_target, 1)
    df.loc[len(df)] = [0, 'Forward Common (Right)', str(forwcom_right),
                       '-', str(tmfcom_right), gccontent(forwcom_right),
                       len(forwcom_right), a, b, c, '-', d]

    a, b, c, d = score(revcom, tm_target, 1)
    df.loc[len(df)] = [0, 'Reverse Common (Both)', str(revcom), '-',
                       str(tmrcom), gccontent(revcom), len(revcom),
                       a, b, c, '-', d]

    a, b, c, d = score(revcom_left, tm_target, 1)
    df.loc[len(df)] = [0, 'Reverse Common (Left)', str(revcom_left),
                       '-', str(tmrcom_left), gccontent(revcom_left),
                       len(revcom_left), a, b, c, '-', d]

    a, b, c, d = score(revcom_right, tm_target, 1)
    df.loc[len(df)] = [0, 'Reverse Common (Right)', str(revcom_right),
                       '-', str(tmrcom_right), gccontent(revcom_right),
                       len(revcom_right), a, b, c, '-', d]


for i in range(1, len(protsequence)):
    primername = namer(i, protsequence)
    primer_seq = (nucsequence[nucseqcount[i] - 18:nucseqcount[i]]
                  + 'NNS' + nucsequence[nucseqcount[i] + 3:nucseqcount[i] + 21])
    pboth, pleft, pright = justtrim(nucsequence, primer_seq,
                                    nucseqcount[i] - 18, nucseqcount[i] + 21, tm_target)
    aboth, aleft, aright = logic(nucsequence, 'NNS',
                                 nucseqcount[i], nucseqcount[i] + 3, tm_target, 20)
    if pleft == None:
        revcomp = reversecomplement(pboth, 1)
        primers_right[primername] = primers_left[primername] = primers[primername] = [pboth, revcomp,
                                                                                      tmcalculator(pboth)]
        a, b, c, d, e = score(pboth, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-BF/R', primers[primername][0],
                           primers[primername][1], str(primers[primername][2]),
                           gccontent(primers[primername][0]), len(primers[primername][0]),
                           a, b, c, d, e]

    else:
        revcomp_left = reversecomplement(pleft, 1)
        revcomp_right = reversecomplement(pright, 1)
        revcomp = reversecomplement(pboth, 1)
        primers[primername] = [pboth, revcomp, tmcalculator(pboth)]
        primers_left[primername] = [pleft, revcomp_left, tmcalculator(pleft)]
        primers_right[primername] = [pright, revcomp_right, tmcalculator(pright)]

        a, b, c, d, e = score(pboth, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-BF/R', primers[primername][0],
                           primers[primername][1], str(primers[primername][2]),
                           gccontent(primers[primername][0]), len(primers[primername][0]),
                           a, b, c, d, e]
        a, b, c, d, e = score(pleft, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-LF/R', primers_left[primername][0],
                           primers_left[primername][1], str(primers_left[primername][2]),
                           gccontent(primers_left[primername][0]), len(primers_left[primername][0]),
                           a, b, c, d, e]

        a, b, c, d, e = score(pright, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-RF/R', primers_right[primername][0],
                           primers_right[primername][1], str(primers_right[primername][2]),
                           gccontent(primers_right[primername][0]), len(primers_right[primername][0]),
                           a, b, c, d, e]

    if aleft == None:
        arevcomp = reversecomplement(aboth, 1)
        aprimers_right[primername] = aprimers_left[primername] = aprimers[primername] = [aboth, arevcomp,
                                                                                         tmcalculator(aboth)]
        a, b, c, d, e = score(aboth, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aBF/R', aprimers[primername][0],
                           aprimers[primername][1], str(aprimers[primername][2]),
                           gccontent(aprimers[primername][0]), len(aprimers[primername][0]),
                           a, b, c, d, e]

    else:
        arevcomp_left = reversecomplement(aleft, 1)
        arevcomp_right = reversecomplement(aright, 1)
        arevcomp = reversecomplement(aboth, 1)
        aprimers[primername] = [aboth, arevcomp, tmcalculator(aboth)]
        aprimers_left[primername] = [aleft, arevcomp_left, tmcalculator(aleft)]
        aprimers_right[primername] = [aright, arevcomp_right, tmcalculator(aright)]
        a, b, c, d, e = score(aboth, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aBF/R', aprimers[primername][0],
                           aprimers[primername][1], str(aprimers[primername][2]),
                           gccontent(aprimers[primername][0]), len(aprimers[primername][0]),
                           a, b, c, d, e]
        a, b, c, d, e = score(aleft, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aLF/R', aprimers_left[primername][0],
                           aprimers_left[primername][1], str(aprimers_left[primername][2]),
                           gccontent(aprimers_left[primername][0]), len(aprimers_left[primername][0]),
                           a, b, c, d, e]
        a, b, c, d, e = score(aright, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aRF/R', aprimers_right[primername][0],
                           aprimers_right[primername][1], str(aprimers_right[primername][2]),
                           gccontent(aprimers_right[primername][0]), len(aprimers_right[primername][0]),
                           a, b, c, d, e]

df.sort_values(by=['Residue #', 'Total Score'],
               ascending=[True, False],
               inplace=True).set_index('Residue #')  #

df.to_csv(genename + '-NNS_Primers.csv')
