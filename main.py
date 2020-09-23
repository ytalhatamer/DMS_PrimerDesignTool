import numpy as np
from globalvars import nucsequence, protsequence, common_fwd_pos, common_rev_pos, \
    tm_target, genename, nucseqcount, df, typedict
from functions import namer, tmcalculator, gccontent, reversecomplement, \
    logic, justtrim, score

primers = dict(right=dict(), left=dict(), both=dict())

aprimers = dict(right=dict(), left=dict(), both=dict())

# Before the CDS starts at least 60 nucleotide needs to be provided.
# This is important for common primer design

# COMMON PRIMER DESIGN
# Forward Common Primers
# Common forward primer for all reactions are send for optimization.
forwcom, forwcom_left, forwcom_right = justtrim(nucsequence,
                                                nucsequence[common_fwd_pos[0]:common_fwd_pos[1]],
                                                common_fwd_pos[0],
                                                common_fwd_pos[1],
                                                tm_target)

if forwcom_left:
    for forward_common, pos in zip([forwcom, forwcom_left, forwcom_right], ['both', 'left', 'right']):
        tm_common_forward = tmcalculator(forward_common)
        primers[pos]['forward'] = [forward_common, None, tm_common_forward]

        a, b, c, d = score(forward_common, tm_target, 1)  # 1 is a flag to score function telling no NNS in primer.
        df.loc[len(df)] = [0, 'Forward Common (' + pos + ')', str(forward_common),
                           reversecomplement(forward_common, 0), str(tm_common_forward),
                           gccontent(forward_common), len(forward_common),
                           a, b, c, np.nan, d]
else:
    tm_common_forward = tmcalculator(forwcom)
    primers['both']['forward'] = \
        primers['right']['forward'] = \
        primers['left']['forward'] = [forwcom, None, tm_common_forward]

    a, b, c, d = score(forwcom, tm_target, 1)
    df.loc[len(df)] = [0, 'Forward Common', str(forwcom),
                       reversecomplement(forwcom, 0), str(tm_common_forward),
                       gccontent(forwcom), len(forwcom), a, b, c, np.nan, d]

# Designing Reverse Common Primers
revcom, revcom_left, revcom_right = justtrim(nucsequence,
                                             nucsequence[common_rev_pos[0]:common_rev_pos[1]],
                                             common_rev_pos[0],
                                             common_rev_pos[1],
                                             tm_target)  # Common reverse primer for all reactions

if revcom_left:
    for reverse_common, pos in zip([revcom, revcom_left, revcom_right], ['both', 'left', 'right']):
        rev_comp_sequence = reversecomplement(reverse_common, 0)
        tm_common_reverse = tmcalculator(reverse_common)
        primers[pos]['reverse'] = [reverse_common, None, tm_common_reverse]

        a, b, c, d = score(reverse_common, tm_target, 1)
        df.loc[len(df)] = [0, 'Reverse Common (' + pos + ')',
                           str(reverse_common), reversecomplement(reverse_common, 0),
                           str(tm_common_reverse), gccontent(reverse_common),
                           len(reverse_common), a, b, c, np.nan, d]
else:
    rev_comp_sequence = reversecomplement(revcom, 0)
    tm_common_reverse = tmcalculator(rev_comp_sequence)
    primers['both']['reverse'] = \
        primers['right']['reverse'] = \
        primers['left']['reverse'] = [rev_comp_sequence, None, tm_common_reverse]

    a, b, c, d = score(revcom, tm_target, 1)
    df.loc[len(df)] = [0, 'Reverse Common', str(revcom),
                       reversecomplement(revcom, 0), str(tm_common_reverse),
                       gccontent(revcom), len(revcom), a, b, c, np.nan, d]

# This part we will find and optimize NNS primers. Next for loop goes over all codons
# starting from second codon and turns each codon to NNS. Takes ~20-25 nucleotides
# before and after NNS tries to optimize primer quality by trimming or adding.
# there are 4 categories of primer quality measure that is evaluated.
# More details on scoring can be found in score function in functions.py

for i in range(1, len(protsequence)):
    primername = namer(i, protsequence)
    primer_seq = (nucsequence[nucseqcount[i] - 18:nucseqcount[i]]
                  + 'NNS' + nucsequence[nucseqcount[i] + 3:nucseqcount[i] + 21])
    pboth, pleft, pright = justtrim(nucsequence, primer_seq,
                                    nucseqcount[i] - 18, nucseqcount[i] + 21, tm_target)
    aboth, aleft, aright = logic(nucsequence, 'NNS',
                                 nucseqcount[i], nucseqcount[i] + 3, tm_target, 20)
    if pleft:
        for curr_primer, p_extension, pos in zip([pboth, pleft, pright],
                                                 ['-BF/R', '-LF/R', '-RF/R'],
                                                 ['both', 'left', 'right']):
            revcomp = reversecomplement(curr_primer, 1)
            primers[pos][primername] = [curr_primer, revcomp, tmcalculator(curr_primer)]
            a, b, c, d, e = score(curr_primer, tm_target)
            df.loc[len(df)] = [i + 1, primername + p_extension,
                               primers[pos][primername][0],
                               primers[pos][primername][1],
                               str(primers[pos][primername][2]),
                               gccontent(primers[pos][primername][0]),
                               len(primers[pos][primername][0]),
                               a, b, c, d, e]
    else:
        revcomp = reversecomplement(pboth, 1)
        primers['right'][primername] = \
            primers['left'][primername] = \
            primers['both'][primername] = [pboth, revcomp, tmcalculator(pboth)]
        a, b, c, d, e = score(pboth, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-BF/R',
                           primers['both'][primername][0],
                           primers['both'][primername][1],
                           str(primers['both'][primername][2]),
                           gccontent(primers['both'][primername][0]),
                           len(primers['both'][primername][0]),
                           a, b, c, d, e]
    if aleft:
        for curr_primer, p_extension, pos in zip([aboth, aleft, aright],
                                                 ['-aBF/R', '-aLF/R', '-aRF/R'],
                                                 ['both', 'left', 'right']):
            arevcomp = reversecomplement(curr_primer, 1)
            aprimers[pos][primername] = [curr_primer, arevcomp, tmcalculator(curr_primer)]
            a, b, c, d, e = score(curr_primer, tm_target)
            df.loc[len(df)] = [i + 1, primername + p_extension,
                               aprimers[pos][primername][0],
                               aprimers[pos][primername][1],
                               str(aprimers[pos][primername][2]),
                               gccontent(aprimers[pos][primername][0]),
                               len(aprimers[pos][primername][0]),
                               a, b, c, d, e]

    else:
        arevcomp = reversecomplement(aboth, 1)
        aprimers['right'][primername] = \
            aprimers['left'][primername] = \
            aprimers['both'][primername] = [aboth, arevcomp, tmcalculator(aboth)]
        a, b, c, d, e = score(aboth, tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aBF/R',
                           aprimers['both'][primername][0],
                           aprimers['both'][primername][1],
                           str(aprimers['both'][primername][2]),
                           gccontent(aprimers['both'][primername][0]),
                           len(aprimers['both'][primername][0]),
                           a, b, c, d, e]

df = df.sort_values(by=['Residue Number', 'Total Score', 'Primer Name'],
                    ascending=[True, False, False])  #

df = df.astype(dtype=typedict)
rounded_cols = ['Tm', 'GC Content (%)', 'GC Score',
                'Tm Score', 'Wing Score', 'Total Score']
df[rounded_cols] = df[rounded_cols].round(2)
df = df.loc[~(df.duplicated(subset=['Residue Number', 'Primer Sequence (Fwd)'], keep='last')), :]
df.to_csv('output/' + genename + '-NNS_Primers.csv', index=False)
