import pandas as pd

# Global codon translation table
translation = {'ACC': 'T', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I',
               'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'CTC': 'L',
               'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
               'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CAC': 'H',
               'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P',
               'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R',
               'CCA': 'P', 'CGC': 'R', 'GAT': 'D', 'CGG': 'R', 'CTT': 'L',
               'TGC': 'C', 'GGG': 'G', 'TAG': '*', 'GGA': 'G', 'TAA': '*',
               'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L',
               'TTT': 'F', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TGG': 'W',
               'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A',
               'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'TGA': '*',
               'TTG': 'L', 'TCC': 'S', 'TCA': 'S', 'TCT': 'S'}

# Input variables
inputfile = 'input/oxb14-tolc.fasta'
print('\n\nFor finding the CDS start position use start indexing with 1.\n'
      'Position example: aaacccATGacacactgtgtct --> POSITION=7')
# Before the CDS starts at least 60 nucleotide needs to be provided.
# This is important for common primer design
posstartcodon = 519
print('Start codon position: {0}'.format(posstartcodon))
# Next two lines provide region where we define common primer positions.
# Forward primer will be selected from sequence starting from 10th position
# to 50th position.(~15 nucleotide before CDS)
# Likewise common reverse primer will be selected from last 50 nucleotides of the sequence
# excluding last 10 nucleotides
common_fwd_pos = [470, 510]
common_rev_pos = [2010, 2050]

# Read gene sequence and translate
genename = inputfile.split('.')[0]
protsequencefilename = ''.join([inputfile.split('.')[0], '-protein.',
                                inputfile.split('.')[-1]]).replace('input','output')
protsequencefile = open(protsequencefilename, 'w')
nucsequence = ''.join(open(inputfile, 'r').read().split('\n')[1:]).lower()
nucstart = posstartcodon - 1  # -1 comes due to 0 indexing

pos = nucstart
protsequence = ''
while translation[nucsequence[pos:pos + 3].upper()] != '*':
    protsequence += translation[nucsequence[pos:pos + 3].upper()]
    pos += 3

posstopcodon = pos
print('Protein Sequence : ', protsequence)
print('Stop codon position: {0}'.format(posstopcodon))
protsequence = protsequence.capitalize()
protsequencefile.write('>'+genename+'\n'+protsequence)
protsequencefile.close()

nucseqcount = list(range(nucstart, nucstart + 3 * len(protsequence) + 3, 3))

# Design Parameters
tm_target = 63  # Target mean tm temperature (˚C) to set for all primers
srange = 3  # Sensitivity range for tm differences (˚C)

typedict = {'Residue Number': 'int64', 'Primer Name': 'object',
            'Primer Sequence (Fwd)': 'object', 'Primer Sequence (Rev-Comp)': 'object',
            'Tm': 'float64', 'GC Content (%)': 'float64', 'Primer Length': 'int64',
            'GC Score': 'float64', 'Tm Score': 'float64', 'Length Score': 'int64',
            'Wing Score': 'float64', 'Total Score': 'float64'}

df = pd.DataFrame([], columns=['Residue Number', 'Primer Name', 'Primer Sequence (Fwd)',
                               'Primer Sequence (Rev-Comp)', 'Tm', 'GC Content (%)',
                               'Primer Length', 'GC Score', 'Tm Score',
                               'Length Score', 'Wing Score', 'Total Score'])
