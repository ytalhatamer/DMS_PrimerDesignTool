translation={'ACC':'T','ATG':'M','AAG':'K','AAA':'K','ATC':'I','AAC':'N','ATA':'I','AGG':'R',
             'CCT':'P','CTC':'L','AGC':'S','ACA':'T','AGA':'R','CAT':'H','AAT':'N','ATT':'I',
             'CTG':'L','CTA':'L','ACT':'T','CAC':'H','ACG':'T','CAA':'Q','AGT':'S','CAG':'Q',
             'CCG':'P','CCC':'P','TAT':'Y','GGT':'G','TGT':'C','CGA':'R','CCA':'P','CGC':'R',
             'GAT':'D','CGG':'R','CTT':'L','TGC':'C','GGG':'G','TAG':'*','GGA':'G','TAA':'*',
             'GGC':'G','TAC':'Y','GAG':'E','TCG':'S','TTA':'L','TTT':'F','GAC':'D','CGT':'R',
             'GAA':'E','TGG':'W','GCA':'A','GTA':'V','GCC':'A','GTC':'V','GCG':'A','GTG':'V',
             'TTC':'F','GTT':'V','GCT':'A','TGA':'*','TTG':'L','TCC':'S','TCA':'S','TCT':'S'}


inputfile = 'tolC.fasta'
genename=inputfile.split('.')[0]
protsequencefilename=''.join([inputfile.split('.')[:-1], '-protein.', inputfile.split('.')[-1]])
protsequencefile=open(protsequencefilename, 'w')
nucsequence = open(inputfile, 'r').read().split('\n')[1].lower()
posstartcodon = 64
nucstart=posstartcodon-1 # -1 comes due to 0 indexing
pos=nucstart
protsequence=''
while translation[nucsequence[pos:pos+3]]!='*':
    protsequence+=translation[nucsequence[pos:pos+3]]
    pos+=3
protsequence=protsequence.capitalize()

posstopcodon = [i for i in range(64,len(nucsequence),3) if translation[nucsequence[i:i+3]]=='*'][0]

tm_target=63 # Target mean tm temperature to set for all primers