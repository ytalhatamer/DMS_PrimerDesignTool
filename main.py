##
# Yusuf Talha TAMER
# UTSW, June 2016
#
# This script can be used to design primers for Saturation Mutagenesis using NEBuilder protocol
##
import os, platform, re
if platform.system()=='Windows':
    os.chdir('C:\\Users\\yusuf\\OneDrive for Business\\Saturation Mutagenesis Primer Design\\NEBuilder_Approach')
else:
    os.chdir('/Users/ytalhatamer/OneDrive - University of Texas Southwestern/Saturation Mutagenesis Primer Design/NEBuilder_Approach')

def reversecomplementer(sequence):
    tmp=sequence.replace('A','1').replace('C','2').replace('G','3').replace('T','4')
    tmp2=tmp.replace('1','T').replace('2','G').replace('3','C').replace('4','A')
    return tmp2[::-1]

def complementer(sequence):
    tmp=sequence.replace('A','1').replace('C','2').replace('G','3').replace('T','4')
    tmp2=tmp.replace('1','T').replace('2','G').replace('3','C').replace('4','A')
    return tmp2
def reverser(sequence):
    return sequence[::-1]

def check5endduplicates(sequence):
    #print(sequence,1)
    if len(sequence.split(sequence[:2]))>2:
        #print (sequence,2)
        return check5endduplicates(sequence[1:])
    else:
        #print(sequence,3)
        return sequence

translation={'ACC':'T','ATG':'M','AAG':'K','AAA':'K','ATC':'I','AAC':'N','ATA':'I','AGG':'R',
             'CCT':'P','CTC':'L','AGC':'S','ACA':'T','AGA':'R','CAT':'H','AAT':'N','ATT':'I',
             'CTG':'L','CTA':'L','ACT':'T','CAC':'H','ACG':'T','CAA':'Q','AGT':'S','CAG':'Q',
             'CCG':'P','CCC':'P','TAT':'Y','GGT':'G','TGT':'C','CGA':'R','CCA':'P','CGC':'R',
             'GAT':'D','CGG':'R','CTT':'L','TGC':'C','GGG':'G','TAG':'*','GGA':'G','TAA':'*',
             'GGC':'G','TAC':'Y','GAG':'E','TCG':'S','TTA':'L','TTT':'F','GAC':'D','CGT':'R',
             'GAA':'E','TGG':'W','GCA':'A','GTA':'V','GCC':'A','GTC':'V','GCG':'A','GTG':'V',
             'TTC':'F','GTT':'V','GCT':'A','TGA':'*','TTG':'L','TCC':'S','TCA':'S','TCT':'S'}

plopen=open('oxb20-tolc.fasta','r')
plseq=('').join(plopen.readlines()[1:])
plseq=plseq.replace('\n','')
plRevComplement=reversecomplementer(plseq)
plComplement=complementer(plseq)
writeresults=open('primersequences49rf.csv','w')
writer=open('primers2see49rf.txt','w')

primerlength=49
plasmidoverlap=[16,13]
primeroverlap=[primerlength-plasmidoverlap[0]-3,primerlength-plasmidoverlap[1]-3]
dupchecklength=10
CDSstart=519
CDSend=2000
codonstartnums=list(range(CDSstart-1,CDSend,3))
CDSseq=plseq[CDSstart-1:CDSend]
proteinsequence=('').join([translation[m] for m in re.findall('...',CDSseq)]) # Protein sequence check
codons=re.findall('...',CDSseq)
writer.write(plseq+'\n'+complementer(plseq)+'\n')

######################## Forward Sequence Maker ####################################
primersfwd={}
for i in range(1,len(codons)-1,2):

    numcodon=codonstartnums[i]
    aanum=((numcodon-CDSstart)/3)+1
    namecodon=str(i+1)+'-'+translation[codons[i]]
    primerseq=plseq[numcodon-primeroverlap[0]:numcodon]+'NNS'+plseq[numcodon+3:numcodon+plasmidoverlap[0]+3]
    duplicatecheckedseq=check5endduplicates(primerseq[:dupchecklength])
    primerseq=duplicatecheckedseq+primerseq[dupchecklength:]
    primerwoNNS=primerseq.split('NNS')
    lens_bfr_after=[len(primerwoNNS[0]) ,len(primerwoNNS[1])]
    primersfwd[namecodon]=[namecodon]
    primersfwd[namecodon].append(i+1)
    primersfwd[namecodon].append(primerseq)
    primersfwd[namecodon].append(reversecomplementer(primerwoNNS[0]))
    primersfwd[namecodon].append(lens_bfr_after[0])
    primersfwd[namecodon].append(lens_bfr_after[1])
    primersfwd[namecodon].append(plseq[2013:2043]) #Plasmid Fwd
    primersfwd[namecodon].append(reversecomplementer(plseq[1983:2043])) #Plasmid Rev
    trimmedseqlen=dupchecklength-len(duplicatecheckedseq)
    writer.write('.'*(numcodon-primeroverlap[0]+trimmedseqlen-len(str(aanum))-3)+str(aanum)+'...'+primerseq+'\n')

    numcodon=codonstartnums[i+1]
    aanum=((numcodon-CDSstart)/3)+1
    namecodon=str(i+2)+'-'+translation[codons[i+1]]
    primerseq=plseq[numcodon-primeroverlap[1]:numcodon]+'NNS'+plseq[numcodon+3:numcodon+plasmidoverlap[1]+3]
    duplicatecheckedseq=check5endduplicates(primerseq[:dupchecklength])
    primerseq=duplicatecheckedseq+primerseq[dupchecklength:]
    primerwoNNS=primerseq.split('NNS')
    lens_bfr_after=[len(primerwoNNS[0]) ,len(primerwoNNS[1])]
    primersfwd[namecodon]=[namecodon]
    primersfwd[namecodon].append(i+2)
    primersfwd[namecodon].append(primerseq)
    primersfwd[namecodon].append(reversecomplementer(primerwoNNS[0]))
    primersfwd[namecodon].append(lens_bfr_after[0])
    primersfwd[namecodon].append(lens_bfr_after[1])
    primersfwd[namecodon].append(plseq[2013:2043]) #Plasmid Fwd
    primersfwd[namecodon].append(reversecomplementer(plseq[1983:2043])) #Plasmid Rev
    trimmedseqlen=dupchecklength-len(duplicatecheckedseq)
    writer.write('.'*(numcodon-primeroverlap[1]+trimmedseqlen-len(str(aanum))-3)+str(aanum)+'...'+primerseq+'\n')


######################## Reverse Sequence Maker ####################################
revcompplseq=reversecomplementer(plseq)
primersrev={}
reversecodons=codons[::-1]
revCDSstart=len(plseq)-CDSend
revCDSend=len(plseq)-CDSstart
revcodonstartnums=list(range(revCDSstart-1,revCDSend,3))
#for i in range(1,len(codons)-1,2):
#    numcodon=revcodonstartnums[i]
#    namecodon=str(len(proteinsequence)-i)+'-'+translation[reversecodons[i]]
#    primerseq=plseq[numcodon-primeroverlap[0]:numcodon]+'NNS'+plseq[numcodon+3:numcodon+plasmidoverlap[0]+3]
#    duplicatecheckedseq=check5endduplicates(primerseq[:dupchecklength])
#    primerseq=duplicatecheckedseq+primerseq[dupchecklength:]
#    primerwoNNS=primerseq.split('NNS')
#    lens_bfr_after=[len(primerwoNNS[0]) ,len(primerwoNNS[1])]
#    primersrev[namecodon]=[namecodon]
#    primersrev[namecodon].append(len(proteinsequence)-i)
#    primersrev[namecodon].append(primerseq)
#    primersrev[namecodon].append(reversecomplementer(primerwoNNS[0]))
#    primersrev[namecodon].append(lens_bfr_after[0])
#    primersrev[namecodon].append(lens_bfr_after[1])
#    primersrev[namecodon].append(reversecomplementer(plseq[472:504])) #Plasmid Rev
#    primersrev[namecodon].append(plseq[471:530]) #Plasmid Fwd
#    trimmedseqlen=dupchecklength-len(duplicatecheckedseq)
#    writer.write('.'*(len(plseq)-numcodon-plasmidoverlap[0])+reverser(primerseq)+'\n')
#
#    numcodon=revcodonstartnums[i+1]
#    namecodon=str(len(proteinsequence)-i-1)+'-'+translation[reversecodons[i+1]]
#    primerseq=plseq[numcodon-primeroverlap[1]:numcodon]+'NNS'+plseq[numcodon+3:numcodon+plasmidoverlap[1]+3]
#    duplicatecheckedseq=check5endduplicates(primerseq[:dupchecklength])
#    primerseq=duplicatecheckedseq+primerseq[dupchecklength:]
#    primerwoNNS=primerseq.split('NNS')
#    lens_bfr_after=[len(primerwoNNS[0]) ,len(primerwoNNS[1])]
#    primersrev[namecodon]=[namecodon]
#    primersrev[namecodon].append(len(proteinsequence)-i-1)
#    primersrev[namecodon].append(primerseq)
#    primersrev[namecodon].append(reversecomplementer(primerwoNNS[0]))
#    primersrev[namecodon].append(lens_bfr_after[0])
#    primersrev[namecodon].append(lens_bfr_after[1])
#    primersrev[namecodon].append(reversecomplementer(plseq[472:504])) #Plasmid Rev
#    primersrev[namecodon].append(plseq[471:530]) #Plasmid Fwd
#    trimmedseqlen=dupchecklength-len(duplicatecheckedseq)
#    writer.write('.'*(len(plseq)-numcodon-plasmidoverlap[1])+reverser(primerseq)+'\n')
for i in range(1,len(codons)-1,2):
    numcodon=revcodonstartnums[i]
    namecodon=str(len(proteinsequence)-i)+'-'+translation[reversecodons[i]]
    primerseq=revcompplseq[numcodon-primeroverlap[0]:numcodon+1]+'NNS'+revcompplseq[numcodon+4:numcodon+plasmidoverlap[0]+3]
    duplicatecheckedseq=check5endduplicates(primerseq[:dupchecklength])
    primerseq=duplicatecheckedseq+primerseq[dupchecklength:]
    primerwoNNS=primerseq.split('NNS')
    lens_bfr_after=[len(primerwoNNS[0]) ,len(primerwoNNS[1])]
    primersrev[namecodon]=[namecodon]
    primersrev[namecodon].append(len(proteinsequence)-i)
    primersrev[namecodon].append(primerseq)
    primersrev[namecodon].append(reversecomplementer(primerwoNNS[0]))
    primersrev[namecodon].append(lens_bfr_after[0])
    primersrev[namecodon].append(lens_bfr_after[1])
    primersrev[namecodon].append(reversecomplementer(plseq[472:504])) #Plasmid Rev
    primersrev[namecodon].append(plseq[471:530]) #Plasmid Fwd
    trimmedseqlen=dupchecklength-len(duplicatecheckedseq)
    writer.write('.'*(len(plseq)-numcodon-plasmidoverlap[0])+reverser(primerseq)+'\n')

    numcodon=revcodonstartnums[i+1]
    namecodon=str(len(proteinsequence)-i-1)+'-'+translation[reversecodons[i+1]]
    primerseq=revcompplseq[numcodon-primeroverlap[1]:numcodon+1]+'NNS'+revcompplseq[numcodon+4:numcodon+plasmidoverlap[1]+3]
    duplicatecheckedseq=check5endduplicates(primerseq[:dupchecklength])
    primerseq=duplicatecheckedseq+primerseq[dupchecklength:]
    primerwoNNS=primerseq.split('NNS')
    lens_bfr_after=[len(primerwoNNS[0]) ,len(primerwoNNS[1])]
    primersrev[namecodon]=[namecodon]
    primersrev[namecodon].append(len(proteinsequence)-i-1)
    primersrev[namecodon].append(primerseq)
    primersrev[namecodon].append(reversecomplementer(primerwoNNS[0]))
    primersrev[namecodon].append(lens_bfr_after[0])
    primersrev[namecodon].append(lens_bfr_after[1])
    primersrev[namecodon].append(reversecomplementer(plseq[472:504])) #Plasmid Rev
    primersrev[namecodon].append(plseq[471:530]) #Plasmid Fwd
    trimmedseqlen=dupchecklength-len(duplicatecheckedseq)
    writer.write('.'*(len(plseq)-numcodon-plasmidoverlap[1])+reverser(primerseq)+'\n')






writeresults.write('ID,Residue #,Mutant Primer,Reverse Primer,\
                Primer Overlap Length,Plasmid Overlap Length,Plasmid Fwd,Plasmid Rev,F/R\n')
for i in primersfwd.values():
    script=str(i).replace('[','').replace(']','').replace('\'','')
    writeresults.write(script+',forward'+'\n')

for i in primersrev.values():
    script=str(i).replace('[','').replace(']','').replace('\'','')
    writeresults.write(script+',reverse'+'\n')

writeresults.close()
