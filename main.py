from functions import *
import os,time,sys
import pandas as pd

'''
If you run this file (main.py) it will ask you two files and two variables.

Prepare two files. 

1. Nucleotide file. Take your gene sequence with at least 60 nucleotide before and after the gene. 
Leading and trailing sequences are important for designing the common primers. (Check figure-1 for 
why we use common primers)

  - You can change the default positions of searching space for best common primer. 
  Now takes 40 nucleotides (position 10:50(line 70) and -50:-10 (line 82)) for forward common and reverse common primers. 
  
  - We try to adjust the Tm of each primer designed as close to TARGET Tm temperature. (You can change target Tm at line 60).  

2. Protein file. Just provide a fasta file of your protein sequence. 

Variables:
1. Location of the fasta files. 
  - I recommend you put them in the same folder with main.py and functions.py files. 

2. Position of your start codon in nucleotide file.(Start counting from 1. (check the example))



This script will design primers for each residue starting from second residue to last residue before stop codon.

It will output 4 alternative primers for each residue. I scored them by 4 criteria.
1. GC score     : 50% GC means 25 points and any deviation decreases the score. (0:25 points)
2. Tm score     : I defined ±3 degrees is OK for Tm but any temperature away from Tm±3 will be punished (-Inf:25 points) 
3. Length score : Optimal primer length is in between 30-40 nucleotides. Longer and shorted primers will be punished (-Inf:25 points)
4. Wing score   : We want NNS triplet to be in the middle of the primer. 
                  Optimum primer has the same number of nucleotides before NNS (left wing) 
                  and after NNS (right wing)  (-Inf:25 points)

  Total score : It is the sum of all four scores above. Scoring scheme can be modified and might improve the design process.

We used the primer that gave highest score for each residue. 

If you want to use forward common primer below are the primer sets you need to use:
    INSERT PCR primers:
          1. Forward Common primer (Forward)
          2. Reverse NNS primer (Reverse)
    PLASMID PCR primers:
          1. Reverse complement of right wing of reverse NNS primer (Forward)
          2. Reverse complement of Forward common primer (Reverse)

If you want to use reverse common primer below are the primer sets you need to use:
    INSERT PCR primers:
          1. Forward NNS primer (Reverse)
          2. Reverse Common primer (Forward)
    PLASMID PCR primers:
          1. Reverse complement of Reverse common primer (Reverse)
          2. Reverse complement of left wing of reverse NNS primer (Forward)

*****************************************************************
Just to give you an example, 
I provided two files. tolC_gene.fasta, tolC_protein.fasta 
You can run main.py and gene starts at position 64.
That will give you the csv file I also provided in the folder.
*****************************************************************
'''
print('You can find location info with right click and get info/properties option')

location=input('Please provide location of folder: ')

print('Prepare the nucleotide sequence of the protein(in fasta format)\
Be careful to put at least 60 nucleotide more before and after the start/stop codon')

print('Prepare the protein sequence of the protein(in fasta format)')


time.sleep(1) 


genename=input('Please provide genename: ').strip()

print ('Change the names of the sequence files to '+genename+'_prot.fasta and '+genename+'_nuc.fasta \
and put these two files into same folder with These python files.')

protfilename=genename+'_prot.fasta'
nucfilename=genename+'_nuc.fasta'
os.chdir(location)
try:
    ppseq=open(protfilename,'r')
    try:
        nnseq=open(nucfilename,'r')
    except:
        exitmessage_nuc='!!!\nMake sure you put the gene sequence file into correct folder with correct name\
        \nGiven Folder Location:'+location+'\nCurrent Directory:'+os.getcwd()+'\n Correct File Name: '+nucfilename
        sys.exit(exitmessage_nuc) 
except:
    exitmessage_prot='!!!\nMake sure you put the protein sequence file into correct folder with correct name\
    \nGiven Folder Location:'+location+'\nCurrent Directory:'+os.getcwd()+'\n Correct File Name: '+protfilename
    sys.exit(exitmessage_prot)

ppseq.readline()
protseq=(''.join(ppseq.readlines())).replace('\n','')
protseq=protseq.capitalize()

nnseq.readline()
nucseq=(''.join(nnseq.readlines())).replace('\n','')
nucseq=nucseq.lower()
print ('\n\nPosition example: aaacccATGacacactgtgtct --> POSITION=7')
nucstart=int(input('Provide the position of A in ...ATG...: '))-1

df=pd.DataFrame([],columns=['Residue #','Primer Name','Primer Seq(Fwd)','Primer Sequence (Rev-Comp)','Tm',\
                            'GC Content(%)','Primer Length','GC Score','Tm Score','Length Score','Wing Score','Total Score'])

tm_target=63 # Target mean tm temperature to set for all primers

primers={}
primers_left={}
primers_right={}

aprimers={}
aprimers_left={}
aprimers_right={}

forwcom,forwcom_left,forwcom_right=justtrim(nucseq,nucseq[10:50],10,50,tm_target) # Common forward primer for all reactions
if forwcom_left==None:
    tmfcom=tmcalculator(forwcom)
    primers_right['forward']=primers_left['forward']=primers['forward']=[forwcom,None,tmfcom]
else:
    tmfcom_left=tmcalculator(forwcom_left)
    tmfcom_right=tmcalculator(forwcom_right)
    tmfcom=tmcalculator(forwcom)
    primers['forward']=[forwcom,None,tmfcom]
    primers_left['forward']=[forwcom_left,None,tmfcom_left]
    primers_right['forward']=[forwcom_right,None,tmfcom_right]

revcom,revcom_left,revcom_right=justtrim(nucseq,nucseq[-50:-10],len(nucseq)-50,len(nucseq)-10,tm_target) # Common reverse primer for all reactions

if revcom_left==None:
    revcom = reversecomplement(revcom,0)
    tmrcom=tmcalculator(revcom)
    primers_right['reverse']=primers_left['reverse']=primers['reverse']=[revcom,None,tmrcom]

else:
    revcom_left = reversecomplement(revcom_left,0)
    revcom_right = reversecomplement(revcom_right,0)
    revcom = reversecomplement(revcom,0)
    tmrcom_left=tmcalculator(revcom_left)
    tmrcom_right=tmcalculator(revcom_right)
    tmrcom=tmcalculator(revcom)
    
    primers['reverse']=[revcom,None,tmrcom]
    primers_left['reverse']=[revcom_left,None,tmrcom_left]
    primers_right['reverse']=[revcom_right,None,tmrcom_right]

if forwcom_left==None:
    a,b,c,d=score(forwcom,tm_target,1)
    df.loc[len(df)]=[0,'Forward Common',str(forwcom),'-',str(tmfcom),percgc(forwcom),str(len(forwcom)),a,b,c,
                '-',d]

    a,b,c,d=score(revcom,tm_target,1)
    df.loc[len(df)] = [0, 'Reverse Common', str(revcom), '-', str(tmrcom), percgc(revcom), str(len(revcom)), a, b, c,
                  '-', d]
else:

    a,b,c,d=score(forwcom,tm_target,1)
    df.loc[len(df)] = [0, 'Forward Common (Both)', str(forwcom), '-', str(tmfcom), percgc(forwcom), str(len(forwcom)), a, b, c,
                  '-', d]

    a,b,c,d=score(forwcom_left,tm_target,1)
    df.loc[len(df)] = [0, 'Forward Common (Left)', str(forwcom_left), '-', str(tmfcom_left), percgc(forwcom_left), str(len(forwcom_left)), a, b, c,
                  '-', d]

    a,b,c,d=score(forwcom_right,tm_target,1)
    df.loc[len(df)] = [0, 'Forward Common (Right)', str(forwcom_right), '-', str(tmfcom_right), percgc(forwcom_right), str(len(forwcom_right)), a, b, c,
                  '-', d]

    a,b,c,d=score(revcom,tm_target,1)
    df.loc[len(df)] = [0, 'Reverse Common (Both)', str(revcom), '-', str(tmrcom), percgc(revcom), str(len(revcom)), a, b, c,
                  '-', d]

    a,b,c,d=score(revcom_left,tm_target,1)
    df.loc[len(df)] = [0, 'Reverse Common (Left)', str(revcom_left), '-', str(tmrcom_left), percgc(revcom_left), str(len(revcom_left)), a, b, c,
                  '-', d]

    a,b,c,d=score(revcom_right,tm_target,1)
    df.loc[len(df)] = [0, 'Reverse Common (Right)', str(revcom_right), '-', str(tmrcom_right), percgc(revcom_right), str(len(revcom_right)), a, b, c,
                  '-', d]

nucseqcount=list(range(nucstart,nucstart+3*len(protseq)+3,3))
for i in range(1,len(protseq)):    
    primername=namer(i,protseq)
    primer_seq=(nucseq[nucseqcount[i]-18:nucseqcount[i]]+'NNS'+nucseq[nucseqcount[i]+3:nucseqcount[i]+21])
    addseq='NNS'
    pboth,pleft,pright=justtrim(nucseq, primer_seq,nucseqcount[i]-18,nucseqcount[i]+21,tm_target)
    aboth,aleft,aright=logic(nucseq,addseq,nucseqcount[i],nucseqcount[i]+3,tm_target,20)
    if pleft==None:
        revcomp = reversecomplement(pboth,1)
        primers_right[primername]=primers_left[primername]=primers[primername]=[pboth, revcomp,tmcalculator(pboth)]
        a,b,c,d,e=score(pboth,tm_target)
        df.loc[len(df)] = [i+1, primername+'-BF/R', primers[primername][0], primers[primername][1], str(primers[primername][2]),
                           percgc(primers[primername][0]), str(len(primers[primername][0])), a, b, c,
                           d, e]

    else:
        revcomp_left = reversecomplement(pleft,1)
        revcomp_right = reversecomplement(pright,1)
        revcomp = reversecomplement(pboth,1)
        primers[primername]=[pboth, revcomp,tmcalculator(pboth)]
        primers_left[primername]=[pleft, revcomp_left,tmcalculator(pleft)]
        primers_right[primername]=[pright, revcomp_right,tmcalculator(pright)]

        a,b,c,d,e=score(pboth,tm_target)
        df.loc[len(df)] = [i + 1, primername + '-BF/R', primers[primername][0], primers[primername][1],
                           str(primers[primername][2]),
                           percgc(primers[primername][0]), str(len(primers[primername][0])), a, b, c,
                           d, e]
        a,b,c,d,e=score(pleft,tm_target)
        df.loc[len(df)] = [i + 1, primername + '-LF/R', primers_left[primername][0], primers_left[primername][1],
                           str(primers_left[primername][2]),
                           percgc(primers_left[primername][0]), str(len(primers_left[primername][0])), a, b, c,
                           d, e]

        a,b,c,d,e=score(pright,tm_target)
        df.loc[len(df)] = [i+1, primername+'-RF/R', primers_right[primername][0], primers_right[primername][1], str(primers_right[primername][2]),
                           percgc(primers_right[primername][0]), str(len(primers_right[primername][0])), a, b, c,
                           d, e]

    if aleft==None:
        arevcomp = reversecomplement(aboth,1)
        aprimers_right[primername]=aprimers_left[primername]=aprimers[primername]=[aboth, arevcomp,tmcalculator(aboth)]
        a,b,c,d,e=score(aboth,tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aBF/R', aprimers[primername][0], aprimers[primername][1],
                           str(aprimers[primername][2]),
                           percgc(aprimers[primername][0]), str(len(aprimers[primername][0])), a, b, c,
                           d, e]

    else:
        arevcomp_left = reversecomplement(aleft,1)
        arevcomp_right = reversecomplement(aright,1)
        arevcomp = reversecomplement(aboth,1)
        aprimers[primername]=[aboth, arevcomp,tmcalculator(aboth)]
        aprimers_left[primername]=[aleft, arevcomp_left,tmcalculator(aleft)]
        aprimers_right[primername]=[aright, arevcomp_right,tmcalculator(aright)]
        a,b,c,d,e=score(aboth,tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aBF/R', aprimers[primername][0], aprimers[primername][1],
                           str(aprimers[primername][2]),
                           percgc(aprimers[primername][0]), str(len(aprimers[primername][0])), a, b, c,
                           d, e]
        a,b,c,d,e=score(aleft,tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aLF/R', aprimers_left[primername][0], aprimers_left[primername][1],
                           str(aprimers_left[primername][2]),
                           percgc(aprimers_left[primername][0]), str(len(aprimers_left[primername][0])), a, b, c,
                           d, e]
        a,b,c,d,e=score(aright,tm_target)
        df.loc[len(df)] = [i + 1, primername + '-aRF/R', aprimers_right[primername][0], aprimers_right[primername][1],
                           str(aprimers_right[primername][2]),
                           percgc(aprimers_right[primername][0]), str(len(aprimers_right[primername][0])), a, b, c,
                           d, e]

df2=df.sort_values(by=['Residue #','Total Score'],ascending=[True,False]).set_index('Residue #')#

df2.to_csv(genename+'-NNS_Primers.csv')
