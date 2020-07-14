# Deep Mutational Scanning Primer Design Toolbox 


Primer Design Tool for making Deep Mutational Scanning Libraries of genes/proteins.

One file and three variables needs to be provided to successfully run this file (main.py).

1. Nucleotide file. Take your gene sequence with at least 60 nucleotide before and after the gene. 
Leading and trailing sequences are important for designing the common primers. (Check figure-1 for 
why we use common primers)

  - The default positions of  for best common primer searching space can be modified. 
  - Adjusted default primer length 40 nucleotides for forward common and reverse common primers.
        (position 10:50(line 70) and -50:-10 (line 82)) These positions can be modified by users choice. 
  - We try to adjust the Tm of each primer designed as close to TARGET Tm temperature. 
        (Target Tm can be modified at line 60).  

2. Protein file. Just provide a fasta file of your protein sequence. 

Variables:
1. Location of the fasta files. 
  - I recommend you put them in the same folder with main.py and functions.py files. 
2. Position of your start codon in nucleotide file.(Start counting from 1.)
3. Position of your stop codon in nucleotide file. 

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
