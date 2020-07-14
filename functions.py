
srange=3 # Sensitivity range for tm differences (˚C)
######################################################
def tmcalculator(seq):
    numgc=seq.count('c')+seq.count('g')
    tm=64.9+41*((numgc-16.4)/len(seq))
    tm=round(tm,2)
    return tm

#####################################################
#####################################################
def reversecomplement(seq,flag):
    tmp=[]
    if flag==0:
        sj=seq
        ssj=sj.replace('a','1').replace('g','2').replace('t','3').replace('c','4')
        k=ssj.replace('1','t').replace('2','c').replace('3','a').replace('4','g')
        revseq=k[::-1]
    else:
        sj=seq.split('NNS')
        for i in sj:
            ssj=i.replace('a','1').replace('g','2').replace('t','3').replace('c','4')
            #print ssj
            k=ssj.replace('1','t').replace('2','c').replace('3','a').replace('4','g')
            tmp.append(k)
        revseq=tmp[1][::-1]+'SNN'+tmp[0][::-1]
    return revseq

####################################################
def checktm(tm_int,tm_target):
    if tm_int<(tm_target-srange):
        return 1
    else:
        if tm_int<(tm_target+srange):
            return 2
        else:
            return 3

####################################################
def namer(count,proteinseq,optname=None):
    l=len(str(len(proteinseq)))
    if optname==None:
        primname='p_'+'0'*(l-len(str(count+1)))+str(count+1)+'_'+proteinseq[count]
    else:
        primname=optname
    return primname

####################################################

def trimleft(nucseq,seq,p_start,p_end,tm_target):
    # This function starts trimming the sequence from left then right consequtively...
    k=0
    while k<10 and checktm(tmcalculator(seq),tm_target)!=2:

        if k%2==0:
            seq=seq[1:]
            p_start+=1
            for m in range(10):
                if seq[0]=='c' or seq[0]=='g':
                    break
                seq=seq[1:]
                p_start+=1
            for m in range(10):
                if seq[-1]=='c' or seq[-1]=='g': 
                    break
                seq=seq[:-1]
                p_end-=1
        else:
            seq=seq[:-1]
            p_end-=1
            for m in range(10):
                if seq[0]=='c' or seq[0]=='g':
                    break
                seq=seq[1:]
                p_start+=1
            for m in range(10):
                if seq[-1]=='c' or seq[-1]=='g': 
                    break
                seq=seq[:-1]
                p_end-=1
        # print 'TRIMLEFT-Trial',k,'SEQUENCE',seq
        k+=1
    return seq

def trimright(nucseq,seq,p_start,p_end,tm_target):
    # This function starts trimming the sequence from right then left consequtively...
    k=0
    while k<10 and checktm(tmcalculator(seq),tm_target)!=2:

        if k%2==0:
            seq=seq[:-1]
            p_end-=1
            for m in range(10):
                if seq[0]=='c' or seq[0]=='g':
                    break
                seq=seq[1:]
                p_start+=1
            for m in range(10):
                if seq[-1]=='c' or seq[-1]=='g': 
                    break
                seq=seq[:-1]
                p_end-=1
        else:
            seq=seq[1:]
            p_start+=1
            for m in range(10):
                if seq[0]=='c' or seq[0]=='g':
                    break
                seq=seq[1:]
                p_start+=1
            for m in range(10):
                if seq[-1]=='c' or seq[-1]=='g': 
                    break
                seq=seq[:-1]
                p_end-=1
        # print 'TRIMRIGHT-Trial',k,'SEQUENCE',seq
        k+=1
    return seq

def trimboth(nucseq,seq,p_start,p_end,tm_target):
    k=0
    while k<10 and checktm(tmcalculator(seq),tm_target)!=2:

        seq=seq[1:-1]
        p_start+=1
        p_end-=1
        for m in range(10):
            if seq[0]=='c' or seq[0]=='g':
                break
            seq=seq[1:]
            p_start+=1
        for m in range(10):
            if seq[-1]=='c' or seq[-1]=='g': 
                break
            seq=seq[:-1]
            p_end-=1
        # print 'TRIMBOTH-Trial',k,'SEQUENCE',seq
        k+=1
    return seq

def addleft(nucseq,seq,p_start,p_end,tm_target,kup=10):
    k=0
    while k<kup and checktm(tmcalculator(seq),tm_target)!=2:
        if len(seq)>=39:
            break
        if k%2==0:
            seq=nucseq[p_start-1]+seq
            p_start-=1
            for m in range(10):
                if seq[0]=='c' or seq[0]=='g':
                    break
                seq=nucseq[p_start-1]+seq
                p_start-=1
            for m in range(10):
                if seq[-1]=='c' or seq[-1]=='g': 
                    break
                seq=seq+nucseq[p_end+1]
                p_end+=1
        else:
            # print ('!!!!!! P_END IS HERE !!!!!!!!! --> %d',p_end)
            seq=seq+nucseq[p_end+1]
            p_end+=1
            for m in range(10):
                if seq[0]=='c' or seq[0]=='g':
                    break
                seq=nucseq[p_start-1]+seq
                p_start-=1
            for m in range(10):
                if seq[-1]=='c' or seq[-1]=='g': 
                    break
                seq=seq+nucseq[p_end+1]
                p_end+=1
        #print 'ADDLEFT-Trial',k,'SEQUENCE',seq
        k+=1
        if len(seq)>=39:
            break
    return seq

def addright(nucseq,seq,p_start,p_end,tm_target,kup=10):

    k=0
    while k<kup and checktm(tmcalculator(seq),tm_target)!=2:
        if len(seq)>=39:
            break
        if k%2==0:
            seq=seq+nucseq[p_end+1]
            p_end+=1
            for m in range(10):
                if seq[0]=='c' or seq[0]=='g':
                    break
                seq=nucseq[p_start-1]+seq
                p_start-=1
            for m in range(10):
                if seq[-1]=='c' or seq[-1]=='g': 
                    break
                seq=seq+nucseq[p_end+1]
                p_end+=1
        else:
            seq=nucseq[p_start-1]+seq
            p_start-=1
            for m in range(10):
                if seq[0]=='c' or seq[0]=='g':
                    break
                seq=nucseq[p_start-1]+seq
                p_start-=1
            for m in range(10):
                if seq[-1]=='c' or seq[-1]=='g': 
                    break
                seq=seq+nucseq[p_end+1]
                p_end+=1
        #print 'ADDRIGHT-Trial',k,'SEQUENCE',seq
        k+=1
        if len(seq)>=39:
            break
    return seq

def addboth(nucseq,seq,p_start,p_end,tm_target,kup=10):
    k=0
    while k<kup and checktm(tmcalculator(seq),tm_target)!=2:
        if len(seq)>=39:
            break
        seq=nucseq[p_start-1]+seq+nucseq[p_end+1]
        p_start-=1
        p_end+=1
        for m in range(10):
            if seq[0]=='c' or seq[0]=='g':
                break
            seq=seq+nucseq[p_start-1]
            p_start-=1
        for m in range(10):
            if seq[-1]=='c' or seq[-1]=='g': 
                break
            seq=seq+nucseq[p_end+1]
            p_end+=1
        #print 'ADDBOTH-Trial',k,'SEQUENCE',seq
        k+=1
        if len(seq)>=39:
            break
    return seq


################################################################
def logic(nucseq,seq,p_start,p_end,tm_target,k=10):
    if tmcalculator(seq)<(tm_target-srange):
        seqleft= addleft(nucseq,seq,p_start,p_end,tm_target,k)
        seqright= addright(nucseq,seq,p_start,p_end,tm_target,k)
        seqboth= addboth(nucseq,seq,p_start,p_end,tm_target,k)
        # print 'ADDRIGHT-->',seqright
        # print 'ADDLEFT-->',seqleft
        # print 'ADDBOTH-->',seqboth
        return seqboth,seqleft,seqright
    else:
        if tmcalculator(seq)<(tm_target+srange):
            seqleft=None
            seqright=None
            # print 'RIGHT-->',seqright
            # print 'LEFT-->',seqleft
            # print 'BOTH-->',seq
            return seq,seqleft,seqright
        else:
            seqleft= trimleft(nucseq,seq,p_start,p_end,tm_target)
            seqright= trimright(nucseq,seq,p_start,p_end,tm_target)
            seqboth= trimboth(nucseq,seq,p_start,p_end,tm_target)
            # print 'TRIMRIGHT-->',seqright
            # print 'TRIMLEFT-->',seqleft
            # print 'TRIMBOTH-->',seqboth
            return seqboth,seqleft,seqright

################################################################

def justtrim(nucseq,seq,p_start,p_end,tm_target):
    # print('JUSTTRIM Bfr-->',seq,len(seq))
    for m in range(10):
        if seq[0]=='c' or seq[0]=='g':
            break
        seq=seq[1:]
        p_start+=1
    for m in range(10):
        if seq[-1]=='c' or seq[-1]=='g': 
            break
        seq=seq[:-1]
        p_end-=1
    #print 'JUSTTRIM-->AfterTrim',seq,len(seq)
    seq,seqleft,seqright=logic(nucseq,seq,p_start,p_end,tm_target)
    return seq,seqleft,seqright
def percgc(seq):
    return str(100*(seq.count('g')+seq.count('c'))/float(len(seq)))

def score(seq,tm_target,flag=0):
    tm=tmcalculator(seq)
    tm_lowlim=tm_target-srange
    tm_uplim=tm_target+srange
    gc_content=((seq.count('g')+seq.count('c'))/float(len(seq)-3))
    gc_score=round((100*gc_content*(1.-gc_content)),2)
    tm_score=round(-((tm-tm_lowlim)*(tm-tm_uplim))/((srange**2)/25.),2)
    len_score=round(-(len(seq)**2-70*len(seq)+1200),2) # This scoring is optimized for primers with length in between 30-40 nucleotides.Score values: [-Inf to 25]
    if flag==0:
        wing=(len(seq.split('NNS')[0])-float(len(seq.split('NNS')[1])))
        wing_score=round((25-wing**2),2)
        totalscore= round(gc_score+tm_score+len_score+wing_score,2)
        return gc_score,tm_score,len_score,wing_score,totalscore
    else:
        totalscore=str(float(gc_score)+float(tm_score)+float(len_score))
        return gc_score,tm_score,len_score,totalscore