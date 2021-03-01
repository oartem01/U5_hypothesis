import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import seaborn as sns

def splice_junctions(file):#Not used by U5_bp_distributions.py
    """Extracts splice junctions out of txt files, that are created by copy
pasting of a gene reference sequence from LOVD database: copy starts exactly
with the line of the first Met and ends with the line of the stop codon.
Returns a list of splice junctions as 11nt strings
(8nt exon end + 3nt exon start)."""
    
    with open(file, 'r') as f:
        file=f.read().splitlines()
    print(file)
    seq_lines=file[::4]
    print(seq_lines)
    seq=''
    for line in seq_lines:
        c=line.count(' | ')
        if c==1:
            seq+=line[1:64]
        elif c==2:
            seq+=line[1:67]
        else:
            seq+=line[1:61]
    print(seq)
    exons=seq.split(' | ')
    print(exons)
    
    splice_junctions=[]
    i=1
    for exon in exons:
        if i>1:
            ex2=exon[:3]
            splice_junctions.append(ex1+ex2)
        ex1=exon[-8:]
        i+=1
    print(splice_junctions)
    return(splice_junctions)

def e_splice_junctions(file):
    """Extracts splice junctions out of fa.txt files downlouded from ensembl.
Returns a list of splice junctions as 11nt strings
(8nt exon end + 3nt exon start)."""
    
    with open(file, 'r') as f:
        file=f.read().splitlines()
    #print(file)
    j=''
    file_joined=j.join(file)
    #print(file_joined)
    file_trimmed=file_joined[1:]
    seq_lines=file_trimmed.split('>')
    #print(seq_lines)
    exons=[]
    for line in seq_lines:
        seq=line.split('coding')
        #print(seq)
        exons.append(seq[1])
    #print(exons)
    
    splice_junctions=[]
    i=1
    for exon in exons:
        if i>1:
            ex2=exon[:3]
            splice_junctions.append(ex1+ex2)
        ex1=exon[-8:]
        i+=1
    #print(splice_junctions)
    return(splice_junctions)

def bp_position(splice_junctions):
    """Pairs U5 with the nucleotide distribution by position from the/
splice junctions of a gene. Returns a dictionary of lists of/
bp type counts by position."""
    u5='UCAUUUUCCGC'
    bp_position={}
    pair=''
    GC=[0,0,0,0,0,0,0,0,0,0,0]
    AU=[0,0,0,0,0,0,0,0,0,0,0]
    GU=[0,0,0,0,0,0,0,0,0,0,0]
    UU=[0,0,0,0,0,0,0,0,0,0,0]
    AC=[0,0,0,0,0,0,0,0,0,0,0]
    CU=[0,0,0,0,0,0,0,0,0,0,0]
    CC=[0,0,0,0,0,0,0,0,0,0,0]
    GG=[0,0,0,0,0,0,0,0,0,0,0]
    AA=[0,0,0,0,0,0,0,0,0,0,0]
    AG=[0,0,0,0,0,0,0,0,0,0,0]
    for j in splice_junctions:
        n=0
        while n<11:
            
            pair=str(j[n])+str(u5[n])
            
            if pair=='GC' or pair=='CG':
                GC[n]+=1
            elif pair=='AU' or pair=='TA':
                AU[n]+=1
            elif pair=='GU' or pair=='TG':
                GU[n]+=1
            elif pair=='TU':
                UU[n]+=1
            elif pair=='AC' or pair=='CA':
                AC[n]+=1
            elif pair=='CU' or pair=='TC':
                CU[n]+=1
            elif pair=='CC':
                CC[n]+=1
            elif pair=='GG':
                GG[n]+=1
            elif pair=='AA':
                AA[n]+=1
            else:
                AG[n]+=1
            n+=1
            if n==11:
                break
            
    bp_position['GC']=GC
    bp_position['AU']=AU
    bp_position['GU']=GU
    bp_position['UU']=UU
    bp_position['AC']=AC
    bp_position['CU']=CU
    bp_position['CC']=CC
    bp_position['GG']=GG
    bp_position['AA']=AA
    bp_position['AG']=AG
    
    print(bp_position)
    return(bp_position)

def N_freq_5pr3_rand(splice_junctions, total):
    """Extracts only exon position -3 (5pr3). Returns a dictionary of/
normalised nucleotide counts - frequencies."""
    N_freq_5pr3_rand={}
    A=0
    C=0
    G=0
    U=0

    for j in splice_junctions:
        if j[5]=='A':
            A+=1
        elif j[5]=='C':
            C+=1
        elif j[5]=='G':
            G+=1
        else:
            U+=1

    freq_5pr3_A=A/total
    freq_5pr3_C=C/total
    freq_5pr3_G=G/total
    freq_5pr3_U=U/total
    
    N_freq_5pr3_rand['A']=freq_5pr3_A
    N_freq_5pr3_rand['C']=freq_5pr3_C
    N_freq_5pr3_rand['G']=freq_5pr3_G
    N_freq_5pr3_rand['U']=freq_5pr3_U    
    
    return(N_freq_5pr3_rand)


def bp_freq_position(bp_position, total):
    bp_freq_position={}

    GC_freq=[n/total for n in bp_position['GC']]
    AU_freq=[n/total for n in bp_position['AU']]
    GU_freq=[n/total for n in bp_position['GU']]
    UU_freq=[n/total for n in bp_position['UU']]
    AC_freq=[n/total for n in bp_position['AC']]
    CU_freq=[n/total for n in bp_position['CU']]
    CC_freq=[n/total for n in bp_position['CC']]
    GG_freq=[n/total for n in bp_position['GG']]
    AA_freq=[n/total for n in bp_position['AA']]
    AG_freq=[n/total for n in bp_position['AG']]

    bp_freq_position['GC']=GC_freq
    bp_freq_position['AU']=AU_freq
    bp_freq_position['GU']=GU_freq
    bp_freq_position['UU']=UU_freq
    bp_freq_position['AC']=AC_freq
    bp_freq_position['CU']=CU_freq
    bp_freq_position['CC']=CC_freq
    bp_freq_position['GG']=GG_freq
    bp_freq_position['AA']=AA_freq
    bp_freq_position['AG']=AG_freq

    print(bp_freq_position)
    return(bp_freq_position)



def geo_position(splice_junctions):
    """Pairs U5 with the nucleotide distribution by position from the/
splice junctions of a gene. Returns a dictionary of lists of/
bp geometry (Watson-Crick, isosteric, different) counts by position."""
    u5='UCAUUUUCCGC'
    geo_position={}
    pair=''
    
    Watson_Crick=[0,0,0,0,0,0,0,0,0,0,0]
    isosteric=[0,0,0,0,0,0,0,0,0,0,0]
    different=[0,0,0,0,0,0,0,0,0,0,0]

    for j in splice_junctions:
        n=0
        while n<11:
            
            pair=str(j[n])+str(u5[n])

            if pair=='GC' or pair=='CG' or pair=='AU' or pair=='TA':
                Watson_Crick[n]+=1
            elif pair=='AA' or pair=='GA' or pair=='GG' or pair=='AG' or pair=='CC':
                different[n]+=1
            else:
                isosteric[n]+=1
            n+=1
            if n==11:
                break
            
    geo_position['Watson_Crick']=Watson_Crick
    geo_position['isosteric']=isosteric
    geo_position['different']=different

    print(geo_position)
    return(geo_position)

def geo_freq_position(geo_position, total):
    geo_freq_position={}

    Watson_Crick_freq=[n/total for n in geo_position['Watson_Crick']]
    isosteric_freq=[n/total for n in geo_position['isosteric']]
    different_freq=[n/total for n in geo_position['different']]

    geo_freq_position['Watson_Crick']=Watson_Crick_freq
    geo_freq_position['isosteric']=isosteric_freq
    geo_freq_position['different']=different_freq

    print(geo_freq_position)
    return(geo_freq_position)

def bp_position_u6(ints):
    """Pairs U6 with the nucleotide distribution by position at the/
start of introns. Returns a dictionary of lists of/
bp type counts by position."""
    u6='CAUAGC'
    bp_position_u6={}
    pair=''
    GC=[0,0,0,0,0,0]
    AU=[0,0,0,0,0,0]
    GU=[0,0,0,0,0,0]
    UU=[0,0,0,0,0,0]
    AC=[0,0,0,0,0,0]
    CU=[0,0,0,0,0,0]
    CC=[0,0,0,0,0,0]
    GG=[0,0,0,0,0,0]
    AA=[0,0,0,0,0,0]
    AG=[0,0,0,0,0,0]
    for i in ints:
        n=0
        while n<6:
            
            pair=str(i[n+4])+str(u6[n])
            
            if pair=='GC' or pair=='CG':
                GC[n]+=1
            elif pair=='AU' or pair=='TA':
                AU[n]+=1
            elif pair=='GU' or pair=='TG':
                GU[n]+=1
            elif pair=='TU':
                UU[n]+=1
            elif pair=='AC' or pair=='CA':
                AC[n]+=1
            elif pair=='CU' or pair=='TC':
                CU[n]+=1
            elif pair=='CC':
                CC[n]+=1
            elif pair=='GG':
                GG[n]+=1
            elif pair=='AA':
                AA[n]+=1
            else:
                AG[n]+=1
            n+=1
            if n==6:
                break
            
    bp_position_u6['GC']=GC
    bp_position_u6['AU']=AU
    bp_position_u6['GU']=GU
    bp_position_u6['UU']=UU
    bp_position_u6['AC']=AC
    bp_position_u6['CU']=CU
    bp_position_u6['CC']=CC
    bp_position_u6['GG']=GG
    bp_position_u6['AA']=AA
    bp_position_u6['AG']=AG
    
    print(bp_position_u6)
    return(bp_position_u6)

def bp_freq_position_u6(bp_position_u6, total):
    bp_freq_position_u6={}

    GC_freq=[n/total for n in bp_position_u6['GC']]
    AU_freq=[n/total for n in bp_position_u6['AU']]
    GU_freq=[n/total for n in bp_position_u6['GU']]
    UU_freq=[n/total for n in bp_position_u6['UU']]
    AC_freq=[n/total for n in bp_position_u6['AC']]
    CU_freq=[n/total for n in bp_position_u6['CU']]
    CC_freq=[n/total for n in bp_position_u6['CC']]
    GG_freq=[n/total for n in bp_position_u6['GG']]
    AA_freq=[n/total for n in bp_position_u6['AA']]
    AG_freq=[n/total for n in bp_position_u6['AG']]

    bp_freq_position_u6['GC']=GC_freq
    bp_freq_position_u6['AU']=AU_freq
    bp_freq_position_u6['GU']=GU_freq
    bp_freq_position_u6['UU']=UU_freq
    bp_freq_position_u6['AC']=AC_freq
    bp_freq_position_u6['CU']=CU_freq
    bp_freq_position_u6['CC']=CC_freq
    bp_freq_position_u6['GG']=GG_freq
    bp_freq_position_u6['AA']=AA_freq
    bp_freq_position_u6['AG']=AG_freq

    print(bp_freq_position_u6)
    return(bp_freq_position_u6)



def geo_position_u6(ints):
    """Pairs U6 with the nucleotide distribution by position at the/
start of introns. Returns a dictionary of lists of/
bp geometry (Watson-Crick, isosteric, different) counts by position."""
    u6='CAUAGC'
    geo_position_u6={}
    pair=''
    
    Watson_Crick=[0,0,0,0,0,0]
    isosteric=[0,0,0,0,0,0]
    different=[0,0,0,0,0,0]

    for i in ints:
        n=0
        while n<6:
            
            pair=str(i[n+4])+str(u6[n])

            if pair=='GC' or pair=='CG' or pair=='AU' or pair=='TA':
                Watson_Crick[n]+=1
            elif pair=='AA' or pair=='GA' or pair=='GG' or pair=='AG' or pair=='CC':
                different[n]+=1
            else:
                isosteric[n]+=1
            n+=1
            if n==6:
                break
            
    geo_position_u6['Watson_Crick']=Watson_Crick
    geo_position_u6['isosteric']=isosteric
    geo_position_u6['different']=different

    print(geo_position_u6)
    return(geo_position_u6)

def geo_freq_position_u6(geo_position_u6, total):
    geo_freq_position_u6={}

    Watson_Crick_freq=[n/total for n in geo_position_u6['Watson_Crick']]
    isosteric_freq=[n/total for n in geo_position_u6['isosteric']]
    different_freq=[n/total for n in geo_position_u6['different']]

    geo_freq_position_u6['Watson_Crick']=Watson_Crick_freq
    geo_freq_position_u6['isosteric']=isosteric_freq
    geo_freq_position_u6['different']=different_freq

    print(geo_freq_position_u6)
    return(geo_freq_position_u6)

#Separate functions for Randoms (_rand) are only needed NOT to print out 
#the distributions 1000 times
    
def geo_position_rand(splice_junctions):
    """Pairs U5 with the nucleotide distribution by position from the/
splice junctions of a gene. Returns a dictionary of lists of/
bp geometry (Watson-Crick, isosteric, different) counts by position."""
    u5='UCAUUUUCCGC'
    geo_position={}
    pair=''
    
    Watson_Crick=[0,0,0,0,0,0,0,0,0,0,0]
    isosteric=[0,0,0,0,0,0,0,0,0,0,0]
    different=[0,0,0,0,0,0,0,0,0,0,0]

    for j in splice_junctions:
        n=0
        while n<11:
            
            pair=str(j[n])+str(u5[n])

            if pair=='GC' or pair=='CG' or pair=='AU' or pair=='TA':
                Watson_Crick[n]+=1
            elif pair=='AA' or pair=='GA' or pair=='GG' or pair=='AG' or pair=='CC':
                different[n]+=1
            else:
                isosteric[n]+=1
            n+=1
            if n==11:
                break
            
    geo_position['Watson_Crick']=Watson_Crick
    geo_position['isosteric']=isosteric
    geo_position['different']=different

    #print(geo_position)
    return(geo_position)



def geo_freq_position_rand(geo_position, total):
    geo_freq_position={}

    Watson_Crick_freq=[n/total for n in geo_position['Watson_Crick']]
    isosteric_freq=[n/total for n in geo_position['isosteric']]
    different_freq=[n/total for n in geo_position['different']]

    geo_freq_position['Watson_Crick']=Watson_Crick_freq
    geo_freq_position['isosteric']=isosteric_freq
    geo_freq_position['different']=different_freq

    #print(geo_freq_position)
    return(geo_freq_position)


def geo_position_u6_rand(ints):
    """Pairs U6 with the nucleotide distribution by position at the/
start of introns. Returns a dictionary of lists of/
bp geometry (Watson-Crick, isosteric, different) counts by position."""
    u6='CAUAGC'
    geo_position_u6={}
    pair=''
    
    Watson_Crick=[0,0,0,0,0,0]
    isosteric=[0,0,0,0,0,0]
    different=[0,0,0,0,0,0]

    for i in ints:
        n=0
        while n<6:
            
            pair=str(i[n+4])+str(u6[n])

            if pair=='GC' or pair=='CG' or pair=='AU' or pair=='TA':
                Watson_Crick[n]+=1
            elif pair=='AA' or pair=='GA' or pair=='GG' or pair=='AG' or pair=='CC':
                different[n]+=1
            else:
                isosteric[n]+=1
            n+=1
            if n==6:
                break
            
    geo_position_u6['Watson_Crick']=Watson_Crick
    geo_position_u6['isosteric']=isosteric
    geo_position_u6['different']=different

    #print(geo_position_u6)
    return(geo_position_u6)

def geo_freq_position_u6_rand(geo_position_u6, total):
    geo_freq_position_u6={}

    Watson_Crick_freq=[n/total for n in geo_position_u6['Watson_Crick']]
    isosteric_freq=[n/total for n in geo_position_u6['isosteric']]
    different_freq=[n/total for n in geo_position_u6['different']]

    geo_freq_position_u6['Watson_Crick']=Watson_Crick_freq
    geo_freq_position_u6['isosteric']=isosteric_freq
    geo_freq_position_u6['different']=different_freq

    #print(geo_freq_position_u6)
    return(geo_freq_position_u6)    

def introns(file):
    """Extracts introns out of fa.txt files downlouded from ensembl./
Returns a list of complete introns."""
    
    with open(file, 'r') as f:
        file=f.read().splitlines()
    #print(file)
    j=''
    file_joined=j.join(file)
    #print(file_joined)
    file_trimmed=file_joined[1:]
    seq_lines=file_trimmed.split('>')
    #print(seq_lines)
    introns=[]
    for line in seq_lines:
        seq=line.split('coding')
        #print(seq)
        introns.append(seq[1])
    #print(introns)
    return(introns)

def intron_starts(gene_introns):
    intron_starts=[]
    for intron in gene_introns:
        intron_start=intron[:10]
        intron_starts.append(intron_start)
    #print(intron_starts)
    return(intron_starts)

def intron_ends(gene_introns):
    intron_ends=[]
    for intron in gene_introns:
        intron_end=intron[-60:]
        intron_ends.append(intron_end)
    #print(intron_ends)
    return(intron_ends)

def sample_random(sample_sequences, total):
    """Generates randoms with replacement of the same size as the sample for bootstrapping"""
    sample_random=[]
    n=1
    while n<(total+1):
        i=random.randint(0, total-1)
        sample_random.append(sample_sequences[i])
        n+=1
        if n==(total+1):
                break
    return(sample_random)
    
def geo_freq_boots_df(geo_freq_sample1_rand, geo_freq_sample2_rand, positions):
    WC1=geo_freq_sample1_rand['Watson_Crick']
    iso1=geo_freq_sample1_rand['isosteric']
    dif1=geo_freq_sample1_rand['different']
    WC2=geo_freq_sample2_rand['Watson_Crick']
    iso2=geo_freq_sample2_rand['isosteric']
    dif2=geo_freq_sample2_rand['different']
    
    geo_freq_boots_df={}
    WC_boots_df=[]
    iso_boots_df=[]
    dif_boots_df=[]
    
    i=0
    while i<positions:
         WC_boots_df.append(WC2[i]-WC1[i])
         iso_boots_df.append(iso2[i]-iso1[i])
         dif_boots_df.append(dif2[i]-dif1[i])
         i+=1
         if i==positions:
             break
         
    geo_freq_boots_df['Watson_Crick']=WC_boots_df
    geo_freq_boots_df['isosteric']=iso_boots_df
    geo_freq_boots_df['different']=dif_boots_df
    
    return(geo_freq_boots_df)
    
def Nfreq_5pr3_boots_df(Nfreq_sample1_rand, Nfreq_sample2_rand):
    A1=Nfreq_sample1_rand['A']
    C1=Nfreq_sample1_rand['C']
    G1=Nfreq_sample1_rand['G']
    U1=Nfreq_sample1_rand['U']
    A2=Nfreq_sample2_rand['A']
    C2=Nfreq_sample2_rand['C']
    G2=Nfreq_sample2_rand['G']
    U2=Nfreq_sample2_rand['U']   
    
    Nfreq_5pr3_boots_df={}        
    
    A_boots_df=(A2-A1)
    C_boots_df=(C2-C1)
    G_boots_df=(G2-G1)
    U_boots_df=(U2-U1)
    
    
    Nfreq_5pr3_boots_df['A']=A_boots_df
    Nfreq_5pr3_boots_df['C']=C_boots_df
    Nfreq_5pr3_boots_df['G']=G_boots_df
    Nfreq_5pr3_boots_df['U']=U_boots_df
    
    return(Nfreq_5pr3_boots_df)    
    
def sj_random1(sjunctions_cons, total_sub, total_cons):
    sj_random1=[]
    n=1
    while n<(total_sub+1):
        i=random.randint(0, total_cons-1)
        if sjunctions_cons[i] not in sj_random1:
            sj_random1.append(sjunctions_cons[i])
            n+=1
        if n==(total_sub+1):
                break
    return(sj_random1)
    
def sj_random2(sjunctions_cons, total_sub, total_cons, sj_rand1):
    sj_random2=[]
    n=1
    while n<(total_sub+1):
        i=random.randint(0, total_cons-1)
        if sjunctions_cons[i] not in (sj_random2+sj_rand1):
            sj_random2.append(sjunctions_cons[i])
            n+=1
        if n==(total_sub+1):
                break
    return(sj_random2) 
    
    
#Stacked barchart of basepair types counts by position for U5/splice_junctions
def bp_stacked_barchart(bp_freq_position):
    
    n_groups=11
    scores1=bp_freq_position['GC']
    scores2=bp_freq_position['AU']
    scores3=bp_freq_position['GU']
    scores4=bp_freq_position['UU']
    scores5=bp_freq_position['AC']
    scores6=bp_freq_position['CU']
    scores7=bp_freq_position['CC']
    scores8=bp_freq_position['GG']
    scores9=bp_freq_position['AA']
    scores10=bp_freq_position['AG']
    
    index=np.arange(n_groups)
    width=0.4
    rects1=plt.bar(index, scores1, width, color='0.1')
    rects2=plt.bar(index, scores2, width, bottom=scores1, color='0.6')
    start3=tuple(i+j for i, j in zip(scores1, scores2))
    rects3=plt.bar(index, scores3, width, bottom=start3, color='0.9', hatch='/////')
    start4=tuple(i+j+k for i, j, k in zip(scores1, scores2, scores3))
    rects4=plt.bar(index, scores4, width, bottom=start4, color='0.7', hatch='..')
    start5=tuple(i+j+k+l for i, j, k, l in zip(scores1, scores2, scores3, scores4))
    rects5=plt.bar(index, scores5, width, bottom=start5, color='0.9', hatch='||||---')
    start6=tuple(i+j+k+l+m for i, j, k, l, m in zip(scores1, scores2, scores3, \
        scores4, scores5))
    rects6=plt.bar(index, scores6, width, bottom=start6, color='0.9', hatch='\\\\\\\\')
    start7=tuple(i+j+k+l+m+n for i, j, k, l, m, n in zip(scores1, scores2, scores3, \
        scores4, scores5, scores6))
    rects7=plt.bar(index, scores7, width, bottom=start7, color='0.9',hatch='----')
    start8=tuple(i+j+k+l+m+n+o for i, j, k, l, m, n, o in zip(scores1, scores2, \
        scores3, scores4, scores5, scores6, scores7))
    rects8=plt.bar(index, scores8, width, bottom=start8, color='0.9', hatch='---')
    start9=tuple(i+j+k+l+m+n+o+p for i, j, k, l, m, n, o, p in zip(scores1, \
        scores2,  scores3, scores4, scores5, scores6, scores7, scores8))
    rects9=plt.bar(index, scores9, width, bottom=start9, color='0.9',hatch='--')
    start10=tuple(i+j+k+l+m+n+o+p+q for i, j, k, l, m, n, o, p, q in zip(scores1, \
        scores2,  scores3, scores4, scores5, scores6, scores7, scores8, scores9))
    rects10=plt.bar(index, scores10, width, bottom=start10, color='0.9', hatch='-')
    
    plt.xticks(index, ('-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', '+1', '+2', '+3'))
    #legend=plt.legend(loc='below')
    
    plt.show()
    return

def bp_stacked_barchart_5pr3(bp_freq_inGsub, bp_freq_inG):
    #imports done above
    n_groups=2
    scores1=(bp_freq_inGsub['AU'][5], bp_freq_inG['AU'][5]) 
    scores2=(bp_freq_inGsub['CU'][5], bp_freq_inG['CU'][5])
    scores3=(bp_freq_inGsub['GU'][5], bp_freq_inG['GU'][5])
    scores4=(bp_freq_inGsub['UU'][5], bp_freq_inG['UU'][5])
    
    f,ax =plt.subplots(figsize=(6,6))
    index=np.arange(n_groups)    
    width=0.6
    rects1=plt.bar(index, scores1, width, color='0.1')
    rects2=plt.bar(index, scores2, width, bottom=scores1, color='0.7')
    start3=tuple(i+j for i, j in zip(scores1, scores2))
    rects3=plt.bar(index, scores3, width, bottom=start3, color='0.9', hatch='///')
    start4=tuple(i+j+k for i, j, k in zip(scores1, scores2, scores3))
    rects4=plt.bar(index, scores4, width, bottom=start4, color='0.7', hatch='..')

    plt.xticks(index, ('+5S', '+5G'))
    #legend=plt.legend(loc='below')
    plt.show()
    return

#Stacked barchart of basepair geometry counts by position for U5/splice_junctions
def geo_stacked_barchart(geo_freq_position):
    
    #sns.set(style="whitegrid")
    #f,ax =plt.subplots(figsize=(15,8.1))    
    #palette="Set1"

    n_groups=11
    scores1=geo_freq_position['Watson_Crick']
    scores2=geo_freq_position['isosteric']
    scores3=geo_freq_position['different']

    index=np.arange(n_groups)
    width=0.6
    rects1=plt.bar(index, scores1, width, color='0.1')
    rects2=plt.bar(index, scores2, width, bottom=scores1, color='0.7')
    start3=tuple(i+j for i, j in zip(scores1, scores2))
    rects3=plt.bar(index, scores3, width, bottom=start3, color='0.9', hatch='///')

    plt.xticks(index, ('-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', '+1', '+2', '+3'))
    #legend=plt.legend(loc='below')
    #sns.despine(left=True, bottom=True)
    plt.show()
    return

#Stacked barchart of basepair geometry counts by position for U6/start of intron
def geo_stacked_barchart_u6(geo_freq_position_u6):
    
    n_groups=6
    scores1=geo_freq_position_u6['Watson_Crick']
    scores2=geo_freq_position_u6['isosteric']
    scores3=geo_freq_position_u6['different']

    index=np.arange(n_groups)
    width=0.6
    rects1=plt.bar(index, scores1, width, color='0.1')
    rects2=plt.bar(index, scores2, width, bottom=scores1, color='0.7')
    start3=tuple(i+j for i, j in zip(scores1, scores2))
    rects3=plt.bar(index, scores3, width, bottom=start3, color='0.9', hatch='///')

    plt.xticks(index, ('+5', '+6', '+7', '+8', '+9', '+10'))
    #legend=plt.legend(loc='below')
    plt.show()
    return


#Symmetrised Kullback-Liebler (sKL) distance
def sKL_SJ(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    geo_1=WC+iso+dif
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    geo_2=WC+iso+dif
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_SJ=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_SJ+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_SJ)

def sKL_5prime(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=0
    while i<8:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==8:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=0
    while i<8:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==8:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_5prime=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_5prime+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_5prime)
    
def sKL_5pr876(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=0
    while i<3:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==3:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=0
    while i<3:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==3:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_5pr876=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_5pr876+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_5pr876)    
    
def sKL_5pr543(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=3
    while i<6:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==6:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=3
    while i<6:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==6:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_5pr543=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_5pr543+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_5pr543)   

def sKL_5pr21(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=6
    while i<8:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==8:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=6
    while i<8:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==8:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_5pr21=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_5pr21+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_5pr21)
    

def sKL_3prime(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=8
    while i<11:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==11:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=8
    while i<11:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==11:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_3prime=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_3prime+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_3prime)
    
def sKL_5pr1(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=7
    geo_1.append(WC[i])
    geo_1.append(iso[i])
    geo_1.append(dif[i])
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=7
    geo_2.append(WC[i])
    geo_2.append(iso[i])
    geo_2.append(dif[i])

    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_5pr1=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_5pr1+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_5pr1)

def sKL_3pr1(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=8   
    geo_1.append(WC[i])
    geo_1.append(iso[i])
    geo_1.append(dif[i])
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=8
    geo_2.append(WC[i])
    geo_2.append(iso[i])
    geo_2.append(dif[i])
   #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_3pr1=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_3pr1+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_3pr1)

def sKL_3pr23(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=9
    while i<11:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==11:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total SJ positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=9
    while i<11:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==11:
            break        
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_3pr23=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_3pr23+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_3pr23)

def mean_sKL(sKL_list):
    sKL=np.array(sKL_list, float)
    mean_sKL=sKL.mean()
    return(mean_sKL)
    
def std_dev_sKL(sKL_list):
    sKL=np.array(sKL_list, float)
    std_dev_sKL=sKL.std()
    return(std_dev_sKL)


def sKL_inU6(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=0
    while i<6:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==6:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total intron positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=0
    while i<6:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==6:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_inU6=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_inU6+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_inU6)

def sKL_in5678(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=0
    while i<4:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==4:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total intron positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=0
    while i<4:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==4:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_in5678=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_in5678+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_in5678)
    
def sKL_in56(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=0
    while i<2:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==2:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total intron positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=0
    while i<2:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==2:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_in56=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_in56+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_in56)
    
def sKL_in78(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=2
    while i<4:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==4:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total intron positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=2
    while i<4:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==4:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_in78=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_in78+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_in78)

def sKL_in910(geo_freq_1, geo_freq_2):
    WC=[]
    iso=[]
    dif=[]
    geo_1=[]
    geo_2=[]
    P=[]
    Q=[]
    epsilon=0.00001
    
    WC=geo_freq_1['Watson_Crick']
    iso=geo_freq_1['isosteric']
    dif=geo_freq_1['different']
    i=4
    while i<6:
        geo_1.append(WC[i])
        geo_1.append(iso[i])
        geo_1.append(dif[i])
        i+=1
        if i==6:
            break
    #print(geo_1)
    P=[p+epsilon for p in geo_1]
    positions=len(P)/3
    #print(P)
    #print('Total intron positions analysied: ', positions)
    
    WC=geo_freq_2['Watson_Crick']
    iso=geo_freq_2['isosteric']
    dif=geo_freq_2['different']
    i=4
    while i<6:
        geo_2.append(WC[i])
        geo_2.append(iso[i])
        geo_2.append(dif[i])
        i+=1
        if i==6:
            break
    #print(geo_2)
    Q=[q+epsilon for q in geo_2]
    #print(Q)
    
    sKL_in910=0
    i=0
    while i<len(P):
        skl=(P[i]*np.log2(P[i]/Q[i]))+(Q[i]*np.log2(Q[i]/P[i]))
        sKL_in910+=skl
        i+=1
        if i==len(P):
            break
    return(sKL_in910)

#Distribution histograms for sKL_RR and sKL_GsubR
#def hist(sKL_list_RR, sKL_list_GsubR, mean_sKL_RR, mean_sKL_GsubR, std_dev_sKL_RR, std_dev_sKL_GsubR):
def hist(sKL_list_RR, sKL_list_GsubR):
    #f,ax =plt.subplots(figsize=(20,6))
    plt.hist(sKL_list_RR, bins=50)
    plt.hist(sKL_list_GsubR, bins=50)
    
    plt.show()
    return

def hist_boots_df(BD):
    
    plt.hist(BD, bins=50)
    plt.show()
    return

def boots_confidence_interval_99(BD):
    BD.sort()
    boots_confidence_interval_99=(BD[49],BD[9949])
    return(boots_confidence_interval_99)

def boots_confidence_interval_999(BD):
    BD.sort()
    boots_confidence_interval_999=(BD[4],BD[9994])
    return(boots_confidence_interval_999)

def boots_dif_mean(BD_list):
    BD=np.array(BD_list, float)
    boots_dif_mean=BD.mean
    return(boots_dif_mean)

def boots_P_H0(BD):
    """Calculates the number of bootstrap differences that are cut off by 0 and 
    returns the probability for the null hypothesis. 
    Bootstrap differences of 0 for non-isoteric pairs of uracil that do not exist
    can be excluded in the main code."""
    
    boots_dif_mean=np.mean(BD)
    #print(boots_dif_mean)

    if boots_dif_mean==0:#Bootstrap differences of 0 for non-isoteric pairs of U.
        boots_P_H0=1#The 'distributions' are the 'same', as such pairs do not exist. 
    
    else: 
        if boots_dif_mean>0:
            BD.sort()
        else:
            BD.sort(reverse=True)
        bd_cut=0        
        for bd in BD:   
            if boots_dif_mean>0:
                if bd<0:
                    bd_cut+=1
                else:
                    break               
            else:
                if bd>0:
                    bd_cut+=1
                else:
                    break 
        
        boots_P_H0=bd_cut/10000
    
    return(boots_P_H0)
    
def violin_plot(BDgeo_positions):#(dictionary)
    df=pd.DataFrame(BDgeo_positions)
    sns.set(style="whitegrid")
    f,ax =plt.subplots(figsize=(15,8.1))#(15,8.1)(10,6))    
    sns.violinplot(data=df, palette="Set3", cut=1, linewidth=1)
    sns.despine(left=True, bottom=True)
    
    plt.show()
    return

def snshist(sKL_RR_GsubR):#(dictionary)   
    df=sKL_RR_GsubR
    
    sns.distplot(df)
    
    plt.show()
    return

def PH0_sKL(MsKL_RR, MsKL_GsubR, sKL_RR, sKL_GsubR):
    
    c_RR=0
    c_GsubR=0
    o_RR=[]
    o_GsubR=[]

    if MsKL_RR<MsKL_GsubR:
        if MsKL_GsubR/MsKL_RR<1.75:
            print('Distributions coincide')
            PH0_sKL='No difference'


        else:
            M=max(sKL_RR)
            m=min(sKL_GsubR)
            if M<m:
                print('No overlap')
                PH0_sKL='0.0000'


            
            else:
                for d in sKL_RR:
                    if d>m:
                        o_RR.append(d)
                        c_RR+=1
                for d in sKL_GsubR:
                    if d<M:
                        o_GsubR.append(d)
                        c_GsubR+=1
                print(c_RR, len(o_RR),c_GsubR, len(o_GsubR))
                o_RR.sort()
                o_GsubR.sort()
            
            
                d=100
                if c_GsubR>=c_RR:
                    F=c_GsubR/c_RR
                    #print(F)
                    i=0
                    while i<(c_RR-1):
                        D=abs(o_RR[c_RR-1-i]-o_GsubR[int(i*F)])
                        if D<d:
                            d=D
                            n=int(i*F)
                            n_RR=c_RR-1-i
                        i+=1
                        if i==c_RR:
                            break
                    N=n
                    N_RR=n_RR
                    if n>=int(F)+1: 
                        i=n-int(F)-1
                    else:
                        i=0
                        while i<n+int(F)+1:
                            if n_RR>0:
                                D1=abs(o_RR[n_RR-1]-o_GsubR[i])
                            else: 
                                D1=100
                            D2=abs(o_RR[n_RR]-o_GsubR[i])
                            D3=abs(o_RR[n_RR+1]-o_GsubR[i])
                            D=min(D1,D2,D3)
                            if D<d:
                                d=D
                                N=i
                                if d==D1:
                                    N_RR=n_RR-1
                                elif d==D3:
                                    N_RR=n_RR+1
                            i+=1
                            if i==n+int(F)+1:
                                break
                                                            
                else:#c_GsubR<c_RR:
                        F=c_RR/c_GsubR
                        #print(F)
                        i=0
                        while i<(c_GsubR-1):
                            D=abs(o_RR[(c_RR-1)-int(i*F)]-o_GsubR[i])
                            if D<d:
                                d=D
                                n=i
                                n_RR=(c_RR-1)-int(i*F)
                            i+=1
                            if i==c_GsubR:
                                    break
                        N=n
                        N_RR=n_RR
                        if n_RR>=int(F)+1: 
                            i=n_RR-int(F)-1
                        else:
                            i=0
                        while i<n_RR+int(F)+1:
                            if n>0:
                                D1=abs(o_RR[i]-o_GsubR[n-1])
                            else: 
                                D1=100
                            D2=abs(o_RR[i]-o_GsubR[n])
                            D3=abs(o_RR[i]-o_GsubR[n+1])
                            D=min(D1,D2,D3)
                            if D<d:
                                d=D
                                if d==D1:
                                    N=n-1
                                elif d==D3:
                                    N=n+1
                                N_RR=i
                            i+=1
                            if i==n+int(F)+1:
                                break
             
                PH0_sKL=N/10000
                #print('M, m, c_RR, c_GsubR, d, n, o_GsubR[n], N, o_GsubR[N], n_RR, o_RR[n_RR], N_RR, o_RR[N_RR], PH0_sKL:')
                #print(M, m, c_RR, c_GsubR, d, n, o_GsubR[n], N, o_GsubR[N], n_RR, o_RR[n_RR], N_RR, o_RR[N_RR], PH0_sKL)                  
    
    else:#MsKL_RR>MsKL_GsubR:
            print('Distributions coincide')
            PH0_sKL='H0 accepted'


    return(PH0_sKL)