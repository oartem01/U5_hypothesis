#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 09:49:59 2020

@author: olga
"""
#splice_sites.py should be run first

import U5 as u5

print('End of exon and start of intron interactions with U5 and U6 snRNAs')
print('The effect of +5G substitutions on the U5 interactions with the exons')
print('Identifying intons with +5G substitutions')
G_sub=[]#List of indeces for introns with +5G substitutions
i=0
for intron in intron_starts_u6:
    #print(i+1, intron)
    if intron[4]!='G':
        G_sub.append(i)
    i+=1
    #print(intron[7], intron[9])
total_G_sub=len(G_sub)
percent_G_sub=((total_G_sub/total_u6)*100)
#print(G_sub)
print('Introns with substitutions of +5G total: ',total_G_sub, '{:6.2f}'.format(percent_G_sub),'%')
#print(G_sub[0])

sjunctions_inG_sub=[]#List of splice junctions with introns missing +5G
sjunctions_inG=[]#List of splice junctions with +5G
#print(splice_junctions[0])
i=0
for junction in splice_junctions_u6:
    if i in G_sub:
        sjunctions_inG_sub.append(splice_junctions_u6[i])
    else:
        sjunctions_inG.append(splice_junctions_u6[i])       
    i+=1
total_inG_sub=len(sjunctions_inG_sub)
#print(sjunctions_inG_sub)
print('Substitutions of +5G, total splice junctions: ',total_inG_sub)
total_inG=len(sjunctions_inG)
print('+5G, total splice junctions: ',total_inG)

print('N=',total_inG_sub,', +5G absent: U5 base pairs')
print('U5 Loop1 bp listed by position for splice junctions with introns missing +5G')
bp_inGsub=u5.bp_position(sjunctions_inG_sub)

print('U5 Loop1 bp frequencies listed by position for splice junctions with introns missing +5G')
bp_freq_inGsub=u5.bp_freq_position(bp_inGsub, total_inG_sub)
#u5.bp_stacked_barchart(bp_freq_position)

print('N=',total_inG,',+5G present: U5 base pairs')
print('U5 Loop1 bp listed by position for all splice junctions with +5G introns')
bp_inG=u5.bp_position(sjunctions_inG)

print('U5 Loop1 bp frequencies listed by position for all splice junctions with +5G introns')
bp_freq_inG=u5.bp_freq_position(bp_inG, total_inG)

print('N=',total_inG_sub,',+5G absent: U5 base pair geometry')
print('U5 Loop1 bp geometry by position for splice junctions with introns missing +5G')
geo_inGsub=u5.geo_position(sjunctions_inG_sub)

print('U5 Loop1 bp geometry: frequency by position for splice junctions with introns missing +5G')
geo_freq_inGsub=u5.geo_freq_position(geo_inGsub, total_G_sub)
u5.geo_stacked_barchart(geo_freq_inGsub)

print('N=',total_inG,',+5G present: U5 base pair geometry')
print('U5 Loop1 bp geometry by position for all splice junctions with +5G introns')
geo_inG=u5.geo_position(sjunctions_inG)

print('U5 Loop1 bp geometry: frequency by position for all splice junctions with +5G introns')
geo_freq_inG=u5.geo_freq_position(geo_inG, total_inG)
u5.geo_stacked_barchart(geo_freq_inG)


print('Symmetrised Kullback-Liebler divergence (sKL) between U5 bp geometry \
      distributions at splice junctions (SJs)')
#total_inG=445
print('Whole 11nt junctions with both exons')
sKL_SJ_inG=u5.sKL_SJ(geo_freq_inGsub, geo_freq_inG)
print('sKL for bp geometry in SJs of introns missing +5G (N=', total_inG_sub, total_G_sub, ') \
      and SJs of +5G introns (N=',total_inG,'): ', sKL_SJ_inG)
print('5prime exons only: SJ positions -8 to -1')
sKL_5prime_inG=u5.sKL_5prime(geo_freq_inGsub, geo_freq_inG)
print('sKL for bp geometry in SJs of introns missing +5G (N=', total_inG_sub, total_G_sub, ') \
      and SJs of +5G introns (N=',total_inG,'): ', sKL_5prime_inG)
print('5prime exons positions -8 to -6')
sKL_5pr876_inG=u5.sKL_5pr876(geo_freq_inGsub, geo_freq_inG)
print('sKL for bp geometry in SJs of introns missing +5G (N=', total_inG_sub, total_G_sub, ') \
      and SJs of +5G introns (N=',total_inG,'): ', sKL_5pr876_inG)
print('5prime exons positions -5 to -3')
sKL_5pr543_inG=u5.sKL_5pr543(geo_freq_inGsub, geo_freq_inG)
print('sKL for bp geometry in SJs of introns missing +5G (N=', total_inG_sub, total_G_sub, ') \
      and SJs of +5G introns (N=',total_inG,'): ', sKL_5pr543_inG)
print('5prime exons positions -2 to -1')
sKL_5pr21_inG=u5.sKL_5pr21(geo_freq_inGsub, geo_freq_inG)
print('sKL for bp geometry in SJs of introns missing +5G (N=', total_inG_sub, total_G_sub, ') \
      and SJs of +5G introns (N=',total_inG,'): ', sKL_5pr21_inG)
print('3prime exons positions +1 to +3')
sKL_3prime_inG=u5.sKL_3prime(geo_freq_inGsub, geo_freq_inG)
print('sKL for bp geometry in SJs of introns missing +5G (N=', total_inG_sub, total_G_sub, ') \
      and SJs of +5G introns (N=',total_inG,'): ', sKL_3prime_inG)
print('SKL total=', sKL_3prime_inG+2*sKL_5prime_inG-sKL_5pr876_inG-sKL_5pr543_inG-sKL_5pr21_inG)

print('Random selections of pairs of distributions (x10000) of ', total_inG_sub, ' splice junctions out of ', total_inG, ' with +5G introns') 
print('Generate lists of 10000 sKLs for Random1/Random2 and inGsub/Random1' )
print('Random2 does not share splice junctions with Random1!')

sKL_SJ_RR=[]
sKL_SJ_GsubR=[]
sKL_5prime_RR=[]
sKL_5prime_GsubR=[]
sKL_5pr876_RR=[]
sKL_5pr876_GsubR=[]
sKL_5pr543_RR=[]
sKL_5pr543_GsubR=[]
sKL_5pr21_RR=[]
sKL_5pr21_GsubR=[]
sKL_3prime_RR=[]
sKL_3prime_GsubR=[]

i=0
while i<10000:
    #Random 1    
    sjunctions_inG_rand1=u5.sj_random1(sjunctions_inG, total_inG_sub, total_inG)
    geo_inG_rand1=u5.geo_position_rand(sjunctions_inG_rand1)
    geo_freq_inG_rand1=u5.geo_freq_position_rand(geo_inG_rand1, total_inG_sub)
    #Random 2
    sjunctions_inG_rand2=u5.sj_random2(sjunctions_inG, total_inG_sub, total_inG, sjunctions_inG_rand1)
    geo_inG_rand2=u5.geo_position_rand(sjunctions_inG_rand2)
    geo_freq_inG_rand2=u5.geo_freq_position_rand(geo_inG_rand2, total_inG_sub)
    #sKLs for Random 1/Random 2 and inGsub/Random1
    sKL_SJ_R1R2=u5.sKL_SJ(geo_freq_inG_rand1, geo_freq_inG_rand2)
    sKL_SJ_RR.append(sKL_SJ_R1R2)
    sKL_SJ_GsubR1=u5.sKL_SJ(geo_freq_inGsub, geo_freq_inG_rand1)
    sKL_SJ_GsubR.append(sKL_SJ_GsubR1)
    
    sKL_5prime_R1R2=u5.sKL_5prime(geo_freq_inG_rand1, geo_freq_inG_rand2)
    sKL_5prime_RR.append(sKL_5prime_R1R2)
    sKL_5prime_GsubR1=u5.sKL_5prime(geo_freq_inGsub, geo_freq_inG_rand1)
    sKL_5prime_GsubR.append(sKL_5prime_GsubR1)
    
    sKL_5pr876_R1R2=u5.sKL_5pr876(geo_freq_inG_rand1, geo_freq_inG_rand2)
    sKL_5pr876_RR.append(sKL_5pr876_R1R2)
    sKL_5pr876_GsubR1=u5.sKL_5pr876(geo_freq_inGsub, geo_freq_inG_rand1)
    sKL_5pr876_GsubR.append(sKL_5pr876_GsubR1)

    sKL_5pr543_R1R2=u5.sKL_5pr543(geo_freq_inG_rand1, geo_freq_inG_rand2)
    sKL_5pr543_RR.append(sKL_5pr543_R1R2)
    sKL_5pr543_GsubR1=u5.sKL_5pr543(geo_freq_inGsub, geo_freq_inG_rand1)
    sKL_5pr543_GsubR.append(sKL_5pr543_GsubR1)

    sKL_5pr21_R1R2=u5.sKL_5pr21(geo_freq_inG_rand1, geo_freq_inG_rand2)
    sKL_5pr21_RR.append(sKL_5pr21_R1R2)
    sKL_5pr21_GsubR1=u5.sKL_5pr21(geo_freq_inGsub, geo_freq_inG_rand1)
    sKL_5pr21_GsubR.append(sKL_5pr21_GsubR1)
    
    sKL_3prime_R1R2=u5.sKL_3prime(geo_freq_inG_rand1, geo_freq_inG_rand2)
    sKL_3prime_RR.append(sKL_3prime_R1R2)
    sKL_3prime_GsubR1=u5.sKL_3prime(geo_freq_inGsub, geo_freq_inG_rand1)
    sKL_3prime_GsubR.append(sKL_3prime_GsubR1)

    i+=1
    if i==10000:
        break
    
print('Mean sKLs and strd.dev. for Random1/Random2 and inGsub/Random1')
print('Check sKL lists')

print('sKL_SJ_RR list is ',len(sKL_SJ_RR), ' long.')
print('sKL_SJ_GsubR list is ',len(sKL_SJ_GsubR), ' long.')
print('sKL_5prime_RR list is ',len(sKL_5prime_RR), ' long.')
print('sKL_5prime_GsubR list is ',len(sKL_5prime_GsubR), ' long.')
print('sKL_5pr876_RR list is ',len(sKL_5pr876_RR), ' long.')
print('sKL_5pr876_GsubR list is ',len(sKL_5pr876_GsubR), ' long.')
print('sKL_5pr543_RR list is ',len(sKL_5pr543_RR), ' long.')
print('sKL_5pr543_GsubR list is ',len(sKL_5pr543_GsubR), ' long.')
print('sKL_5pr21_RR list is ',len(sKL_5pr21_RR), ' long.')
print('sKL_5pr21_GsubR list is ',len(sKL_5pr21_GsubR), ' long.')
print('sKL_3prime_RR list is ',len(sKL_3prime_RR), ' long.')
print('sKL_3prime_GsubR list is ',len(sKL_3prime_GsubR), ' long.')

MsKL_SJ_RR=u5.mean_sKL(sKL_SJ_RR)
MsKL_SJ_GsubR=u5.mean_sKL(sKL_SJ_GsubR)
MsKL_5prime_RR=u5.mean_sKL(sKL_5prime_RR)
MsKL_5prime_GsubR=u5.mean_sKL(sKL_5prime_GsubR)
MsKL_5pr876_RR=u5.mean_sKL(sKL_5pr876_RR)
MsKL_5pr876_GsubR=u5.mean_sKL(sKL_5pr876_GsubR)
MsKL_5pr543_RR=u5.mean_sKL(sKL_5pr543_RR)
MsKL_5pr543_GsubR=u5.mean_sKL(sKL_5pr543_GsubR)
MsKL_5pr21_RR=u5.mean_sKL(sKL_5pr21_RR)
MsKL_5pr21_GsubR=u5.mean_sKL(sKL_5pr21_GsubR)
MsKL_3prime_RR=u5.mean_sKL(sKL_3prime_RR)
MsKL_3prime_GsubR=u5.mean_sKL(sKL_3prime_GsubR)

StDsKL_SJ_RR=u5.std_dev_sKL(sKL_SJ_RR)
StDsKL_SJ_GsubR=u5.std_dev_sKL(sKL_SJ_GsubR)
StDsKL_5prime_RR=u5.std_dev_sKL(sKL_5prime_RR)
StDsKL_5prime_GsubR=u5.std_dev_sKL(sKL_5prime_GsubR)
StDsKL_5pr876_RR=u5.std_dev_sKL(sKL_5pr876_RR)
StDsKL_5pr876_GsubR=u5.std_dev_sKL(sKL_5pr876_GsubR)
StDsKL_5pr543_RR=u5.std_dev_sKL(sKL_5pr543_RR)
StDsKL_5pr543_GsubR=u5.std_dev_sKL(sKL_5pr543_GsubR)
StDsKL_5pr21_RR=u5.std_dev_sKL(sKL_5pr21_RR)
StDsKL_5pr21_GsubR=u5.std_dev_sKL(sKL_5pr21_GsubR)
StDsKL_3prime_RR=u5.std_dev_sKL(sKL_3prime_RR)
StDsKL_3prime_GsubR=u5.std_dev_sKL(sKL_3prime_GsubR)

print('Positions_RR/GsubR   MEAN sKL          Std Dev')
print('SJ_RR             ', MsKL_SJ_RR, StDsKL_SJ_RR)
print('SJ_GsubR          ', MsKL_SJ_GsubR, StDsKL_SJ_GsubR)
print('5prime_RR         ', MsKL_5prime_RR, StDsKL_5prime_RR)
print('5prime_GsubR      ', MsKL_5prime_GsubR, StDsKL_5prime_GsubR)
print('5pr876_RR         ', MsKL_5pr876_RR, StDsKL_5pr876_RR)
print('5pr876_GsubR      ', MsKL_5pr876_GsubR, StDsKL_5pr876_GsubR)
print('5pr543_RR         ', MsKL_5pr543_RR, StDsKL_5pr543_RR)
print('5pr543_GsubR      ', MsKL_5pr543_GsubR, StDsKL_5pr543_GsubR)
print('5pr21_RR          ', MsKL_5pr21_RR, StDsKL_5pr21_RR)
print('5pr21_GsubR       ', MsKL_5pr21_GsubR, StDsKL_5pr21_GsubR)
print('3prime_RR         ', MsKL_3prime_RR, StDsKL_3prime_RR)
print('3prime_GsubR      ', MsKL_3prime_GsubR, StDsKL_3prime_GsubR)



print('Plots for sKLs')
sKL_SJ_RR_GsubR={'sKL_SJ_RR':sKL_SJ_RR, 'sKL_SJ_GsubR':sKL_SJ_GsubR}
#u5.snshist(sKL_SJ_RR_GsubR)
#u5.hist(sKL_SJ_RR, sKL_SJ_GsubR, MsKL_SJ_RR, MsKL_SJ_GsubR, StDsKL_SJ_RR, StDsKL_SJ_GsubR)
u5.hist(sKL_SJ_RR, sKL_SJ_GsubR)
print('sKL_SJ_RR and sKL_SJ_GsubR')


u5.hist(sKL_5prime_RR, sKL_5prime_GsubR)
print('sKL_5prime_RR and sKL_5prime_GsubR')


u5.hist(sKL_5pr876_RR, sKL_5pr876_GsubR)
print('sKL_5pr876_RR and sKL_5pr876_GsubR')


u5.hist(sKL_5pr543_RR, sKL_5pr543_GsubR)
print('sKL_5pr543_RR and sKL_5pr543_GsubR')


u5.hist(sKL_5pr21_RR, sKL_5pr21_GsubR)
print('sKL_5pr21_RR and sKL_5pr21_GsubR')


u5.hist(sKL_3prime_RR, sKL_3prime_GsubR)
print('sKL_3prime_RR and sKL_3prime_GsubR')



print('Isolating the exon/intron boundaries missing both conserved Gs (end of exon G and +5G at the stat of te intron):')
n=0
for i in G_sub:
    if splice_junctions_u6[i][7]!='G':
        print(i+1, splice_junctions_u6[i], intron_starts_u6[i], intron_ends_u6[i])
        n+=1
total_noGs=n
percent_noGs=((n/total_u6)*100)
print("Total exon/intron boundaries missing both conserved Gs: ",total_noGs, '{:6.2f}'.format(percent_noGs),"%")