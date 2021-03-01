#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 10:17:29 2020

@author: olga
"""

#splice_sites.py should be run first

import U5 as u5

print('Isolating splice junctions of major GU-AG introns with substitutions of the end of exon G') 
introns_exGsub=[]#List of indices with end of exon G substitutions
introns_exG=[]#List of indices with end of exon G

i=0
while i<total_u6:
    if splice_junctions_u6[i][7]=='G':
        introns_exG.append(intron_starts_u6[i])
    else:
        introns_exGsub.append(intron_starts_u6[i])
 
    i+=1
    if i==total_u6:
        break

      
total_exGsub=len(introns_exGsub)
total_exG=len(introns_exG)
percent_exGsub=(total_exGsub/total_u6)*100

print('Total major GU-AG introns: ', total_u6)
print('Total exons missing the end of exon G (-1G): ', total_exGsub, ', ','{:6.2f}'.format(percent_exGsub),'%')
print('Total exons with conserved end of exon G (-1G): ', total_exG, '(adds up to ',total_exGsub+total_exG,')')

print('U6 base pairs with the start of the intron, -1G substitutions')

print('U6 base pairs with the start of the intron after exons missing -1G, N=',total_exGsub)
bp_exGsub=u5.bp_position_u6(introns_exGsub)

print('Frequency of U6 base pairs with the start of the intron after exons missing -1G, N=',total_exGsub)
bp_freq_exGsub=u5.bp_freq_position_u6(bp_exGsub, total_exGsub)

print('U6 base pair geometry: with introns after exons missing -1G, N=', total_exGsub)
geo_exGsub=u5.geo_position_u6(introns_exGsub)

print('U6 base pair geometry: frequencies by position, introns after exons missing -1G, N=', total_exGsub)
geo_freq_exGsub=u5.geo_freq_position_u6(geo_exGsub, total_exGsub)

print('U6 base pairs with the start of the intron after exons with conserved -1G')

print('U6 base pairs with the start of the intron after exons with -1G, N=',total_exG)
bp_exG=u5.bp_position_u6(introns_exG)

print('Frequency of U6 base pairs with the start of the intron after exons with -1G, N=',total_exG)
bp_freq_exG=u5.bp_freq_position_u6(bp_exG, total_exG)

print('U6 base pair geometry: with the start of the intron after exons with -1G, N=', total_exG)
geo_exG=u5.geo_position_u6(introns_exG)

print('U6 base pair geometry: frequencies by position, introns after exons with -1G, N=', total_exG)
geo_freq_exG=u5.geo_freq_position_u6(geo_exG, total_exG)


print('Barcharts for bp geometry at U6 binding site at the start of the intron')

u5.geo_stacked_barchart_u6(geo_freq_exGsub)
print('U6 base pair geometry: frequencies by position, introns after exons missing -1G, N=', total_exGsub)
print('')
print('')
      
u5.geo_stacked_barchart_u6(geo_freq_exG)
print('U6 base pair geometry: frequencies by position, introns after exons with the conserved -1G, N=', total_exG)
print('')
print('')


print('Symmetrised Kullback-Liebler distance (sKL) between bp geometry \
      distributions at the start of intron')

print('sKLs for bp geometry in introns after exons missing -1G (N=', total_exGsub,') in introns after exons with -1G (N=',total_exG,'):')

sKL_inU6_exG=u5.sKL_inU6(geo_freq_exGsub, geo_freq_exG)
print('Intron positions +5 to +10:', sKL_inU6_exG)

sKL_in5678_exG=u5.sKL_in5678(geo_freq_exGsub, geo_freq_exG)
print('Intron positions +5 to +8:', sKL_in5678_exG)

sKL_in56_exG=u5.sKL_in56(geo_freq_exGsub, geo_freq_exG)
print('Intron positions +5 to +6:', sKL_in56_exG)
      
sKL_in78_exG=u5.sKL_in78(geo_freq_exGsub, geo_freq_exG)
print('Intron positions +7 to +8:', sKL_in78_exG)

sKL_in910_exG=u5.sKL_in910(geo_freq_exGsub, geo_freq_exG)
print('Intron positions +9 to +10:', sKL_in910_exG)


print('Random selections of pairs of distributions (x10000) of ', total_exGsub, 'introns out of',total_exG, 'preceded by exon end G') 
print('Generate lists of 10000 sKLs for Random1/Random2 and exGsub/Random1' )
print('Random2 does not share splice junctions with Random1!')

print('sKL distributions for the start of intron U6 site')
print('Influence of the end of exon G')

sKL_inU6_RR=[]
sKL_inU6_exGsubR=[]
sKL_in5678_RR=[]
sKL_in5678_exGsubR=[]
sKL_in56_RR=[]
sKL_in56_exGsubR=[]
sKL_in78_RR=[]
sKL_in78_exGsubR=[]
sKL_in910_RR=[]
sKL_in910_exGsubR=[]

i=0
while i<10000:
    #Using the same functions u5.sj_random1 and u5.sj_random2
    #Using geo_..._rand functions NOT to print out 1000 distribution 
    #Random 1     
    introns_exG_rand1=u5.sj_random1(introns_exG, total_exGsub, total_exG)
    geo_exG_rand1=u5.geo_position_u6_rand(introns_exG_rand1)
    geo_freq_exG_rand1=u5.geo_freq_position_u6_rand(geo_exG_rand1, total_exGsub)
    #Random 2
    introns_exG_rand2=u5.sj_random2(introns_exG, total_exGsub, total_exG, introns_exG_rand1)
    geo_exG_rand2=u5.geo_position_u6_rand(introns_exG_rand2)
    geo_freq_exG_rand2=u5.geo_freq_position_u6_rand(geo_exG_rand2, total_exGsub)
    #sKLs for Random 1/Random 2 and exGsub/Random1
    sKL_inU6_R1R2=u5.sKL_inU6(geo_freq_exG_rand1, geo_freq_exG_rand2)
    sKL_inU6_RR.append(sKL_inU6_R1R2)
    sKL_inU6_exGsubR1=u5.sKL_inU6(geo_freq_exGsub, geo_freq_exG_rand1)
    sKL_inU6_exGsubR.append(sKL_inU6_exGsubR1)
    
    sKL_in5678_R1R2=u5.sKL_in5678(geo_freq_exG_rand1, geo_freq_exG_rand2)
    sKL_in5678_RR.append(sKL_in5678_R1R2)
    sKL_in5678_exGsubR1=u5.sKL_in5678(geo_freq_exGsub, geo_freq_exG_rand1)
    sKL_in5678_exGsubR.append(sKL_in5678_exGsubR1)

    sKL_in56_R1R2=u5.sKL_in56(geo_freq_exG_rand1, geo_freq_exG_rand2)
    sKL_in56_RR.append(sKL_in56_R1R2)
    sKL_in56_exGsubR1=u5.sKL_in56(geo_freq_exGsub, geo_freq_exG_rand1)
    sKL_in56_exGsubR.append(sKL_in56_exGsubR1)

    sKL_in78_R1R2=u5.sKL_in78(geo_freq_exG_rand1, geo_freq_exG_rand2)
    sKL_in78_RR.append(sKL_in78_R1R2)
    sKL_in78_exGsubR1=u5.sKL_in78(geo_freq_exGsub, geo_freq_exG_rand1)
    sKL_in78_exGsubR.append(sKL_in78_exGsubR1)

    sKL_in910_R1R2=u5.sKL_in910(geo_freq_exG_rand1, geo_freq_exG_rand2)
    sKL_in910_RR.append(sKL_in910_R1R2)
    sKL_in910_exGsubR1=u5.sKL_in910(geo_freq_exGsub, geo_freq_exG_rand1)
    sKL_in910_exGsubR.append(sKL_in910_exGsubR1)    

    i+=1
    if i==10000:
        break
    
print('Mean sKLs and strd.dev. for Random1/Random2 and exGsub/Random1')

MsKL_inU6_RR=u5.mean_sKL(sKL_inU6_RR)
MsKL_inU6_exGsubR=u5.mean_sKL(sKL_inU6_exGsubR)
MsKL_in5678_RR=u5.mean_sKL(sKL_in5678_RR)
MsKL_in5678_exGsubR=u5.mean_sKL(sKL_in5678_exGsubR)
MsKL_in56_RR=u5.mean_sKL(sKL_in56_RR)
MsKL_in56_exGsubR=u5.mean_sKL(sKL_in56_exGsubR)
MsKL_in78_RR=u5.mean_sKL(sKL_in78_RR)
MsKL_in78_exGsubR=u5.mean_sKL(sKL_in78_exGsubR)
MsKL_in910_RR=u5.mean_sKL(sKL_in910_RR)
MsKL_in910_exGsubR=u5.mean_sKL(sKL_in910_exGsubR)

StDsKL_inU6_RR=u5.std_dev_sKL(sKL_inU6_RR)
StDsKL_inU6_exGsubR=u5.std_dev_sKL(sKL_inU6_exGsubR)
StDsKL_in5678_RR=u5.std_dev_sKL(sKL_in5678_RR)
StDsKL_in5678_exGsubR=u5.std_dev_sKL(sKL_in5678_exGsubR)
StDsKL_in56_RR=u5.std_dev_sKL(sKL_in56_RR)
StDsKL_in56_exGsubR=u5.std_dev_sKL(sKL_in56_exGsubR)
StDsKL_in78_RR=u5.std_dev_sKL(sKL_in78_RR)
StDsKL_in78_exGsubR=u5.std_dev_sKL(sKL_in78_exGsubR)
StDsKL_in910_RR=u5.std_dev_sKL(sKL_in910_RR)
StDsKL_in910_exGsubR=u5.std_dev_sKL(sKL_in910_exGsubR)


print('Positions_RR/GsubR   MEAN sKL          Std Dev')
print('inU6_RR             ', MsKL_inU6_RR, StDsKL_inU6_RR)
print('inU6_exGsubR          ', MsKL_inU6_exGsubR, StDsKL_inU6_exGsubR)
print('in5678_RR         ', MsKL_in5678_RR, StDsKL_in5678_RR)
print('in5678_exGsubR      ', MsKL_in5678_exGsubR, StDsKL_in5678_exGsubR)
print('in56_RR         ', MsKL_in56_RR, StDsKL_in56_RR)
print('in56_exGsubR      ', MsKL_in56_exGsubR, StDsKL_in56_exGsubR)
print('in78_RR         ', MsKL_in78_RR, StDsKL_in78_RR)
print('in78_exGsubR      ', MsKL_in78_exGsubR, StDsKL_in78_exGsubR)
print('in910_RR          ', MsKL_in910_RR, StDsKL_in910_RR)
print('in910_exGsubR       ', MsKL_in910_exGsubR, StDsKL_in910_exGsubR)



print('Plots for sKLs')

u5.hist(sKL_inU6_RR, sKL_inU6_exGsubR)
print('sKL_inU6_RR and sKL_inU6_exGsubR')

u5.hist(sKL_in5678_RR, sKL_in5678_exGsubR)
print('sKL_in5678_RR and sKL_in5678_exGsubR')

u5.hist(sKL_in56_RR, sKL_in56_exGsubR)
print('sKL_in56_RR and sKL_in56_exGsubR')

u5.hist(sKL_in78_RR, sKL_in78_exGsubR)
print('sKL_in78_RR and sKL_in78_exGsubR')

u5.hist(sKL_in910_RR, sKL_in910_exGsubR)
print('sKL_in910_RR and sKL_in910_exGsubR')
