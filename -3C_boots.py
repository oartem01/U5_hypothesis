#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 23:48:00 2021

@author: olga
"""

import U5 as u5

print('Start of exon and end of intron interactions with U5 and U2(?) snRNAs')      
print('Isolating introns with -3C and missing -3C (major, U2 introns only (total: ', total_u2, ')')
c_sub=[]#List of indeces for introns missing -3C
c=[]#List of indeces for introns with -3C
i=0
for intron in intron_ends_u2:
    if intron[-3]=='C':
        c.append(i)
    else:
        c_sub.append(i)
    i+=1
total_c_sub=len(c_sub)
total_c=len(c)

percent_c_sub=(total_c_sub/total_u2)*100
percent_c=(total_c/total_u2)*100

print('Introns with -3C: ', total_c, '{:6.2f}'.format(percent_c), '%')
print('introns missing -3C: ', total_c_sub, '{:6.2f}'.format(percent_c_sub), '%')

print('Extracting splice junctions of the introns with -3C and missing -3C (major, U2 introns only')
sjunctions_c=[]#List of splice junctions with introns with -3C
sjunctions_c_sub=[]#List of splice junctions with introns missing -3C

i=0
for junction in splice_junctions_u2:
    if i in c:
        sjunctions_c.append(splice_junctions_u2[i])
    else:
        sjunctions_c_sub.append(splice_junctions_u2[i])
    i+=1
total_sj_c=len(sjunctions_c)
total_sj_c_sub=len(sjunctions_c_sub)

print('Splice junctions of -3C introns: ',total_sj_c,', -3subC introns: ',total_sj_c_sub)
print('Check: these numbers should equal the numbers above')

print('In relation to intron position -3: Distribution of base pairs at splice junctions')  

print('U5 base pairs with splice junctions of -3C introns, N=',total_c)
bp_c=u5.bp_position(sjunctions_c)

print('Frequency of U5 Loop1 base pairs for -3C introns, N=',total_c)
bp_c_freq=u5.bp_freq_position(bp_c, total_c)

print('U5 base pair geometry: splice junctions of -3C introns, N=', total_c)
geo_c=u5.geo_position(sjunctions_c)

print('U5 base pair geometry: frequencies for -3C intons, N=',total_c)
geo_c_freq=u5.geo_freq_position(geo_c, total_c)


print('U5 base pairs with splice junctions of introns missing -3C, N=', total_c_sub)
bp_c_sub=u5.bp_position(sjunctions_c_sub)

print('Frequency of U5 Loop1 base pairs for introns missing -3C, N=',total_c_sub)
bp_c_sub_freq=u5.bp_freq_position(bp_c_sub, total_c_sub)

print('U5 base pair geometry: splice junctions of introns missing -3C, N=', total_c_sub)
geo_c_sub=u5.geo_position(sjunctions_c_sub)

print('U5 base pair geometry: frequencies for intons missing -3C, N=',total_c_sub)
geo_c_sub_freq=u5.geo_freq_position(geo_c_sub, total_c_sub)


print('Stacked barcharts for U5 bp geometry frequencies by position')
u5.geo_stacked_barchart(geo_c_freq)
print('U5 base pair geometry: frequencies for -3C intons, N=',total_c)
print(' ')
u5.geo_stacked_barchart(geo_c_sub_freq)
print('U5 base pair geometry: frequencies for intons missing -3C, N=',total_c_sub)

print('Hypothesis testing: difference of U5/exons bp geometry frequencies by bootstrapping of -3C and -3Csub')
print('33 categories: 3 types of bp geometry in 11 splice junction positions')
print('33 lists of 10000 bootstrap differencies (BD) of frequencies each for bp geometries at positions of splice junctions for introns with conserved -3C (N=',total_c,'), and lacking -3C (N=', total_c_sub,')')

BD_5pr8_WC=[]
BD_5pr8_iso=[]
BD_5pr8_dif=[]
BD_5pr7_WC=[]
BD_5pr7_iso=[]
BD_5pr7_dif=[]
BD_5pr6_WC=[]
BD_5pr6_iso=[]
BD_5pr6_dif=[]
BD_5pr5_WC=[]
BD_5pr5_iso=[]
BD_5pr5_dif=[]
BD_5pr4_WC=[]
BD_5pr4_iso=[]
BD_5pr4_dif=[]
BD_5pr3_WC=[]
BD_5pr3_iso=[]
BD_5pr3_dif=[]
BD_5pr2_WC=[]
BD_5pr2_iso=[]
BD_5pr2_dif=[]
BD_5pr1_WC=[]
BD_5pr1_iso=[]
BD_5pr1_dif=[]
BD_3pr1_WC=[]
BD_3pr1_iso=[]
BD_3pr1_dif=[]
BD_3pr2_WC=[]
BD_3pr2_iso=[]
BD_3pr2_dif=[]
BD_3pr3_WC=[]
BD_3pr3_iso=[]
BD_3pr3_dif=[]

i=0
while i<10000:
    
    sjunctions_c_rand=u5.sample_random(sjunctions_c, total_c)
    geo_c_rand=u5.geo_position_rand(sjunctions_c_rand)
    geo_freq_c_rand=u5.geo_freq_position_rand(geo_c_rand, total_c)
    if i==0:
        print('geo_freq_c_rand', geo_freq_c_rand)
    
    sjunctions_c_sub_rand=u5.sample_random(sjunctions_c_sub, total_c_sub)
    geo_c_sub_rand=u5.geo_position_rand(sjunctions_c_sub_rand)
    geo_freq_c_sub_rand=u5.geo_freq_position_rand(geo_c_sub_rand, total_c_sub)
    if i==0:
        print('geo_freq_c_sub_rand', geo_freq_c_sub_rand)
    #Bootstrap difference in bp geometry frequencies between c and ic_sub
    #Temporariry dictionary is returned by the function and re-written 10000 times with each iteration
    positions=11
    geo_freq_boots_df=u5.geo_freq_boots_df(geo_freq_c_rand, geo_freq_c_sub_rand, positions)
    if i==0:
        print('geo_freq_boots_df', geo_freq_boots_df)    
    #Bootstrap differences 
    #33 lists to append to: 3 types of bp geometry in 11 splice junction positions
    BD_5pr8_WC.append(geo_freq_boots_df['Watson_Crick'][0])
    BD_5pr8_iso.append(geo_freq_boots_df['isosteric'][0])    
    BD_5pr8_dif.append(geo_freq_boots_df['different'][0])
    BD_5pr7_WC.append(geo_freq_boots_df['Watson_Crick'][1])
    BD_5pr7_iso.append(geo_freq_boots_df['isosteric'][1])    
    BD_5pr7_dif.append(geo_freq_boots_df['different'][1])
    BD_5pr6_WC.append(geo_freq_boots_df['Watson_Crick'][2])
    BD_5pr6_iso.append(geo_freq_boots_df['isosteric'][2])    
    BD_5pr6_dif.append(geo_freq_boots_df['different'][2]) 
    BD_5pr5_WC.append(geo_freq_boots_df['Watson_Crick'][3])
    BD_5pr5_iso.append(geo_freq_boots_df['isosteric'][3])    
    BD_5pr5_dif.append(geo_freq_boots_df['different'][3])
    BD_5pr4_WC.append(geo_freq_boots_df['Watson_Crick'][4])
    BD_5pr4_iso.append(geo_freq_boots_df['isosteric'][4])    
    BD_5pr4_dif.append(geo_freq_boots_df['different'][4])
    BD_5pr3_WC.append(geo_freq_boots_df['Watson_Crick'][5])
    BD_5pr3_iso.append(geo_freq_boots_df['isosteric'][5])    
    BD_5pr3_dif.append(geo_freq_boots_df['different'][5])
    BD_5pr2_WC.append(geo_freq_boots_df['Watson_Crick'][6])
    BD_5pr2_iso.append(geo_freq_boots_df['isosteric'][6])    
    BD_5pr2_dif.append(geo_freq_boots_df['different'][6])
    BD_5pr1_WC.append(geo_freq_boots_df['Watson_Crick'][7])
    BD_5pr1_iso.append(geo_freq_boots_df['isosteric'][7])    
    BD_5pr1_dif.append(geo_freq_boots_df['different'][7])
    BD_3pr1_WC.append(geo_freq_boots_df['Watson_Crick'][8])
    BD_3pr1_iso.append(geo_freq_boots_df['isosteric'][8])    
    BD_3pr1_dif.append(geo_freq_boots_df['different'][8])
    BD_3pr2_WC.append(geo_freq_boots_df['Watson_Crick'][9])
    BD_3pr2_iso.append(geo_freq_boots_df['isosteric'][9])    
    BD_3pr2_dif.append(geo_freq_boots_df['different'][9])    
    BD_3pr3_WC.append(geo_freq_boots_df['Watson_Crick'][10])
    BD_3pr3_iso.append(geo_freq_boots_df['isosteric'][10])    
    BD_3pr3_dif.append(geo_freq_boots_df['different'][10])    
        
    i+=1
    if i==10000:
        break

print('Violin plots: bootstrap difference of U5 bp geometry between -3Csub and -3C introns')
#Dictionaries of BD lists to create dfs for violin plots
BD_WC={'BD_5pr8_WC':BD_5pr8_WC, 'BD_5pr7_WC':BD_5pr7_WC, 
       'BD_5pr6_WC':BD_5pr6_WC, 'BD_5pr5_WC':BD_5pr5_WC, 
       'BD_5pr4_WC':BD_5pr4_WC, 'BD_5pr3_WC':BD_5pr3_WC,
       'BD_5pr2_WC':BD_5pr2_WC, 'BD_5pr1_WC':BD_5pr1_WC, 
       'BD_3pr1_WC':BD_3pr1_WC, 'BD_3pr2_WC':BD_3pr2_WC, 
       'BD_3pr3_WC':BD_3pr3_WC}
u5.violin_plot(BD_WC)
print('Watson-Crick pairs')

BD_iso={'BD_5pr8_iso':BD_5pr8_iso, 'BD_5pr7_iso':BD_5pr7_iso, 
        'BD_5pr6_iso':BD_5pr6_iso, 'BD_5pr5_iso':BD_5pr5_iso, 
        'BD_5pr4_iso':BD_5pr4_iso, 'BD_5pr3_iso':BD_5pr3_iso,
        'BD_5pr2_iso':BD_5pr2_iso, 'BD_5pr1_iso':BD_5pr1_iso, 
        'BD_3pr1_iso':BD_3pr1_iso, 'BD_3pr2_iso':BD_3pr2_iso, 
        'BD_3pr3_iso':BD_3pr3_iso}
u5.violin_plot(BD_iso)
print('Isosteric pairs (G--U/U--G, U--U, A--C/C--A and C--U/U--C)')

BD_dif={'BD_5pr8_dif':BD_5pr8_dif, 'BD_5pr7_dif':BD_5pr7_dif, 
        'BD_5pr6_dif':BD_5pr6_dif, 'BD_5pr5_dif':BD_5pr5_dif, 
        'BD_5pr4_dif':BD_5pr4_dif, 'BD_5pr3_dif':BD_5pr3_dif,
        'BD_5pr2_dif':BD_5pr2_dif, 'BD_5pr1_dif':BD_5pr1_dif, 
        'BD_3pr1_dif':BD_3pr1_dif, 'BD_3pr2_dif':BD_3pr2_dif, 
        'BD_3pr3_dif':BD_3pr3_dif}        
u5.violin_plot(BD_dif)
print('Non-isosteric pairs (A--G/G--A, A--A, G--G and C--C)')

print('Histograms of bootstrap difference for U5 bp geometry between -3C and -3Csub introns')
u5.hist_boots_df(BD_5pr8_WC)
print('BD_5pr8_WC')
P_H0_5pr8_WC=u5.boots_P_H0(BD_5pr8_WC)
print('P_H0_5pr8_WC', P_H0_5pr8_WC)
print(' ')
u5.hist_boots_df(BD_5pr8_iso)
print('BD_5pr8_iso')
P_H0_5pr8_iso=u5.boots_P_H0(BD_5pr8_iso)
print('P_H0_5pr8_iso', P_H0_5pr8_iso)
print(' ')
u5.hist_boots_df(BD_5pr8_dif)
print('BD_5pr8_dif')
print(' ')
u5.hist_boots_df(BD_5pr7_WC)
print('BD_5pr7_WC')
P_H0_5pr7_WC=u5.boots_P_H0(BD_5pr7_WC)
print('P_H0_5pr7_WC', P_H0_5pr7_WC)
print(' ')
u5.hist_boots_df(BD_5pr7_iso)
print('BD_5pr7_iso')
P_H0_5pr7_iso=u5.boots_P_H0(BD_5pr7_iso)
print('P_H0_5pr7_iso', P_H0_5pr7_iso)
print(' ')
u5.hist_boots_df(BD_5pr7_dif)
print('BD_5pr7_dif')
P_H0_5pr7_dif=u5.boots_P_H0(BD_5pr7_dif)
print('P_H0_5pr7_dif', P_H0_5pr7_dif)
print(' ')
u5.hist_boots_df(BD_5pr6_WC)
print('BD_5pr6_WC')
P_H0_5pr6_WC=u5.boots_P_H0(BD_5pr6_WC)
print('P_H0_5pr6_WC', P_H0_5pr6_WC)
print(' ')
u5.hist_boots_df(BD_5pr6_iso)
print('BD_5pr6_iso')
P_H0_5pr6_iso=u5.boots_P_H0(BD_5pr6_iso)
print('P_H0_5pr6_iso', P_H0_5pr6_iso)
print(' ')
u5.hist_boots_df(BD_5pr6_dif)
print('BD_5pr6_dif')
P_H0_5pr6_dif=u5.boots_P_H0(BD_5pr6_dif)
print('P_H0_5pr6_dif', P_H0_5pr6_dif)
print(' ')
u5.hist_boots_df(BD_5pr5_WC)
print('BD_5pr5_WC')
P_H0_5pr5_WC=u5.boots_P_H0(BD_5pr5_WC)
print('P_H0_5pr5_WC', P_H0_5pr5_WC)
print(' ')
u5.hist_boots_df(BD_5pr5_iso)
print('BD_5pr5_iso')
P_H0_5pr5_iso=u5.boots_P_H0(BD_5pr5_iso)
print('P_H0_5pr5_iso', P_H0_5pr5_iso)
print(' ')
u5.hist_boots_df(BD_5pr5_dif)
print('BD_5pr5_dif')
print(' ')
u5.hist_boots_df(BD_5pr4_WC)
print('BD_5pr4_WC')
P_H0_5pr4_WC=u5.boots_P_H0(BD_5pr4_WC)
print('P_H0_5pr4_WC', P_H0_5pr4_WC)
print(' ')
u5.hist_boots_df(BD_5pr4_iso)
print('BD_5pr4_iso')
P_H0_5pr4_iso=u5.boots_P_H0(BD_5pr4_iso)
print('P_H0_5pr4_iso', P_H0_5pr4_iso)
print(' ')
u5.hist_boots_df(BD_5pr4_dif)
print('BD_5pr4_dif')
print(' ')
u5.hist_boots_df(BD_5pr3_WC)
print('BD_5pr3_WC')
P_H0_5pr3_WC=u5.boots_P_H0(BD_5pr3_WC)
print('P_H0_5pr3_WC', P_H0_5pr3_WC)
print(' ')
u5.hist_boots_df(BD_5pr3_iso)
print('BD_5pr3_iso')
P_H0_5pr3_iso=u5.boots_P_H0(BD_5pr3_iso)
print('P_H0_5pr3_iso', P_H0_5pr3_iso)
print(' ')
u5.hist_boots_df(BD_5pr3_dif)
print('BD_5pr3_dif')
print(' ')
u5.hist_boots_df(BD_5pr2_WC)
print('BD_5pr2_WC')
P_H0_5pr2_WC=u5.boots_P_H0(BD_5pr2_WC)
print('P_H0_5pr2_WC', P_H0_5pr2_WC)
print(' ')
u5.hist_boots_df(BD_5pr2_iso)
print('BD_5pr2_iso')
P_H0_5pr2_iso=u5.boots_P_H0(BD_5pr2_iso)
print('P_H0_5pr2_iso', P_H0_5pr2_iso)
print(' ')
u5.hist_boots_df(BD_5pr2_dif)
print('BD_5pr2_dif')
print(' ')
u5.hist_boots_df(BD_5pr1_WC)
print('BD_5pr1_WC')
P_H0_5pr1_WC=u5.boots_P_H0(BD_5pr1_WC)
print('P_H0_5pr1_WC', P_H0_5pr1_WC)
print(' ')
u5.hist_boots_df(BD_5pr1_iso)
print('BD_5pr1_iso')
P_H0_5pr1_iso=u5.boots_P_H0(BD_5pr1_iso)
print('P_H0_5pr1_iso', P_H0_5pr1_iso)
print(' ')
u5.hist_boots_df(BD_5pr1_dif)
print('BD_5pr1_dif')
P_H0_5pr1_dif=u5.boots_P_H0(BD_5pr1_dif)
print('P_H0_5pr1_dif', P_H0_5pr1_dif)
print(' ')
u5.hist_boots_df(BD_3pr1_WC)
print('BD_3pr1_WC')
P_H0_3pr1_WC=u5.boots_P_H0(BD_3pr1_WC)
print('P_H0_3pr1_WC', P_H0_3pr1_WC)
print(' ')
u5.hist_boots_df(BD_3pr1_iso)
print('BD_3pr1_iso')
P_H0_3pr1_iso=u5.boots_P_H0(BD_3pr1_iso)
print('P_H0_3pr1_iso', P_H0_3pr1_iso)
print(' ')
u5.hist_boots_df(BD_3pr1_dif)
print('BD_3pr1_dif')
P_H0_3pr1_dif=u5.boots_P_H0(BD_3pr1_dif)
print('P_H0_3pr1_dif', P_H0_3pr1_dif)
print(' ')
u5.hist_boots_df(BD_3pr2_WC)
print('BD_3pr2_WC')
P_H0_3pr2_WC=u5.boots_P_H0(BD_3pr2_WC)
print('P_H0_3pr2_WC', P_H0_3pr2_WC)
print(' ')
u5.hist_boots_df(BD_3pr2_iso)
print('BD_3pr2_iso')
P_H0_3pr2_iso=u5.boots_P_H0(BD_3pr2_iso)
print('P_H0_3pr2_iso', P_H0_3pr2_iso)
print(' ')
u5.hist_boots_df(BD_3pr2_dif)
print('BD_3pr2_dif')
P_H0_3pr2_dif=u5.boots_P_H0(BD_3pr2_dif)
print('P_H0_3pr2_dif', P_H0_3pr2_dif)
print(' ')
u5.hist_boots_df(BD_3pr3_WC)
print('BD_3pr3_WC')
P_H0_3pr3_WC=u5.boots_P_H0(BD_3pr3_WC)
print('P_H0_3pr3_WC', P_H0_3pr3_WC)
print(' ')
u5.hist_boots_df(BD_3pr3_iso)
print('BD_3pr3_iso')
P_H0_3pr3_iso=u5.boots_P_H0(BD_3pr3_iso)
print('P_H0_3pr3_iso', P_H0_3pr3_iso)
print(' ')
u5.hist_boots_df(BD_3pr3_dif)
print('BD_3pr3_dif')
P_H0_3pr3_dif=u5.boots_P_H0(BD_3pr3_dif)
print('P_H0_3pr3_dif', P_H0_3pr3_dif)