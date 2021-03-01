#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 23:09:55 2021

@author: olga
"""

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

print('Hypothesis testing: difference of U6/intron bp geometry frequencies by bootstrapping of -1G and -1Gsub')
print('18 categories: 3 types of bp geometry in 6 intron positions: +5 to +10')
print('18 lists of 10000 bootstrap differencies (BD) of frequencies each for bp geometries at positions of introns preceded by exons with conserved +1G (N=',total_exG,'), and lacking +1G (N=', total_exGsub,')')

BD_in5_WC=[]
BD_in5_iso=[]
BD_in5_dif=[]
BD_in6_WC=[]
BD_in6_iso=[]
BD_in6_dif=[]
BD_in7_WC=[]
BD_in7_iso=[]
BD_in7_dif=[]
BD_in8_WC=[]
BD_in8_iso=[]
BD_in8_dif=[]
BD_in9_WC=[]
BD_in9_iso=[]
BD_in9_dif=[]
BD_in10_WC=[]
BD_in10_iso=[]
BD_in10_dif=[]

i=0
while i<10000:
    
    introns_exG_rand=u5.sample_random(introns_exG, total_exG)
    geo_exG_rand=u5.geo_position_u6_rand(introns_exG_rand)
    geo_freq_exG_rand=u5.geo_freq_position_u6_rand(geo_exG_rand, total_exG)
    if i==0:
        print('geo_freq_exG_rand', geo_freq_exG_rand)
    
    introns_exGsub_rand=u5.sample_random(introns_exGsub, total_exGsub)
    geo_exGsub_rand=u5.geo_position_u6_rand(introns_exGsub_rand)
    geo_freq_exGsub_rand=u5.geo_freq_position_u6_rand(geo_exGsub_rand, total_exGsub)
    if i==0:
        print('geo_freq_exGsub_rand', geo_freq_exGsub_rand)
    #Bootstrap difference in bp geometry frequencies between exG and exGsub
    #Temporariry dictionary is returned by the function and re-written 10000 times with each iteration
    positions=6
    geo_freq_boots_df=u5.geo_freq_boots_df(geo_freq_exG_rand, geo_freq_exGsub_rand, positions)
    if i==0:
        print('geo_freq_boots_df', geo_freq_boots_df)    
    #Bootstrap differences 
    #18 lists to append to: 3 types of bp geometry in 6 intron positions (+5 to +10)
    BD_in5_WC.append(geo_freq_boots_df['Watson_Crick'][0])
    BD_in5_iso.append(geo_freq_boots_df['isosteric'][0])    
    BD_in5_dif.append(geo_freq_boots_df['different'][0])
    BD_in6_WC.append(geo_freq_boots_df['Watson_Crick'][1])
    BD_in6_iso.append(geo_freq_boots_df['isosteric'][1])    
    BD_in6_dif.append(geo_freq_boots_df['different'][1])
    BD_in7_WC.append(geo_freq_boots_df['Watson_Crick'][2])
    BD_in7_iso.append(geo_freq_boots_df['isosteric'][2])    
    BD_in7_dif.append(geo_freq_boots_df['different'][2]) 
    BD_in8_WC.append(geo_freq_boots_df['Watson_Crick'][3])
    BD_in8_iso.append(geo_freq_boots_df['isosteric'][3])    
    BD_in8_dif.append(geo_freq_boots_df['different'][3])
    BD_in9_WC.append(geo_freq_boots_df['Watson_Crick'][4])
    BD_in9_iso.append(geo_freq_boots_df['isosteric'][4])    
    BD_in9_dif.append(geo_freq_boots_df['different'][4])
    BD_in10_WC.append(geo_freq_boots_df['Watson_Crick'][5])
    BD_in10_iso.append(geo_freq_boots_df['isosteric'][5])    
    BD_in10_dif.append(geo_freq_boots_df['different'][5])   
        
    i+=1
    if i==10000:
        break

print('Violin plots: bootstrap difference of U6 bp geometry following -1G and -1Gsub exons')
#Dictionaries of BD lists to create dfs for violin plots
BD_WC={'BD_in5_WC':BD_in5_WC, 'BD_in6_WC':BD_in6_WC, 
       'BD_in7_WC':BD_in7_WC, 'BD_in8_WC':BD_in8_WC, 
       'BD_in9_WC':BD_in9_WC, 'BD_in10_WC':BD_in10_WC}
u5.violin_plot(BD_WC)
print('Watson-Crick pairs')

BD_iso={'BD_in5_iso':BD_in5_iso, 'BD_in6_iso':BD_in6_iso, 
        'BD_in7_iso':BD_in7_iso, 'BD_in8_iso':BD_in8_iso, 
        'BD_in9_iso':BD_in9_iso, 'BD_in10_iso':BD_in10_iso}
u5.violin_plot(BD_iso)
print('Isosteric pairs (G--U/U--G, U--U, A--C/C--A and C--U/U--C)')

BD_dif={'BD_in5_dif':BD_in5_dif, 'BD_in6_dif':BD_in6_dif, 
        'BD_in7_dif':BD_in7_dif, 'BD_in8_dif':BD_in8_dif, 
        'BD_in9_dif':BD_in9_dif, 'BD_in10_dif':BD_in10_dif}        
u5.violin_plot(BD_dif)
print('Non-isosteric pairs (A--G/G--A, A--A, G--G and C--C)')

print('Histograms of bootstrap difference for U6 bp geometry in introns preceded by exon-end Gs and lacking exon-end Gs')

u5.hist_boots_df(BD_in5_WC)
print('BD_in5_WC')
P_H0_in5_WC=u5.boots_P_H0(BD_in5_WC)
print('P_H0_in5_WC', P_H0_in5_WC)
print(' ')

u5.hist_boots_df(BD_in5_iso)
print('BD_in5_iso')
P_H0_in5_iso=u5.boots_P_H0(BD_in5_iso)
print('P_H0_in5_iso', P_H0_in5_iso)
print(' ')

u5.hist_boots_df(BD_in5_dif)
print('BD_in5_dif')
P_H0_in5_dif=u5.boots_P_H0(BD_in5_dif)
print('P_H0_in5_dif', P_H0_in5_dif)
print(' ')

u5.hist_boots_df(BD_in6_WC)
print('BD_in6_WC')
P_H0_in6_WC=u5.boots_P_H0(BD_in6_WC)
print('P_H0_in6_WC', P_H0_in6_WC)
print(' ')

u5.hist_boots_df(BD_in6_iso)
print('BD_in6_iso')
P_H0_in6_iso=u5.boots_P_H0(BD_in6_iso)
print('P_H0_in6_iso', P_H0_in6_iso)
print(' ')

u5.hist_boots_df(BD_in6_dif)
print('BD_in6_dif')
P_H0_in6_dif=u5.boots_P_H0(BD_in6_dif)
print('P_H0_in6_dif', P_H0_in6_dif)
print(' ')

u5.hist_boots_df(BD_in7_WC)
print('BD_in7_WC')
P_H0_in7_WC=u5.boots_P_H0(BD_in7_WC)
print('P_H0_in7_WC', P_H0_in7_WC)
print(' ')

u5.hist_boots_df(BD_in7_iso)
print('BD_in7_iso')
P_H0_in7_iso=u5.boots_P_H0(BD_in7_iso)
print('P_H0_in7_iso', P_H0_in7_iso)
print(' ')

u5.hist_boots_df(BD_in7_dif)
print('BD_in7_dif')
print(' ')
#Uracil pair: non-isosteric (dif) pairs do not exist

u5.hist_boots_df(BD_in8_WC)
print('BD_in8_WC')
P_H0_in8_WC=u5.boots_P_H0(BD_in8_WC)
print('P_H0_in8_WC', P_H0_in8_WC)
print(' ')

u5.hist_boots_df(BD_in8_iso)
print('BD_in8_iso')
P_H0_in8_iso=u5.boots_P_H0(BD_in8_iso)
print('P_H0_in8_iso', P_H0_in8_iso)
print(' ')

u5.hist_boots_df(BD_in8_dif)
print('BD_in8_dif')
P_H0_in8_dif=u5.boots_P_H0(BD_in8_dif)
print('P_H0_in8_dif', P_H0_in8_dif)
print(' ')

u5.hist_boots_df(BD_in9_WC)
print('BD_in9_WC')
P_H0_in9_WC=u5.boots_P_H0(BD_in9_WC)
print('P_H0_in9_WC', P_H0_in9_WC)
print(' ')

u5.hist_boots_df(BD_in9_iso)
print('BD_in9_iso')
P_H0_in9_iso=u5.boots_P_H0(BD_in9_iso)
print('P_H0_in9_iso', P_H0_in9_iso)
print(' ')

u5.hist_boots_df(BD_in9_dif)
print('BD_in9_dif')
P_H0_in9_dif=u5.boots_P_H0(BD_in9_dif)
print('P_H0_in9_dif', P_H0_in9_dif)
print(' ')

u5.hist_boots_df(BD_in10_WC)
print('BD_in10_WC')
P_H0_in10_WC=u5.boots_P_H0(BD_in10_WC)
print('P_H0_in10_WC', P_H0_in10_WC)
print(' ')

u5.hist_boots_df(BD_in10_iso)
print('BD_in10_iso')
P_H0_in10_iso=u5.boots_P_H0(BD_in10_iso)
print('P_H0_in10_iso', P_H0_in10_iso)
print(' ')

u5.hist_boots_df(BD_in10_dif)
print('BD_in10_dif')
P_H0_in10_dif=u5.boots_P_H0(BD_in10_dif)
print('P_H0_in10_dif', P_H0_in10_dif)