#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 10:23:25 2020

@author: olga
"""

#splice_sites.py should be run first

import U5 as u5

  
print('End of exon and start of intron interactions with U5 and U6 snRNAs')
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
geo_freq_inGsub=u5.geo_freq_position(geo_inGsub, total_inG_sub)

print('N=',total_inG,',+5G present: U5 base pair geometry')
print('U5 Loop1 bp geometry by position for all splice junctions with +5G introns')
geo_inG=u5.geo_position(sjunctions_inG)

print('U5 Loop1 bp geometry: frequency by position for all splice junctions with +5G introns')
geo_freq_inG=u5.geo_freq_position(geo_inG, total_inG)

print('Stacked barcharts for U5 bp geometry frequencies by position')
u5.geo_stacked_barchart(geo_freq_inGsub)
print('U5 Loop1 bp geometry: frequency by position for splice junctions with introns missing +5G')
print('N=',total_inG_sub)
print(' ')
u5.geo_stacked_barchart(geo_freq_inG)
print('U5 Loop1 bp geometry: frequency by position for all splice junctions with +5G introns')
print('N=',total_inG)
print(' ')

print('Hypothesis testing: difference of U5/exons bp geometry frequencies by bootstrapping of +5G and +5Gsub')
print('33 categories: 3 types of bp geometry in 11 splice junction positions')
print('33 lists of 10000 bootstrap differencies (BD) of frequencies each for bp geometries at positions of splice junctions for introns with conserved +5G (N=',total_inG,'), and lacking +5G (N=', total_inG_sub,')')

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
    
    sjunctions_inG_rand=u5.sample_random(sjunctions_inG, total_inG)
    geo_inG_rand=u5.geo_position_rand(sjunctions_inG_rand)
    geo_freq_inG_rand=u5.geo_freq_position_rand(geo_inG_rand, total_inG)
    if i==0:
        print('geo_freq_inG_rand', geo_freq_inG_rand)
    
    sjunctions_inGsub_rand=u5.sample_random(sjunctions_inG_sub, total_inG_sub)
    geo_inGsub_rand=u5.geo_position_rand(sjunctions_inGsub_rand)
    geo_freq_inGsub_rand=u5.geo_freq_position_rand(geo_inGsub_rand, total_inG_sub)
    if i==0:
        print('geo_freq_inGsub_rand', geo_freq_inGsub_rand)
    #Bootstrap difference in bp geometry frequencies between inG and inGsub
    #Temporariry dictionary is returned by the function and re-written 10000 times with each iteration
    positions=11
    geo_freq_boots_df=u5.geo_freq_boots_df(geo_freq_inG_rand, geo_freq_inGsub_rand, positions)
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
print('Violin plots: bootstrap difference of U5 bp geometry between +5G and +5Gsub introns')
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

print('Histograms of bootstrap difference for U5 bp geometry between +5G and +5Gsub introns')
u5.hist_boots_df(BD_5pr8_WC)
P_H0_5pr8_WC=u5.boots_P_H0(BD_5pr8_WC)
print('BD_5pr8_WC')
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

print('U1 snRNA interaction with position 5pr3 in the presence and absence of intron +5G')

print('Substitutions of +5G (inGsub), total introns: ',total_inG_sub)
print('Conserved +5G (inG), total introns: ',total_inG)

print('Nucleotide frequencies at position -3 of the 5prime exon (5pr3) followed by introns missing +5G (inGsub)')
print('inGsub 5pr3 A:', bp_freq_inGsub['AU'][5])
print('inGsub 5pr3 C:', bp_freq_inGsub['CU'][5])
print('inGsub 5pr3 G:', bp_freq_inGsub['GU'][5]) 
print('inGsub 5pr3 U:', bp_freq_inGsub['UU'][5]) 

print('Nucleotide frequencies at position -3 of the 5prime exon (5pr3) followed by introns with conserved +5G (inG)')
print('inG 5pr3 A:', bp_freq_inG['AU'][5])
print('inG 5pr3 C:', bp_freq_inG['CU'][5])
print('inG 5pr3 G:', bp_freq_inG['GU'][5]) 
print('inG 5pr3 U:', bp_freq_inG['UU'][5]) 

u5.bp_stacked_barchart_5pr3(bp_freq_inGsub, bp_freq_inG)
print('Nucleotide frequencies at exon position 5pr3 for inG, N=', total_inG, 'and inGsub, N=', total_inG_sub) 
print(' ')

print('U1 vs U5: increased occurrence of Watson-Crick pairs to compensate for substitutions of +5G')
print('Hypothesis testing: difference of nucleotide frequencies at exon position 5pr3 by bootstrapping of +5G and +5Gsub')
print('4 categories: 4 nucleotides in 1 exon position')
print('4 lists of 10000 bootstrap differencies (BD) of frequencies for each nucleotide at exon position 5pr3 for introns with conserved +5G (N=',total_inG,'), and lacking +5G (N=', total_inG_sub,')')

BD_5pr3_A=[]
BD_5pr3_C=[]
BD_5pr3_G=[]
BD_5pr3_U=[]

i=0
while i<10000:
    
    sjunctions_inG_rand=u5.sample_random(sjunctions_inG, total_inG)
    Nfreq_5pr3_inG_rand=u5.N_freq_5pr3_rand(sjunctions_inG_rand, total_inG)
    if i==0:
        print('Nfreq_5pr3_inG_rand', Nfreq_5pr3_inG_rand)
    elif i==1:
        print('Nfreq_5pr3_inG_rand1', Nfreq_5pr3_inG_rand)
    elif i==2:
        print('Nfreq_5pr3_inG_rand2', Nfreq_5pr3_inG_rand)        
    sjunctions_inGsub_rand=u5.sample_random(sjunctions_inG_sub, total_inG_sub)
    Nfreq_5pr3_inGsub_rand=u5.N_freq_5pr3_rand(sjunctions_inGsub_rand, total_inG_sub)
    if i==0:
        print('Nfreq_5pr3_inGsub_rand', Nfreq_5pr3_inGsub_rand)
    elif i==1:
        print('Nfreq_5pr3_inGsub_rand1', Nfreq_5pr3_inGsub_rand)
    elif i==2:
        print('Nfreq_5pr3_inGsub_rand2', Nfreq_5pr3_inGsub_rand)
    #Bootstrap difference in bp geometry frequencies between inG and inGsub
    #Temporariry dictionary is returned by the function and re-written 10000 times with each iteration
    Nfreq_5pr3_boots_df=u5.Nfreq_5pr3_boots_df(Nfreq_5pr3_inG_rand, Nfreq_5pr3_inGsub_rand)
    if i==0:
        print('Nfreq_5pr3_boots_df', Nfreq_5pr3_boots_df)
    elif i==1:
        print('Nfreq_5pr3_boots_df1', Nfreq_5pr3_boots_df)
    elif i==2:
        print('Nfreq_5pr3_boots_df2', Nfreq_5pr3_boots_df)
    #Bootstrap differences 
    #4 lists to append to: 4 nucleotides in exon position 5pr3
    BD_5pr3_A.append(Nfreq_5pr3_boots_df['A'])
    BD_5pr3_C.append(Nfreq_5pr3_boots_df['C'])   
    BD_5pr3_G.append(Nfreq_5pr3_boots_df['G'])
    BD_5pr3_U.append(Nfreq_5pr3_boots_df['U'])
        
    i+=1
    if i==10000:
        break


print('Histograms of bootstrap difference for Nfreq at 5pr3 between +5G and +5Gsub introns')
u5.hist_boots_df(BD_5pr3_A)
#CI99_BD_5pr8_WC=u5.boots_confidence_interval_99(BD_5pr8_WC)
#CI999_BD_5pr8_WC=u5.boots_confidence_interval_999(BD_5pr8_WC)
#print('BD_5pr8_WC, CI 99%:',CI99_BD_5pr8_WC,'CI 99.9%:',CI999_BD_5pr8_WC)
print('BD_5pr3_A')
P_H0_5pr3_A=u5.boots_P_H0(BD_5pr3_A)
print('P_H0_5pr3_A', P_H0_5pr3_A)
print(' ')
u5.hist_boots_df(BD_5pr3_C)
print('BD_5pr3_C')
P_H0_5pr3_C=u5.boots_P_H0(BD_5pr3_C)
print('P_H0_5pr3_C', P_H0_5pr3_C)
print(' ')
u5.hist_boots_df(BD_5pr3_G)
print('BD_5pr3_G')
P_H0_5pr3_G=u5.boots_P_H0(BD_5pr3_G)
print('P_H0_5pr3_G', P_H0_5pr3_G)
print(' ')
u5.hist_boots_df(BD_5pr3_U)
print('BD_5pr3_U')
P_H0_5pr3_U=u5.boots_P_H0(BD_5pr3_U)
print('P_H0_5pr3_U', P_H0_5pr3_U)
print(' ')
