#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 10:20:16 2020

@author: olga
"""

import U5 as u5
import human_genes_saved as g


print('Extracting slice junctions from 132 human genes')
files=g.exon_files()
splice_junctions=[]
for file in files:
    print(file)
    s_junctions=u5.e_splice_junctions(file)
    print('Splice junctions in this gene: ',len(s_junctions))
    splice_junctions+=s_junctions
total=len(splice_junctions)
print('Splice junctions total: ',total)

print('Extracting introns (starts and ends) from 132 human genes')

files=g.intron_files()
intron_starts=[]
intron_ends=[]
for file in files:
    gene_introns=[]
    print(file)
    gene_introns=u5.introns(file)
    print('Introns in this gene: ',len(gene_introns))
    gene_intron_starts=u5.intron_starts(gene_introns)
    intron_starts+=gene_intron_starts
    gene_intron_ends=u5.intron_ends(gene_introns)
    intron_ends+=gene_intron_ends
total1=len(intron_starts)
total2=len(intron_ends)
total=0
if total1==total2:
    total=total1
print('Introns total: ',total)
#print(intron_starts)
#print(intron_ends)

print('Unusual introns')

print('Isolating all minor (U12) introns')
print('Looking for U6atac binding site motif ATT(+5C)C')
i=0
u12=[]#List of indeces for U12 introns
for intron in intron_starts:
    if 'ATCC' in intron and intron[4]=='C':
        u12.append(i)
        print(i+1, splice_junctions[i], intron, intron_ends[i])
    i+=1
print(u12)
percent_u12=(len(u12)/total)*100
print('Total U12 introns: ',len(u12), '{:6.2f}'.format(percent_u12),'%')
print('Looking for U12 binding site')
i=0
for intron in intron_ends:
    if 'CCTTAAC' in intron:
        print(i+1, splice_junctions[i], intron_starts[i], intron)
    i+=1
print('Caution: U12 motif CCTTAAC and U2 motif TAAC are not easy to separate')
print('Isolating GC-AG introns')
i=0
GC_AG=[]#List of indeces for GC(A)_AG introns
for intron in intron_starts:
    if intron[1]!='T':
        GC_AG.append(i)
        print(i+1, splice_junctions[i], intron, intron_ends[i])
    i+=1
print(GC_AG)
percent_GC_AG=(len(GC_AG)/total)*100
print('Total GC(A)_AG introns: ',len(GC_AG), '{:6.2f}'.format(percent_GC_AG),'%')
print('All GC(A)-AC introns have multiple W-C pairs in U5 and U6 helices at the ex/in boundary.')   

print('Isolating AU-AC introns')
i=0
print('Testing for the first intron G')
for intron in intron_starts:
    if intron[0]!='G':
        print(i+1, splice_junctions[i], intron, intron_ends[i])
    i+=1
print('Testing for the last intron G')
found=0
for intron in intron_ends:
    if intron[-1]!='G':
        print(i+1, splice_junctions[i], intron_starts[i], intron)
        found=1
    i+=1
if found==0:
    print('All the examined introns end with a G')

print('Excluding the minor spliceosome (U12) introns from the analysis')
i=0
splice_junctions_u2=[]#Splice junctions of the major spliceosome only
for junction in splice_junctions:
    if i not in u12:
        splice_junctions_u2.append(splice_junctions[i])
    i+=1
total_u2=len(splice_junctions_u2)    
print('Total splice junctions of major introns only: ',total_u2)
i=0
intron_starts_u2=[]#Start sequences (10nt) of the major introns only 
for intron in intron_starts:
    if i not in u12:
        intron_starts_u2.append(intron_starts[i])
    i+=1
total_u2_introns=len(intron_starts_u2)    
if total_u2!=total_u2_introns:
    print(total_u2_introns)
i=0
intron_ends_u2=[]#End sequences (60nt) of the major introns only 
for intron in intron_ends:
    if i not in u12:
        intron_ends_u2.append(intron_ends[i])
    i+=1
total_u2_intends=len(intron_ends_u2)    
if total_u2!=total_u2_intends:
    print(total_u2_intends)

print('Data for all major (U2) splicesome introns, N=2003')
print('U5 Loop1 base pairs listed by position')
bp_position=u5.bp_position(splice_junctions_u2)

print('Frequency of U5 Loop1 base pairs listed by position')
bp_freq_position=u5.bp_freq_position(bp_position, total_u2)

print('U5 Loop1 base pair geometry by position')
geo_position=u5.geo_position(splice_junctions_u2)

print('U5 base pair geometry: frequencies by position')
geo_freq_position=u5.geo_freq_position(geo_position, total_u2)
u5.geo_stacked_barchart(geo_freq_position)

print('Excluding GC(A)-AG introns for the analysis of the U5-U6 compensation')
i=0
splice_junctions_u6=[]#Splice junctions of the major GU-AG introns only
for junction in splice_junctions_u2:
    if i not in GC_AG:
        splice_junctions_u6.append(splice_junctions_u2[i])
    i+=1
total_u6=len(splice_junctions_u6)    
print('Total splice junctions of the major GU-AG introns only: ',total_u6)
i=0
intron_starts_u6=[]#Start sequences (10nt) of the major GU-AG introns only 
for intron in intron_starts_u2:
    if i not in GC_AG:
        intron_starts_u6.append(intron_starts_u2[i])
    i+=1
total_u6_introns=len(intron_starts_u6)    
if total_u6!=total_u6_introns:
    print(total_u6_introns)
i=0
intron_ends_u6=[]#End sequences (60nt) of the major GU-AG introns only 
for intron in intron_ends_u2:
    if i not in GC_AG:
        intron_ends_u6.append(intron_ends_u2[i])
    i+=1
total_u6_intends=len(intron_ends_u6)    
if total_u6!=total_u6_intends:
    print(total_u6_intends)