#This takes protein information from a FASTA file and then uses protparam to analyse protein physical properities
#It will output a lsit of protein IDs, and lists of physical properties
#These can be joined to make a dataframe.

#Importing the necesarry packages
import Bio

from Bio.SeqUtils import ProtParam

import itertools
import pandas as pd
import numpy as np
import regex
from collections import OrderedDict


#set infile name
infile = "aa_freqs.txt"

#create empty dataframe to update with dipeptide frequencies
cns=["peptide","pept_length"]
aa_list = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
protparam=cns+["aromatic","gravy","insta","isoelectric","flexibility"]+aa_list
    
dipept_list2 = []
dipept_list1 = itertools.product("ACDEFGHIKLMNPQRSTVWY",repeat=2) #get all dipeptides

#cleanup dipeptide list
for dipept in dipept_list1:
    dipept = str(dipept)
    dipept = dipept.replace('\'','')
    dipept = dipept.replace(')','')
    dipept = dipept.replace(', ','')
    dipept = dipept.replace('(','')
    dipept_list2.append(dipept)
cns = cns + dipept_list2

#opening an output file along with the fasta file containing the protein ID and sequence
f = open(infile, 'r')

#set up empty dataframes with appropriate length from infile    
#df_dipeptides = pd.DataFrame(data=pd.read_csv(infile),columns=cns)
#df_dipeptides[["peptide"]] = df_dipeptides[["peptide"]].to_string

df_all = pd.DataFrame(data=pd.read_csv(infile),columns=protparam)
df_all[["peptide"]] = df_all[["peptide"]].to_string
print("Starting")

#set counters to keep track of df rows
counter= 0
count_hun = 0
count_mil = 0

#if length is over 1 mil, split up dataframe
if len(df_all)>1000000:
    df_protparam = pd.DataFrame(data=df_all.iloc[0:999999])
else:
    df_protparam = df_all

#getting all possible dipeptides
for line in f.readlines(): # reading each line in the file
    seq=line.strip() #removing the hidden new line characters in each line
    dipeptides = regex.findall("..",seq,overlapped=True) #get all dipeptides from a given sequence
    unique_dipeptides = list(set(dipeptides)) # get all unique dipeptides for a given sequence
    pept_length = len(seq) # get peptide length
    
    #update df to include peptide and peptide length
    #df_dipeptides.set_value(counter, "peptide", str(seq))
    #df_dipeptides.set_value(counter, "pept_length", pept_length)
    
    #for dipept in unique_dipeptides: #loop through unique dipeptides
        #df_dipeptides.set_value(counter, dipept, dipeptides.count(dipept)) #update dipeptide counts per peptide
    
    
    aa_seq = regex.findall(".",seq) #get all aa from a given sequence
    unique_aa = list(set(aa_seq)) # get all unique aa for a given sequence
    
    for aas in unique_aa: #loop through unique aa
        df_protparam.set_value(counter, aas, aa_seq.count(aas)) #update aa counts per peptide

    #Add properties to df
    X = ProtParam.ProteinAnalysis(seq) #analysing with protparam
    
    #update df to include peptide and peptide length
    df_protparam.set_value(counter, "peptide", str(seq))
    df_protparam.set_value(counter, "pept_length", pept_length)
    
    df_protparam.set_value(counter, "aromatic", X.aromaticity()) #adding aromaticity to df
    df_protparam.set_value(counter, "gravy", X.gravy()) #adding gravy score to df
    df_protparam.set_value(counter, "insta", X.instability_index()) #adding instability to df
    df_protparam.set_value(counter, "isoelectric", X.isoelectric_point()) #adding isoelectric point to df
    df_protparam.set_value(counter, "flexibility", np.mean(X.flexibility())) #adding flexibility score to df. I took the mean here as it returns several values
    
    counter +=1
    if counter % 100000 == 0: #print out a number for each 100,000 peptides done
        count_hun += 1
        print(count_hun)
    
    #added this to avoid an error when trying to update with an index >1 mil.
    #I split the dataset into 1 mil samples at a time, get the values, then append to the end of
    #the next 1 mil set    
    if counter % 1000000 == 0 :
        count_mil += 1 # count how many millions
        start = count_mil * 1000000 #get start of next set
        stop = min(len(df_all),((count_mil+1)*1000000)-1) #get end of next set
        counter = 0 # restart counter
        df_protparam = pd.concat([df_all.iloc[start:stop],df_protparam],ignore_index=True) #merge new dataset to before old dataset

f.close()	

#convert NAs to 0s
#df_dipeptides = df_dipeptides.fillna(value=0)
df_protparam.fillna(value=0)

#combine dataframes to have one with all values
#df3 = pd.merge(df_protparam,df_dipeptides,left_on=["peptide","pept_length"],right_on=["peptide","pept_length"])


#write dataframes to tsv files
df_protparam.to_csv('aa_freq_props.tsv',sep='\t', index=False)

#df_dipeptides.to_csv('random_peptide_dipeptides.tsv',sep='\t', index=False)
#df3.to_csv('random_peptide_all_properties.tsv',sep='\t', index=False)


