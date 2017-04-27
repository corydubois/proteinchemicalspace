#!/usr/bin/python 
#this program will create a 100K peptide library with the extreme values for each peptide library
#from an input set of 5 million or 10 million peptides. The input file has the peptide library
#and the chemical properties of each peptide

#importing the necessary packages
import pandas as pd
import numpy as np

#getting the input file (should be a large peptide library with their chemical properties)
inputfile=raw_input("What is the file name? It should be in tab separated format. ")
#input of the file
df1 = pd.read_csv(inputfile, sep='\t', header=0)
#creating a list of all columns in the input dataframe
properties=list(df1.columns.values)
#sorting the dataframe by the aromaticity score
df1 = df1.sort_values(by='aromatic')
#taking the 1500 peptides with the highest aromaticity score
df5=df1.head(1500)
#taking the 1500 peptides with the lowest aromaticity score
df4=df1.tail(1500)
#adding the highest and lowest peptides to a new dataframe
df5=df5.append(df4,ignore_index=True)

#remove the aromaticity and peptide from the list of pwptide properties
properties.remove('peptide')
properties.remove('aromatic')

#take top 1500 and bottom 1500 for the other chemical properties, and add these to a new dataframe
for columnname in properties:
	df1=df1.sort_values(by=columnname)
	df3=df1.head(1500)
	df4=df1.tail(1500)
	df5=df5.append(df3, ignore_index=True)
	df5=df5.append(df4, ignore_index=True)

#remove duplicates from the peptides in the datafrane containing the top and bottom peptides
df5=df5.drop_duplicates()

#add random peptides from the 5 million peptide library to the table of top/bottom peptides until we have 100,000 peptides
number_rows=100000-len(df5.index)
df6=df1.sample(n=number_rows, replace=False)
frames = [df5,df6]
df7=pd.concat(frames)
#export just the top and bottom peptides to a tsv file
df5.to_csv('top_peptides5mil_mostvar.tsv', sep='\t',index=False)
#export the 100K library to a tsv file
df7.to_csv('select100K_5mil_mostvar.tsv',sep='\t', index=False)



