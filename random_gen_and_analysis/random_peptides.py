#!/usr/bin/python 
import random

aromatic_max = "HFYWP"
aliphatic_max = "AVIL"
gravy_min = "RKNDEQ"
gravy_max = "IVLCF"
insta_min = "GTVNF"
insta_max = "SPMCE"
iso_max = "CKRY"
iso_min = "DEH"
flexi_max = "EDNSP"
flexi_min = "CWFIY"
all_aa = 'ACDEFGHIKLMNPQRSTVWY'
double_aa = "RRRFFFWWWIIIYYYEEECDHKLNPQSV"
all_aa = all_aa+double_aa
#all_aa= "CRILVAMFSEDKHQN"
list_in = [aromatic_max,aliphatic_max,gravy_min,gravy_max,insta_min,insta_max,iso_max,iso_min,flexi_max,flexi_min,all_aa]


#function that will generate a random protein sequence
def protein(length,aa_list):
	return ''.join(random.choice(aa_list) for _ in xrange(length))


number_of_sequences=raw_input("How many random peptide sequences do you want? ") #asks for the number of random sequences you want to generate
peptide_length=raw_input("How long do you want the peptide sequences to be? ")# asks for the length of the random peptide
extremes=raw_input("Do you want to include proteins with extreme phenotypes? ")
output_name=raw_input("What file do you want to write to? ")
peptide_length=int(peptide_length) #making sure this variable is an integer
peptide_length=peptide_length-1
number_of_sequences=int(number_of_sequences) #making sure this variable is an integer
counter=0 #setting counter to 0

f=open(output_name,'w') #output file


while counter<number_of_sequences: #while counter is less than the number of sequences that you want
    
    if extremes == "yes":
        if counter / 10000 < 10:
            list_count = counter / 10000
        else:
            list_count = 10
    else:
        list_count = 10
        
    if counter % 10000 == 0:
        print(counter)
    
    protein_seq=protein(peptide_length,list_in[list_count]) #get a random protein sequence
    protein_seq='M'+protein_seq #adding M to the protein as translated peptides start with M
    f.write(protein_seq+'\n') #write this random sequence to a file
    counter+=1 #add one to the counter

f.close() #closing up the file
