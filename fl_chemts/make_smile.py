
# coding: utf-8

# In[29]:


import csv
import itertools
import operator
import numpy as np
import nltk
import os
from rdkit import Chem
from rdkit.Chem import Draw
from IPython import display
#import matplotlib.pyplot as plt
from rdkit.Chem import Descriptors



def zinc_processed_with_bracket(sen_space):
    #print sen_space
    all_smile=[]
    length=[]
    end="\n"
    element_table=["C","N","B","O","P","S","F","Cl","Br","I","(",")","=","#","Si"]
    ring=["1","2","3","4","5","6","7","8","9","10"]

    for i in range(len(sen_space)):
        #word1=sen_space[i]
        word_space=sen_space[i]
        word=[]
        #word_space.insert(0,end)
        j=0
        while j<len(word_space):
            word_space1=[]
            #word_space1.append(word_space[j])
            if word_space[j]=="[":
                word_space1.append(word_space[j])
                j=j+1
                while word_space[j]!="]":
                    word_space1.append(word_space[j])
                    j=j+1
                word_space1.append(word_space[j])
                word_space2=''.join(word_space1)
                
                    
                word.append(word_space2)
                j=j+1
            else:
                word_space1.append(word_space[j])

                if j+1<len(word_space):
                    word_space1.append(word_space[j+1])
                    word_space2=''.join(word_space1)
                else:
                    word_space1.insert(0,word_space[j-1])
                    word_space2=''.join(word_space1)

                if word_space2 not in element_table:
                    word.append(word_space[j])
                    j=j+1
                else:
                    word.append(word_space2)
                    j=j+2
                    
        
        word.append(end)
        word.insert(0,"&")
        len1=len(word)
        
        length.append(len1)
        if '[SiH2]' not in list(word):
            if '[SiH3]' not in list(word):
                if '[SiH]' not in list(word):
                    all_smile.append(list(word))
       
            
        #si=['[SiH2]','[SiH3]','[SiH]']
    #print all_smile[4]
            
            
        
    after_all_smile=all_smile
    print len(after_all_smile)
            
    
    
    val=["\n"]
    delid=[]
    all_smile_go=[]
    for i in range(len(after_all_smile)):
        for j in range(len(after_all_smile[i])):
            if after_all_smile[i][j] not in val:
                val.append(after_all_smile[i][j])
    #print delid
    #print val
    #val.remove("\n")
    #val.insert(0,"\n")
    #print val
    #print all_smile[0]
    #print all_smile[1]
    #print all_smile[2]
    #print len(all_smile)
    print max(length)
    #print len(val)

    return val, after_all_smile





def zinc_logp(smile):
    logp_value=[]
    compound=[]

    for i in range(len(smile)):
        m = Chem.MolFromSmiles(smile[i])
        if m!=None:
            compound.append(smile[i])
            logp=Descriptors.MolLogP(m)
            logp_value.append(logp)
            

    #print max(logp_value)
    #print logp_value
    
    return compound


def zinc_data_with_bracket_original():

    sen_space=[]
    #f = open('/Users/yang/smiles.csv', 'rb')
    f = open('/home/yang/DP-ChemTS/leaf_parallel_wv_less/124k.csv', 'rb')
    #f = open('/Users/yang/ChemTS/data/250k_rndm_zinc_drugs_clean.smi', 'rb')

    reader = csv.reader(f)
    for row in reader:
        #word_space[row].append(reader[row])
        #print word_sapce
        sen_space.append(row)
    #print sen_space
    f.close()
    
    print len(sen_space)
    #print sen_space[0]
    #print sen_space[1][0]
    #sen_space=sen_space[10:]
    #print sen_space[0][0]
    #print sen_space[0][1]
    #print sen_space[0]
    
    

    word1=sen_space
    #word_space=list(word1[0])
    end="\n"

    zinc_processed=[]
    organic_smile=[]
    t=0
    for i in range(len(sen_space)):
        word1=sen_space[i]
        #print word1[0]
        #m = Chem.MolFromSmiles(word1[0])
        #Chem.Kekulize(m)
        #s=Chem.MolToSmiles(m,kekuleSmiles=True)
        if word1!=[]:
            zinc_processed.append(word1[0])
        #word_space=list(word1[0])
    #print len(zinc_processed)

    #while t <len(zinc_processed):
    #    #print t
    #    word2=zinc_processed[t]
    #    word_space=list(word2)
    #    word=[]

    #    organic_smile.append(word_space)
    #    t=t+1

    #print len(organic_smile)
    #print organic_smile
    #print zinc_processed[0]
    print zinc_processed[1]
    return zinc_processed


#smile_old=zinc_data_with_bracket_original()
#com=zinc_logp(smile_old)
#print len(com)
#val,smile=zinc_processed_with_bracket(com)
#print val




#hi=organic()
#organic_logp(hi)

#hi=zinc_data_with_bracket_original()
#good=zinc_logp(hi)
#zinc_processed_with_bracket(hi)




