#!/usr/bin/env python
# coding: utf-8

# In[9]:


from pyopenms import *
dig = ProteaseDigestion()
dig.getEnzymeName()

entries = []
f = FASTAFile()
f.load("uniprot-id-P12688+OR+id-P16521+OR+id-P25491+OR+id-Q04437+OR+id-P41800+--.fasta", entries)
print( len(entries) )

bsa2=[]
for e in entries:
    bsa2.append(AASequence.fromString( e.sequence))
   # print ( e.sequence)    
   # print (100*"=")
    
resultAll = []
resultOneSeq=[]
result=[]
for u in bsa2:
    print ("SEQUENCE",'\n') 
    print(u)
    print (10*"=",'\n','\n',"After Digestion SEQUENCE",'\n',10*"=",'\n') 
    resultOneSeq.append( dig.digest(u, result))
    for re in result:
        print(re)
    print (120*"_")
    resultAll.append(resultOneSeq)


print (100*"*")


# In[ ]:




