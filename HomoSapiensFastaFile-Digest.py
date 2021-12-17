#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pyopenms import *

 
dig = ProteaseDigestion()
dig.getEnzymeName() # Trypsin

entries = []
f = FASTAFile()
f.load("uniprot-reviewed_yes+taxonomy_9606 (1).fasta", entries)
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
    print ("SEQUENCE") 
    print(u)
    print (10*"=",'\n',"After Digestion SEQUENCE",'\n',10*"=") 
    resultOneSeq.append( dig.digest(u, result))
    for re in result:
        print(re)
    print (100*"=")
    resultAll.append(resultOneSeq)


print (100*"*")

   


# In[ ]:





# In[ ]:




