#!/usr/bin/env python
# coding: utf-8

# In[23]:


from math import*
from pyopenms import *
seq = AASequence.fromString("VAKA")

print("The peptide", str(seq), "consists of the following amino acids:")
for s in seq:
    print(s.getName(), ":", s.getMonoWeight())
    
listtt = [] 
for r in seq:
    listtt.append(r.getMonoWeight())  

print(listtt)
print(sum(listtt))    
    


# In[16]:


seq.getMonoWeight()


# In[22]:


isEqual=(seq.getMonoWeight()==sum(listtt))
isEqual


# In[ ]:




