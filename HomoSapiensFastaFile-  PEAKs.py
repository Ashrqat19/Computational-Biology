#!/usr/bin/env python
# coding: utf-8

# In[3]:


from pyopenms import *
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = "".join([l.strip() for l in open("uniprot-reviewed_yes+taxonomy_9606 (1).fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
# create all digestion products
result = []
dig.digest(bsa, result)
peptides = [AASequence.fromString(s.toString()) for s in result]

 for peptide in peptides:
    tsg = TheoreticalSpectrumGenerator()
    spec = MSSpectrum()

     p = Param()
    p.setValue("add_b_ions", "false")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec, peptide, 1, 1) # charge range 1:1
    print("Spectrum 1 of", peptide, "has", spec.size(), "peaks.")
    for ion, peak in zip(spec.getStringDataArrays()[0], spec):
        print(ion.decode(), "is generated at m/z", peak.getMZ())


# In[ ]:




