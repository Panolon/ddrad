#!/usr/bin/env python
# coding: utf-8

# ## DDrad sequencing
# 
# In this project we try to implement DDrad sequencing. In order to eliminate the cost of sequencing we utilize restriction enzymes. Those enzymes cut the dna sequence in specific fragments related to their site and we end up with subsequences having significant smaller length than the original dna sequence. 
# Furthermore, we analyse the SNPs on those fragments trying to investigate which couple of restriction enzymes is giving the largest fragments with the maximum number of SNPs.

# In[371]:


# For analysis purpose we use the restriction enzyme package from Bio.Restriction
#
# load libraries

import pandas as pd
from Bio import SeqIO, Restriction
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq
import os
import sys
import seaborn as sn


# In[202]:


# Load information about restriction enzymes
renAll = Restriction.AllEnzymes


# Here are the attributes of AllEnzymes
# 
# 'add',
#  'add_nocheck',
#  'add_supplier',
#  'already_mapped',
#  'as_string',
#  'clear',
#  'copy',
#  'current_suppliers',
#  'difference',
#  'difference_update',
#  'discard',
#  'elements',
#  'format',
#  'get',
#  'intersection',
#  'intersection_update',
#  'is_restriction',
#  'isdisjoint',
#  'issubset',
#  'issuperset',
#  'lambdasplit',
#  'mapping',
#  'pop',
#  'remove',
#  'search',
#  'show_codes',
#  'split',
#  'suppl_codes',
#  'suppliers',
#  'symmetric_difference',
#  'symmetric_difference_update',
#  'union',
#  'update'
#  

# In[203]:


#number of restriction enzymes
print("Restriction enzymes are",len(renAll.elements()))
for e in renAll.elements()[:10]:
    print(e, end = " ")
print(" and so on...")


# In[204]:


#attributes of enzyme class
for at in dir(renAll.get("BcoDI"))[59:]:
    print(at, end = "   ")


# In[206]:


# Read additional information about restriction enzymes from xlsx file
infofile = pd.read_csv("/home/loukas/Documents/bioinformatics/ddrad RE double digest.xlsx - RE characteristics.tsv", sep="\t")
infofile


# In[207]:


infofile.columns


# In[208]:


palindromic_indexes = infofile["Palindromic enzymes"] == "Yes"

# Select palindromic enzymes
restr = infofile.loc[palindromic_indexes, "RE"]

print("These are the {} palindromic enzymes\n".format(len(restr)))
for e in restr:
    print(e, end="  ")


# In[209]:


# Extract recognition sites of palindromic enzymes
restrsites = infofile.loc[palindromic_indexes, "recognition site /cut"]
restrsites = restrsites.str.replace("/", "")
restrsites = pd.Series(restrsites.values, index=infofile.loc[palindromic_indexes, "RE"]); restrsites


# In[210]:


# Keep all the palindromic enzymes in RestrictionBatch class

from Bio.Restriction import *

renz = RestrictionBatch([])

for i in restr.index:
    print(i, restr[i], end="  ")
    if i not in [39,205]:
        renz.add(restr[i])

print("\n\nThere are {} palindromic restriction enzymes".format(len(renz)))


# In[283]:


# Categorize recognition sites into "frequent" if len ==4 otherwise "rare" 
frequent = RestrictionBatch([])
rare = RestrictionBatch([])

for e in renz:
    if e.size == 4:
        frequent.add(e)
    else:
        rare.add(e)

print("Enzymes frequent : {} \n{}".format(len(frequent),frequent))
print("\nEnzymes rare : {} \n{}".format(len(rare),rare))


# In[314]:


# Load 19 chromosomes of Vitis Vinifera
file = "/home/loukas/Documents/bioinformatics/Vitis_vinifera.fa.gz"
seqs = list(SeqIO.parse(file, "fasta"))

for s in seqs:
    print("Chromosome {}, size {}".format(s.id,sys.getsizeof(s)))

print("\nTotal size:",sys.getsizeof(seqs))


# In[18]:


# attributes and methods for dna sequence
for at in dir(seqs[0])[38:]:
    print(at,end="  ")


# In[212]:


# fcut the chromosome with a pair of enzymes
# returns the positions of the double digest fragments

def double_digest(chromosome,enz1, enz2):
    enzpos1 = enz1.search(chromosome.seq)
    enzpos2 = enz2.search(chromosome.seq)
    ddpos = []; start = []; end = []
    i=0; j=0
    
    if enzpos1[0] < enzpos2[0]:
        start = enzpos1
        end = enzpos2
    else:
        start = enzpos2
        end = enzpos1

    while i < len(start) and j < len(end):
        while i < len(start) and start[i] <= end[j]:
            i+=1
        
        if i == len(start):
            ddpos.append(start[i-1])
            ddpos.append(end[j])
            return ddpos
        
        ddpos.append(start[i-1])
        ddpos.append(end[j])
        
        while  j < len(end) and  end[j] <= start[i]:
            j+=1
        
        if j == len(end):
            ddpos.append(end[j-1])
            ddpos.append(start[i])
            return ddpos
        
        ddpos.append(end[j-1])
        ddpos.append(start[i])
        
    return ddpos


# We cut each chromosome with all the possible couples of restriction enzymes and store them to fasta file

# In[338]:



def ddcut():
    os.chdir("/home/loukas/Documents/bioinformatics")

    for ch_num in range(1,len(seqs)):
        ch = seqs[ch_num]
        path = "chromosome" + str(ch.id)
        os.mkdir(path)
        os.chdir(path)
        for i in range(len(renz.elements()) - 1):
            e1 = renz.get(renz.elements()[i])
            for j in range(i+1, len(renz.elements())):
                e2 = renz.get(renz.elements()[j])
                dd = double_digest(ch,e1,e2)

                #write to file
                filename = "chromosome" + ch.id + str(e1) + str(e2) + ".fasta"
                with open(filename, "w") as file1:
                    for k in range(0,len(dd),2):
                        file1.write(" ".join([">", str(e1), str(e2), str(k), "chromosome_id", ch.id, "length", str(len(ch.seq[dd[k]:dd[k+1]])), "\n"]))
                        file1.write(str(ch.seq[dd[k]-1:dd[k+1]-1]) + "\n")

        os.chdir("/home/loukas/Documents/bioinformatics")

#ddcut()


# In[339]:


# function that yields every line in fasta file without storing in RAM
import os
import os.path

def read_fasta(chromosome, enz1,enz2):
    os.chdir("/home/loukas/Documents/bioinformatics/chromosome"+str(chromosome.id))
    
    filename = "chromosome"+str(chromosome.id)+str(enz1)+str(enz2)+".fasta"
    
    if os.path.exists(filename)==False:
        filename = "chromosome"+str(chromosome.id)+str(enz2)+str(enz1)+".fasta"
    
    with open(filename, "r") as handle:
        for line in handle:
            yield line


# In[446]:


# matrices for stats

resMat = pd.DataFrame(0, index=renz.elements(), columns=renz.elements())

# number of fragments for each couple
resMatNumFrag = resMat.copy()
# summation of cut lengths for each combination
resMatLengths = resMat.copy()
# mean length values foreach combination and frequent/rare value
resMeanLengths = pd.DataFrame(index = resMat.index, columns=["MeanLength","std","frequent"])
# length fragments for each couple
resMatDistr = pd.DataFrame(list(), index=renz.elements(), columns=renz.elements())


# In[447]:


# read the .fasta file from chromosome

def build_matrices(chrom_num):
    os.chdir("/home/loukas/Documents/bioinformatics/")
    ch = seqs[chrom_num]
    mnl = []
    std = []
    frq = []

    for i in range(len(renz.elements())):
        e1 = renz.elements()[i]
        for j in range(len(renz.elements())):
            e2 = renz.elements()[j]
            num_frag = 0
            len_frag = []
            if e1 != e2:
                for l in read_fasta(ch,e1,e2):
                    if l[0]==">":
                        num_frag += 1
                        len_frag.append(int(l.split()[-1]))
                resMatNumFrag.loc[e1,e2] = num_frag
                resMatLengths.loc[e1,e2] = sum(len_frag)
                resMatDistr.loc[e1,e2] = len_frag
            else:
                continue
                
    
    for e in resMatLengths.columns:
            mnl.append(resMatLengths[e].mean())
            std.append(resMatLengths[3].std())

            if len(renz.get(e).site)==4:
                frq.append("freq")
            else:
                frq.append("rare")

    resMeanLengths["MeanLength"] = mnl
    resMeanLengths["std"] = std
    resMeanLengths["frequent"] = frq


build_matrices(chrom_num=1)


# In[448]:


resMatLengths


# In[449]:


resMatNumFrag


# In[468]:


resMatDistr = resMatDistr.fillna(0)
resMatDistr.iloc[:5,:10]


# In[451]:


resMeanLengths = resMeanLengths.sort_values(by="frequent")
resMeanLengths


# In[458]:


# find enzymes with 200 < fragment < 600 ------ let's call it fine cut!

def fineCut(dataframe):
    
    d = pd.DataFrame(0, index = dataframe.index, columns = dataframe.columns)
    
    for i in range(len(dataframe.index)):
        for j in range(len(dataframe.columns)):
            if dataframe.loc[dataframe.index[i],dataframe.columns[j]] != 0:
                count=0
                for l in dataframe.loc[dataframe.index[i],dataframe.columns[j]]:
                    if l >=200 and l <=600:
                        count+=1
                d.loc[dataframe.index[i],dataframe.columns[j]] = count
    
    return d
    
finecut = fineCut(resMatDistr)
finecut


# In[455]:


# histogram for fine cuts
rareFreqColor = {}
for e in resMeanLengths.index:
    rareFreqColor[e] = resMeanLengths.loc[e,"frequent"]
    
colors = ["green"  if rareFreqColor[v]=="freq" else "purple" for v in finecut.index]


# Create subplots with specified colors
fig, axes = pl.subplots(11, 6, figsize=(15, 15), tight_layout=True)

k=0
for i in range(11):
    for j in range(6):
        col = finecut.columns[k]
        axes[i,j].hist(finecut[col], bins=10, color=colors[k], log=True)
        axes[i,j].set_title(f'{col}')
        k+=1
        if k == len(finecut.columns):
            break

pl.suptitle('Fine Cuts for Enzymes\n')
pl.tight_layout()
pl.title("Fine Cuts")
pl.savefig("/home/loukas/Documents/bioinformatics/fineCutsHist.pdf",format="pdf",dpi=100,orientation='portrait')
pl.show()


# In[456]:



#sort lengths matrix by rare/freq
resMatLengths = resMatLengths.reindex(resMeanLengths.index)
resMatLengths = resMatLengths[resMeanLengths.index]
lengthsSimple = resMatLengths.copy()

for i in range(len(lengthsSimple.index)-1):
    for j in range(i+1,len(lengthsSimple.index)):
        lengthsSimple.iloc[i,j] = 0

def heatmap_(dataframe,fname):
    pl.figure(figsize=(18,12))
    hm = sn.heatmap(dataframe, cmap="Spectral",xticklabels="auto", yticklabels="auto",linewidths=0)
    for lb in hm.get_xticklabels():
        if resMeanLengths.loc[lb.get_text(),"frequent"]=="freq":
            lb.set_color("darkgreen")
    for lb in hm.get_yticklabels():
        if resMeanLengths.loc[lb.get_text(),"frequent"]=="freq":
            lb.set_color("darkgreen")

    pl.tight_layout()
    pl.title("Fragment Lengths\n")
    pl.savefig("/home/loukas/Documents/bioinformatics/"+fname+".pdf",format="pdf",dpi=100,orientation='portrait')

heatmap_(lengthsSimple,"fraglen")


# In[460]:


# finecut with zeroes above diagonial

finecutSimple = finecut.copy()

for i in range(len(finecut.index)):
    for j in range(i,len(finecut.columns)):
        finecutSimple.iloc[i,j] = 0
        
heatmap_(finecutSimple,"frag200_600")


# In[ ]:


# Top-10 enzymes with the largest cuts



# In[ ]:





# Observing from another perspective, if we sort the resMatLengths by column we could see areas of high and low density

# In[366]:


# sort resMatLengths by column

sort_resMatLengths = resMatLengths.sort_values(by=list(resMatLengths.columns))

heatmap_(sort_resMatLengths)


# In[463]:


# scatter for mean lengths

pl.figure(figsize=(15,3))
color_mapping = {'freq':"green","rare":"blue"}
colors = [color_mapping[cut] for cut in resMeanLengths["frequent"]]
pl.scatter(x=resMeanLengths.index,y=resMeanLengths["MeanLength"],c=colors, label='freq', alpha=0.9)
pl.grid(True)
pl.title("Mean Length cut for enzymes\n")
pl.xticks(rotation=90)
pl.legend()
pl.tight_layout()

pl.savefig("/home/loukas/Documents/bioinformatics/scatterLengths.pdf",format="pdf",dpi=600)
pl.show()


# We observe that frequent enzymes have low mean length whereas rare enzymes have high, which is expected. In contrast, two frequent enzymes seems to sut the chromosome in rather large fragments (HinP1I, HhaI) and three rare enzymes in small fragments (AsiSI, AscI, FseI)

# In[281]:


sort_resMatLengths


# In[464]:


resMatLengths.hist(bins=10, grid=True,sharey=True,figsize=(40,20))
pl.tight_layout()
pl.savefig("/home/loukas/Documents/bioinformatics/hist_all.pdf",format="pdf",dpi=1600)

pl.show()


# In[163]:


#enzymes with much small cut lengths
print(renz.get(NotI).site)
print(renz.get(AsiSI).site)
print(renz.get(NotI).site)


# In[161]:


#enzymes with large cut lengths
print(renz.get(BamHI).site)


# In[170]:


# Find cases without a conflict
for enzyme1 in frequent:
    for enzyme2 in rare:
        if enzyme1.site in enzyme2.site:
            resMatConf.loc[str(enzyme1), str(enzyme2)] = 1
        else:
            resMatConf.loc[str(enzyme1), str(enzyme2)] = 0
        

resMatConf

