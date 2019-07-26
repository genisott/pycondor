#!/usr/bin/env python
# coding: utf-8

# In[28]:


import numpy as np
import pandas as pd
import time
from igraph import *


# In[29]:


def matrices(net,genes,regs):
    """Builds the adjacency matrix and the modularity matrix."""
    # Building matrix A
    ls = {(genes[net.iloc[i,0]],regs[net.iloc[i,1]]):net.iloc[i,2] for i in range(0,len(net))}
    A = np.matrix(np.zeros((len(genes),len(regs))))
    for i in ls:
        A[i]=ls[i]
      
    ki = A.sum(1)
    dj = A.sum(0)   
    m = float(sum(ki))
    B = A - ((ki@dj)/m)
    
    
    return A,B,m


# In[34]:


def bipartite_modularity(B,m,R,T):  
    """Computation of the bipartite modularity as described in ""Modularity and community detection in bipartite networks" by Michael J. Barber." """
    RtBT = R.transpose().dot(B.dot(T))
    Q = (1/m)*(np.trace(RtBT))
    return Q


# In[31]:


def max_modularity(B,m,R0,T0,c,deltaQmin):
    """Implementation of the BRIM algorithm to iteratively maximize bipartite modularity."""
    Qnow = 0
    deltaQ = 1
    p,q = B.shape
    while(deltaQ>deltaQmin):
        #Right sweep
        Tp = B.dot(T0)
        R = np.zeros((p,c))
        am = np.array(np.argmax(Tp, axis=1))
        for i in range(0,len(am)):
            R[i,am[i][0]]=1
        #Left sweep
        Rp = B.transpose().dot(R)
        T = np.zeros((q,c))
        am = np.array(np.argmax(Rp, axis=1))
        for i in range(0,len(am)):
            T[i,am[i][0]]=1
        T0 = T
        
        Qthen = Qnow
        Qnow = bipartite_modularity(B,m,R,T)
        deltaQ = Qnow-Qthen
        print(Qnow)
    return R,T    


# In[32]:


def initial_communities(A,gn,rg,c):
    """Uses multilevel community detection from the package python-igraph to create the initial community structure."""
    
    p,q = len(gn),len(rg)
    
    Am = np.zeros((p+q,p+q))
    Am[:p, p:] = A
    Am[p:,:p] = A.transpose()
    
    grn = Graph.Weighted_Adjacency(Am.tolist(), mode=ADJ_UNDIRECTED, attr="weight")
    n = p+q
    vc = Graph.community_multilevel(grn,weights="weight")
    
    
    M = np.zeros((n,c))
    for i in range(0,len(vc)):
        for j in vc[i]:
            M[j,i]=1
    R,T = M[:len(gn),],M[len(gn):,]
    
    return R,T


# In[35]:


def condor(filename,c,deltaQmin="def"):
    """Main function of the pyCONDOR package.
    This function reads a filename which should be a tab-separated file containing edges and weights for a bipartite graph.
    c is the max number of communities that can be found.
    deltaQmin is the tolerance for the iterative maximization of the modularity.
    
    IMPORTANT: The network MUST be bipartite and connected, otherwise the results are meaningless.
    
    The function outputs two files with the communities to which the vertices of each of the bipartite sets belong. It also 
    returns the matrices Rf,Tf of membership to communities.    
    
    This method uses the algorithm descibed in "Modularity and community detection in bipartite networks" by Michael J. Barber.
    and is somewhat translated from the R package CONDOR https://github.com/jplatig/condor described in the paper
    "Bipartite Community Structure of eQTLs" by John Platig , Peter J. Castaldi, Dawn DeMeo, John Quackenbush. 
    
    This python version is a bit more limited than the R package in terms of options.
    
    """
    tt = time.time()
    #Loads the network file.
    t = time.time()
    net = pd.read_csv(filename,sep="\t",index_col=0)
    net = net[["tar","reg","weight"]]
    
    if(deltaQmin=="def"):
        deltaQmin = min(1/len(net),1e-4)
    
    gn = list(sorted(set(net["tar"])))
    genes = {gn[i]:i for i in range(0,len(gn))}
    rg = list(sorted(set(net["reg"])))
    regs = {rg[i]:i for i in range(0,len(rg))}
    print("Data loaded.",time.time()-t)
    
    #Builds matrices A (tilde), B (tilde):
    t = time.time()
    A,B,m = matrices(net,genes,regs)
    print("Matrices computed.",time.time()-t)
    
    #Computes initial community T0:
    t = time.time()
    R0,T0 = initial_communities(A,gn,rg,c)
    print("Initial communities computed.",time.time()-t)
    
    #Maximize modularity:
    t=time.time()
    Rf,Tf = max_modularity(B,m,R0,T0,c,deltaQmin)
    print("Modularity maximized.",time.time()-t)
    
    #Output to file:
    t = time.time()
    _,j = np.where(Rf>0)
    pd.DataFrame({"gene":gn,"comm":list(j)}).to_csv("com_gene_"+filename+".txt")
    _,j = np.where(Tf>0)
    pd.DataFrame({"reg":rg,"comm":list(j)}).to_csv("com_reg_"+filename+".txt")
    print("Communities exported.",time.time()-t)
    
    print("Total time:",time.time()-tt)
    return Rf,Tf


# In[ ]:




