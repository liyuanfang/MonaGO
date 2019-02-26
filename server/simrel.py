import networkx as nx
import json
import math
import unicodedata
import csv
from Bio.UniProt.GOA import gafiterator
import gzip
import time


def makeTree():
    # makes gene ontology tree using networkx and data from GO.js
    G=nx.DiGraph()
    for key in go_hier.keys():
        if 'GO' in unicodedata.normalize('NFKD',key).encode('ascii'):
            G.add_node(key)
    for key in G.nodes():
        parents=go_hier[key]['p']
        for parentGO in parents:
            G.add_edges_from([(parentGO,key)])
    return(G)

with open('js/GO.js','r') as fr_GO:
    for GO in fr_GO:
        go_hier=json.loads(str(GO))
G=makeTree()


PDictFile= open('PDict.json').read()
PDict=json.loads(PDictFile)


def simRel(nodeA,nodeB):
    # returns the simrel of two GO IDs
    ancestorsNodeA=nx.ancestors(G,nodeA)
    ancestorsNodeA.add(nodeA)
    ancestorsNodeB=nx.ancestors(G,nodeB)
    ancestorsNodeB.add(nodeB)
    commonAncestors=ancestorsNodeA.intersection(ancestorsNodeB)
    if len(commonAncestors)==0:
        return(0,'')
    if len(commonAncestors)==0:
        return(0,'')
    maxIC=0
    MICA=''
    for ancestor in commonAncestors:
        try:
            ICanc=-math.log10(PDict[ancestor])
        except ValueError:
            return(1,ancestor)
        if ICanc>=maxIC:
            maxIC=ICanc
            MICA=ancestor
    try:
        numerator=2*math.log10(PDict[MICA])
    except ValueError:
        return(1,MICA)
    try: 
        denominator=math.log10(PDict[nodeA])+math.log10(PDict[nodeB])
    except ValueError:
        return(0,MICA)
    simrel=numerator * (1-PDict[MICA])/denominator
    return(simrel,MICA)







def SimRelDis(iGOids, x, yGOids, clustComp):
    '''
    calculate the distance between clusters based on the average resnik distance between each GO
    '''
    simrelDistances=[]
    for iGOid in iGOids:
        for yGOid in yGOids:
            simrelDistances.append(simRel(iGOid,yGOid)[0])
    if clustComp=='ave':
        averageSimRel=float(sum(simrelDistances))/float(len(simrelDistances))
        return(averageSimRel)
    elif clustComp=='min':
        minSimRel=float(min(simrelDistances))
        return(minSimRel)
    elif clustComp=='max':
        maxSimRel=float(max(simrelDistances))
        return(maxSimRel)



def SimRelDisPerc(iGOids,x,yGOids, clustComp):
    scalar=1#-math.log10(float(1)/float(maxFreq))
    normalised=SimRelDis(iGOids,x,yGOids, clustComp)/scalar
    # use cuberoot to adjust scale making it slightly more linear... might need to think about using a different scalar or writing a function that checks for linearity.
    return((normalised**(1.)*100))





def SimRelDisMIA(iGOids, x, yGOids):
    sim=0
    resnikDistances=[]
    for iGOid in iGOids:
        for yGOid in yGOids:
            sim,MIA=simRel(iGOid,yGOid)
    MIAName=go_hier[MIA]['n']
    return(sim,MIAName)

