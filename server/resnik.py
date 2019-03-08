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


ICDictFile= open('ICDict.json').read()
ICDict=json.loads(ICDictFile)


def simRes(nodeA,nodeB):
    # returns the Resnik similarity of two GO IDs
    ancestorsNodeA=nx.ancestors(G,nodeA)
    ancestorsNodeA.add(nodeA)
    ancestorsNodeB=nx.ancestors(G,nodeB)
    ancestorsNodeB.add(nodeB)
    commonAncestors=ancestorsNodeA.intersection(ancestorsNodeB)
    if len(commonAncestors)==0:
        return(0,'')
    maxIC=0
    MICA=''
    for ancestor in commonAncestors:
        if type(ancestor)=='unicode':
            ancestor=unicodedata.normalize('NFKD',ancestor).encode('ascii')
        ICanc=ICDict[ancestor]
        if ICanc>=maxIC:
            maxIC=ICanc
            MICA=ancestor
    return(maxIC, MICA)







def ResnikDis(iGOids, x, yGOids, clustComp):
    '''
    calculate the distance between clusters based on the average resnik distance between each GO
    '''
    resnikDistances=[]
    for iGOid in iGOids:
        for yGOid in yGOids:
            resnikDistances.append(simRes(iGOid,yGOid)[0])
    if clustComp=='ave':
        averageResnik=float(sum(resnikDistances))/float(len(resnikDistances))
        return(averageResnik)
    elif clustComp=='min':
        minResnik=float(min(resnikDistances))
        return(minResnik)
    elif clustComp=='max':
        maxResnik=float(max(resnikDistances))
        return(maxResnik)




def ResnikDisPerc(iGOids,x,yGOids, clustComp):
    scalar=1
    normalised=ResnikDis(iGOids,x,yGOids, clustComp)/scalar
    # use cuberoot to adjust scale making it slightly more linear... might need to think about using a different scalar or writing a function that checks for linearity.
    return((normalised**(1.)*100))





def ResnikDisMIA(iGOids, x, yGOids):
    sim=0
    resnikDistances=[]
    for iGOid in iGOids:
        for yGOid in yGOids:
            sim,MIA=simRes(iGOid,yGOid)
    MIAName=go_hier[MIA]['n']
    return(sim,MIAName)

print(simRes(u'GO:0006355',u'GO:0006351'))
nodeA=u'GO:0006355'
nodeB=u'GO:0006351'
if True:
    ancestorsNodeA=nx.ancestors(G,nodeA)
    ancestorsNodeA.add(nodeA)
    ancestorsNodeB=nx.ancestors(G,nodeB)
    ancestorsNodeB.add(nodeB)
    commonAncestors=ancestorsNodeA.intersection(ancestorsNodeB)
    print(ancestorsNodeA)
    print(ancestorsNodeB)
