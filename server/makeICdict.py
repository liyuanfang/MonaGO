import networkx as nx
import json
import math
import unicodedata
import csv
from Bio.UniProt.GOA import gafiterator
import gzip


def makeTree():
    # makes gene ontology tree using networkx and data from GO.js
    with open('js/GO.js','r') as fr_GO:
        for GO in fr_GO:
            go_hier=json.loads(str(GO))
    G=nx.DiGraph()
    for key in go_hier.keys():
        if 'GO' in unicodedata.normalize('NFKD',key).encode('ascii'):
            G.add_node(key)
    for key in G.nodes():
        parents=go_hier[key]['p']
        for parentGO in parents:
            G.add_edges_from([(parentGO,key)])
    return(G)


G=makeTree()



annotDict= open('annotFile.json').read()
getCount=json.loads(annotDict)





def getRoot(node):
    # returns the root of each GO ID
    roots=[u'GO:0003674', u'GO:0008150', u'GO:0005575']
    for root in roots:
        if root in nx.ancestors(G,node):
            return(root)
        elif root==node:
            return(root)






def getFrequency(node):
    # returns the frequency of a GO ID (ie. the sum of it and its descendants annotations )
    frequency=getCount[node]
    for descendant in nx.descendants(G,node):
        frequency+=getCount[descendant]
    return(frequency)




def IC(node):
    # returns IC of a GO ID
    f_node=getFrequency(node)
    f_root=getFrequency(getRoot(node))
    p_node=float(f_node)/float(f_root)
    # use try to deal with terms with 0 frequency
    try:
        IC_node=-math.log10(p_node)
    except ValueError:
        IC_node=-math.log10(float(1)/float(f_root))
    #scale_max=math.log(2,f_root)
    #IC_uniform=IC_node/scale_max
    #return(IC_uniform)
    return(IC_node)


ICDict={}

print(G.number_of_nodes())

count=0
for node in G.nodes():
    count+=1
    ICNode=IC(node)
    ICDict[node]=ICNode
    if count==1000:
        print('hi')
        count=0

ICdictFile=open('ICDict.json','w')
json.dump(ICDict,ICdictFile)


    
