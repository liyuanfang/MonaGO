import json
import networkx as nx
import unicodedata
import gzip
from Bio.UniProt.GOA import gafiterator




def createGOjsFile():
    # recreates GO.js file so that ontology is up to date
    import json
    from goatools import obo_parser

    file=obo_parser.GODag('data/go-basic.obo','relationship')

    GOjsDict={}
    for goID in file:
        parents=[]
	uppers=[]
	for parent in file[goID].parents:
	    parents.append(parent.id)
	    uppers.append((parent.id,'is_a'))

	if len(list(file[goID].relationship))==1:
            for relationships in list(file[goID].relationship.values()):
                relation=list(relationships)[0].id
	    parents.append(relation)
	    uppers.append([relation,list(file[goID].relationship)[0]])
	if len(list(file[goID].relationship))>1:
	    i=0
            for relationships in list(file[goID].relationship.values()):
                relation=list(relationships)[0].id
	        parents.append(relation)
	        uppers.append([relation,list(file[goID].relationship)[i]])
		i+=1
        name=file[goID].name
        namespace=file[goID].namespace
        if not file[goID].is_obsolete:
            GOjsDict[unicode(file[goID].id)]={'p':parents,'c':namespace,'n':name, 'u':uppers}
    GOjs=open('js/GO.js','w')
    json.dump(GOjsDict,GOjs)


createGOjsFile()







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


def getAnnotations():
    # uses database from uniprot to count the number of annotations each GO term has, and returns a dictionary of GO IDs and their number of annotations
    annotDict={}
    for goID in G.nodes():
        annotDict[goID]=0
    filename = 'data/goa_uniprot_all.gaf.gz'
    with gzip.open(filename, 'rt') as fp:
        count=0
        for annotation in gafiterator(fp):
            count+=1
            percentage=float(count)*100/424606000
            if percentage>5:
                print('hi')
                count=0
            goID=unicode(annotation['GO_ID'])
            try:
                annotDict[goID]+=1
            except KeyError:
                print(goID)
    return annotDict



def createAnnotationDict():
    # creates annotation dictionary as it takes a while to go through 8GB file from uniprot
    getCount=getAnnotations()
    with open('annotFile.json','w') as annotFile:
        json.dump(getCount,annotFile)


createAnnotationDict()
