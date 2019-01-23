import json
import random
import csv
from resnik import simRes
from simrel import simRel

def makeRandomList():

    with open('js/GO.js','r') as fr_GO:
        for GO in fr_GO:
            go_hier=json.loads(str(GO))


    sampleList=random.sample(go_hier.keys(),50)


    sampleString=''
    for GOTerm in sampleList:
        sampleString+=GOTerm.encode('ascii')+' '

    print(sampleString)



#makeRandomList()

def readTable():
    with open('table.csv') as f:
        data=[]
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            data+=[[row['GO term 1'],row['GO term 2'],float(row['Res']), float(row['simRel'])]]
    return(data)



table=readTable()

errorRes=[]
errorRel=[]
for row in table:
    GO1=row[0]
    GO2=row[1]
    funsimmatRes=row[2]
    mysimRes,MICA=simRes(unicode(GO1),unicode(GO2))
    if funsimmatRes!=0:
        errorRes+=[(funsimmatRes-mysimRes)/funsimmatRes]
    elif funsimmatRes==0 and mysimRes==0:
        errorRes+=[0]
    elif funsimmatRes==0:
        print('hi')
        print(GO1)
        print(GO2)
        errorRes+=[1]


    funsimmatRel=row[3]
    mysimRel,MICA=simRel(unicode(GO1),unicode(GO2))
    if funsimmatRel!=0:
        errorRel+=[(funsimmatRel-mysimRel)/funsimmatRel]
    elif funsimmatRel==0 and mysimRel==0:
        errorRel+=[0]
    elif funsimmatRel==0:
        print('hi')
        print(GO1)
        print(GO2)
        errorRel+=[1]

print(sum(errorRes)/len(errorRes))
print(sum(errorRel)/len(errorRel))
