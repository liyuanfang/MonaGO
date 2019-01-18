# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from HTMLParser import HTMLParser
import re
import requests
import threading
import logging

from logTime import logTime

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

class myThread (threading.Thread):
    def __init__(self, threadID, request, url):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.url = url
        self.request = request

    def run(self):
        self.res = self.request.get(self.url,verify=False)


    def getRes(self):
        return self.res.content

class DavidDataScrawler(object):

    def setParams(self,inputIds,idType,annotCat,pVal):
        self.inputIds,self.idType,self.annotCat,self.pVal = inputIds,idType,annotCat,pVal

    def run(self):

        s = requests.session()
        #s.cert = '/path/client.cert'

        #d = Deferred()#init defer

        res = self._uploadGene(s,self.idType,self.inputIds)

        with open('test.html',"w") as fr_html:
            fr_html.write(res) 

        if self._checkSuccess(res):
            logger.debug("upload gene success")

            url_1 = 'https://david.ncifcrf.gov/chartReport.jsp?annot={0}&currentList=0'.format(self.annotCat)
            url_2 = 'https://david.ncifcrf.gov/list.jsp'

            thread1 = myThread(1, s, url_1)
            thread2 = myThread(2, s, url_2)

            thread1.start()
            thread2.start()

            threads = []
            threads.append(thread1)
            threads.append(thread2)

            for t in threads:
                t.join()


            try:
                getGO_response = thread1.getRes()
                geneList_response = thread2.getRes()
            except:
                raise Exception("get response from david failed") ##unlikely to happen but good to have

            # with open('test1.html',"w") as fr_html:
            #     fr_html.write(getGO_response) 

            go,geneIds = self._parseGO(getGO_response, s)

            geneList = self._parseGenes(geneList_response)

                
            # with open("go.txt","w") as fw:
            #     fw.write(str(go))

            go_filtered = self._filterGO(self.pVal,go)

            # with open("go_filtered.txt","w") as fw:
            #     fw.write(str(go_filtered))

            geneIds = self._getUniqueGeneIds(geneIds)


            #logger.debug("geneIds:{}".format(geneIds))

            geneIdNameMapping = self._getGenesNamesByIds(geneIds,geneList)

            #change the gene id into gene name in go
            go_processed = self._changeGeneIdToNameInGO(go_filtered,geneIdNameMapping)


            if not go_processed:
                raise Exception("get go_processed failed")
            

            #before return, close the connection to DAVID explicitly

            return go_processed


        else:
            logger.info("get chartReport failed")
            raise Exception("upload genes to DAVID failed")


    def _checkSuccess(self,res):
        if res.find("DAVID: Functional Annotation Tools")==-1:
            return False
        else:
            return True


    def _parseGO(self,getGO_response,request):
        parser = DavidDataScrawler.GOParser(request)
        parser.feed(getGO_response)#get go
        go = parser.getGO_inf()
        geneIds = parser.getGeneLists()

        #logger.debug('go:{}'.format(go))

        if not go:
            raise Exception("get GO terms failed") 

        #logger.debug('geneIds:{}'.format(geneIds))

        if not geneIds:
            raise Exception("get gene lists failed") 

        return go,geneIds


    def _parseGenes(self,geneList_response):
        parser = DavidDataScrawler.geneParser()
        parser.feed(geneList_response)
        genesIdName = parser.getParsedGeneNameWithId()

        if len(genesIdName) == 0:
            raise Exception("parse gene list failed") 

        return genesIdName




    @logTime
    def _uploadGene(self,s,idType,inputIds):

        data = {'idType':(None, idType), 'uploadType': (None, 'list'), 'multiList': (None,'false'),'Mode':(None,'paste'),
                         'useIndex':(None,'null'),'usePopIndex': (None,'null'),'demoIndex': (None,'null'),'ids': (None,inputIds),
                         'removeIndex': (None,'null'),'renameIndex': (None,'null'),'renamePopIndex': (None,'null'),
                         'newName': (None,'null'),'combineIndex': (None,'null'),'selectedSpecies': (None,'null'),'uploadHTML': (None,'null'),
                         'managerHTML': (None,'null'), 'sublist': (None,''),'rowids': (None,''),'convertedListName': (None,'null'),
                         'convertedPopName': (None,'null'),'pasteBox': (None,inputIds),'Identifier': (None,idType) ,
                         'selectedSpecies': (None,',0') , 'speciesList': (None,'0'), 'rbUploadType':(None,'list')}

        r = s.post('https://david.ncifcrf.gov/tools.jsp', files=data,  verify=False)

        return r.content




    def _changeGeneIdToNameInGO(self,go,geneIdNameMapping):
        '''
        transfrom the gene ids in GO terms to gene name according to the mapping

        Args:
            go:a list of GO terms
            geneIdNameMapping:mapping between gene id and gene name

        return:
            processed GO terms
        '''

        for i in range(0,len(go)):

            geneNames = go[i]["genes"].split(",")
            geneNameString = ";".join([geneIdNameMapping[geneName.strip().lower()] for geneName in geneNames]) 

            go[i]["genes"] = geneNameString[:-1]# strip the last '|'

        return go


    def _getGenesNamesByIds(self,geneIds,genesIdName):
        '''
        create a dict(geneIdNameMapping) to map gene id to gene name

        Args:
            pcHelper: pycurl helper object
            geneIds: an array of gene ids

        return:
            a dict mapping gene id to gene name
        '''

        geneIdNameMapping = {}

        for i in range(0,len(genesIdName)-1,2):#ugly way, need improve
            geneIds = genesIdName[i].split(",")
            if len(geneIds) > 1:#if multiple gene ids
                for j in geneIds:
                    geneIdNameMapping[j.strip().lower()] = genesIdName[i+1]
            else:
                geneIdNameMapping[genesIdName[i].lower()] = genesIdName[i+1]

        return geneIdNameMapping

    def _getUniqueGeneIds(self,geneListsDic):
        rowids = set([])
        rowidstr = ""

        #map(rowids.add, [genes for genes in geneListsDic.values()]) #list is not hashable
        for genes in geneListsDic.values():
            if isinstance(genes,list):
                map(rowids.add,[gene for gene in genes])
            else:
                rowids.add(genes)

        rowidstr = ','.join(rowids)

        return rowidstr



    def _filterGO(self,pVal,go_inf):

        filterGO_inf = [GO for GO in go_inf if float(GO["pVal"]) < float(pVal)]

        #logger.debug('go_filtered:{}'.format(filterGO_inf))

        if not filterGO_inf:
            raise Exception("get go_filtered failed")

        return filterGO_inf


    class GOParser(HTMLParser):
        '''
        parse a list of GO terms(name and a list of associated gene ids)
        '''

        def __init__(self,request):
            HTMLParser.__init__(self)

            self.request = request
            self.go_inf = []
            self.geneLists = {}
            self.metacount = 0

        def handle_starttag(self, tag, attrs):

            #get GO_id,GO_name,p-value,count
            if tag == "a":
                m = re.search('(data/download/chart_\w+.txt)',attrs[0][1])
                
                if m!=None:
                #if exists
                    url = 'https://david.ncifcrf.gov/'+m.group(0)
                    res = self.request.get(url,verify = False)
                    self._parseGO(res.content)

            #get gene rowid
            if tag == "img":
                if attrs[0][1] == 'graphics/two_tone_2_a.jpg':
                    genes = attrs[6][1].split(";")[1]
                    self.geneLists[self.metacount] = self._parseStringIntoList(genes)
                    self.metacount+=1

        def _parseGO(self,GO):
            line = GO.encode('ascii','ignore').split("\n")
            for i in range(1,len(line)-1):
                cell = line[i].split("\t")
                GO_term = cell[1].split("~")
                self.go_inf.append({"cat":cell[0],"GO_id":GO_term[0],"GO_name":GO_term[1],"count":cell[2],"pVal":cell[4],"genes":cell[5].lower().strip()})

        def getGO_inf(self):
            return self.go_inf

        def getGeneLists(self):
            return self.geneLists

        def _parseStringIntoList(self,genes):
            index = genes.find("geneReport")
            geneList = genes[index+10:][2:-2].encode('ascii','ignore').split(",")
            return geneList


    class geneParser(HTMLParser):
        '''
        parse gene id and gene name
        '''

        table = 0
        tr = 0
        td = 0
        genesIdName = []
        # function to handle an opening tag in the doc

        def handle_starttag(self, tag, attrs):
            if tag == "table":
                if len(attrs) > 0:
                    if attrs[0][1] == "dataTable":
                        self.table = 1

            if self.table == 1 and tag == "tr":
                if len(attrs) > 0:
                    if attrs[0][0] == "class":
                        self.tr = 1
            if self.table == 1 and self.tr == 1 and tag == "td":
                self.td+=1

        def handle_endtag(self, tag):
            if tag == "table":
                self.table = 0

            if tag == "tr" and self.table == 1:
                self.tr = 0

            if self.table == 1 and self.tr == 0:
                self.td = 0


        def handle_data(self, data):

            if self.tr == 1:
                if data != "RG" and self.td!= 4 and "\n" not in data:
                    self.genesIdName.append(data.encode('ascii','ignore'))

        def getParsedGeneNameWithId(self):
            return self.genesIdName



    # class sessionIdParser(HTMLParser):
    #     def __init__(self):
    #         HTMLParser.__init__(self)
    #         self.sessionId = 0

    #     def handle_starttag(self, tag, attrs):
    #         if tag == 'input':
    #             if attrs[1][0] == 'name' and attrs[1][1] == 'SESSIONID':
    #                 self.sessionId = attrs[2][1]

    #     def returnSessionId(self):
    #         return self.sessionId