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
import json

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

class MetaScapeScrawler(object):
    URLWSGI = "http://metascape.org/gp_server"

    def createNewSession(self):
        url = self.URLWSGI + "/get_session_id"
        r = requests.get(url)
        self.session_id = r.content;
        return self.session_id

    def guessIDType(self):
        url = self.URLWSGI + "/id_type_guessing"
        inputId_array = self.inputIds.split(",")
  
        data = {"data":inputId_array}
        try:
            r = requests.post(url, json = data)
        except:
            logger.error("unable to guess id type")
            raise Exception("unable to guess id type")
        else:
            self.resultType = json.loads(r.content)["result_type"]
            return self.resultType

        
    def saveTextBoxToServer(self):
        url = self.URLWSGI + "/save_textbox_to_server"
        inputId_array = self.inputIds.split(",")
        data = {"session_id":self.session_id, "ids":inputId_array}
        try:
            r = requests.post(url, json = data)
        except:
            logger.error("unable to save textbox to server")
            raise Exception("unable to save textbox to server")
        else:
            self.result = json.loads(r.content)["result"]
            return self.result

    def guessSpeciesFromInput(self):
        url = self.URLWSGI + "/guess_species_from_input"
        inputId_array = self.inputIds.split(",")
        data = {"isMultipleList":"false","session_id":self.session_id, "ids":inputId_array}
        try:
            r = requests.post(url, json = data)
        except:
            logger.error("unable to guess species from input")
            raise Exception("unable to guess species from input")
        else:
            self.tax2name = json.loads(r.content)["tax2name"]
            return self.tax2name

    def applySpecies(self,analysisSpecies, inputSpecies):
        url = self.URLWSGI + "/apply_species"
        data = {"specifiedSpeciesOption":{"session_id":self.session_id,"multipleList":{"value":"false"},"specifiedSpecies":{"analysisSpecies":analysisSpecies,"inputSpecies":inputSpecies}}}
        try:
            r = requests.post(url, json = data)
        except:
            logger.error("unable to apply species")
            raise Exception("unable to apply species")
        else:
            self.content = json.loads(r.content)["content"]
            return self.content

    def idConversion(self):
        url = self.URLWSGI + "/idconversion"
        inputId_array = self.inputIds.split(",")
        data = {"data":inputId_array,"option":{"wherePutResult":"CURRENT_SHEET","hasHeader":"false","one2many":"KEEP_FIRST_ONLY",
        "fromId":"symbol,RefSeq_Proteins,RefSeq_RNAs,gene_synonym,Gene_History,uniprot,ensembl_gene_id,ensembl_peptide_id,ensembl_transcript_id,ucsc,ensembl_array_id,bioconductor_array_id,locus_tag,dbxref",
        "toIds":["6239_gid"],"inputSpecies":"10090","session_id":self.session_id}}

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


