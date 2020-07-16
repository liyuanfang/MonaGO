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


import sys
import json
# import yaml

from datetime import datetime

#for praser
from flask import Flask,render_template,request,send_from_directory,Response

from DavidDataScrawler import DavidDataScrawler
import logging
import time

from logTime import logTime

# creating specific error when multiple graphs are used for semantic calculations
class SemanticsError(Exception):
        pass

logging.basicConfig(filename="debug.txt",level=logging.INFO)
logger = logging.getLogger(__name__)

config = {}

with open("config.ini","r") as fr:
    for line in fr:
        opition = line.split(":")
        config.update({opition[0]:opition[1].strip()})

if(config["remote_server"]=="true"):
    root_dir = "/home/ubuntu/monaGoWebsite/server/"
else:
    root_dir = ""

GO_dict = {}

app = Flask(__name__)
@app.route('/', methods=['POST','GET'])
def index():

    annotCatDict = {
        'GOTERM_BP':'30',#biological process direct
        'GOTERM_CC':'38',#celluar component direct
        'GOTERM_MF':'46'#melocular function direct
    }





    if request.method == 'GET':
        return render_template("index.html")
    else:





        if request.form['type'] == "manual":
            pVal = request.form['pVal']
            if request.form['inputGOs'] == "":
                go = parseInputGOsFromCSV(request.files['files2'],float(str(pVal)))
            else:
                go = parseInputGOs(request.form['inputGOs'],float(str(pVal)))

            try:
                matrix_count, array_order, go_hier, go_inf_reord, clusterHierData, simDict = processedData(go, request.form['similarity'], request.form['clustComp'])
            except SemanticsError:
                return "Cannot measure semantic similarity across multiple roots. Use a single function annotation, or measure similarity between terms using number of overlapping genes"
            myList=simDict
            if not matrix_count:
                return "Failure to process data"
            data = "<script>"+"var go_inf_input="+str(go)+";"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+str(matrix_count)+';'+'var simDict='+str(myList)+";"+"var array_order="+str(array_order)+";"\
            +"var clusterHierData="+str(clusterHierData) +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+";"+"var clustComp='"+str(request.form['clustComp'])+"';"+"var similarity='"+str(request.form['similarity'])+"'</script>"

        if request.form['type'] == "david":
            #parameters needed for querying DAVID
            if request.form['inputIds'] == "":
                inputIds = parseInputIdsFromCSV(request.files['files1'])
            else:
                inputIds = request.form['inputIds']

            idType = request.form['idType']
            annotCat = request.form['annotCat']
            pVal = request.form['pVal']
            #transform annotation name to number recognized by DAVID(e.g. GOTERM_BP_FAT to 25) .
            annotCat = ','.join([annotCatDict[cat] for cat in annotCat.split(",")])

            go,status = getDataFromDavid(inputIds,idType,annotCat,pVal)

            if status == False:
                return "Failure to get data, please make sure the identifier is correct and try again.\nOtherwise, get enriched GO terms at https://david.ncifcrf.gov/summary.jsp and then use option 1"
            try:
                matrix_count, array_order, go_hier, go_inf_reord, clusterHierData, simDict = processedData(go, request.form['similarity'],request.form['clustComp'])
            except SemanticsError:
                return "Cannot measure semantic similarity across multiple roots. Use a single function annotation, or measure similarity between terms using number of overlapping genes"
            myList=simDict
            if not matrix_count:
                return "Failure to process data"
            data = "<script>"+"var go_inf_input="+str(go)+";"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+str(matrix_count)+';'+'var simDict='+str(myList)+";"+"var array_order="+str(array_order)+";"\
            +"var clusterHierData="+str(clusterHierData) +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+";"+"var clustComp='"+str(request.form['clustComp'])+"';"+"var similarity='"+str(request.form['similarity'])+"'</script>"


        if request.form['type'] == "MonaGO":
            pVal = request.form['pVal']
            content = request.files['files3'].read()
            content_dict = json.loads(content, encoding='utf-8')
            content_dict["go_inf"] = parseInputMonagos(content_dict, float(str(pVal)))
            similarity=content_dict['similarity'].encode('ascii', 'ignore')
            clustComp=content_dict['clustComp'].encode('ascii', 'ignore')

            # content_dict = yaml.safe_load(content_dict, encoding='utf-8')
            matrix_count, array_order, go_hier, go_inf_reord, clusterHierData, simDict = processedData(content_dict["go_inf"],similarity, clustComp)
            myList=simDict

            if not matrix_count:
                return "Failure to process data"
            data = "<script>"+"var go_inf_input="+str(content_dict["go_inf"])+";"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+str(matrix_count)+';'+'var simDict='+str(myList)+";"+"var array_order="+str(array_order)+";"\
            +"var clusterHierData="+str(clusterHierData) +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+";"+"var clustComp='"+str(clustComp)+"';"+"var similarity='"+str(similarity)+"'</script>"
            #data = "<script>"+ "var size = 0"+";"+"var content ="+content+";"+"</script>"




        with open(root_dir+'templates/chord_layout.html',"r") as fr_html:
            html = "".join(fr_html.readlines())


        return data+html


@app.route('/css/<fileName>')
def getCss(fileName):
    return send_from_directory(root_dir+'css', fileName)

@app.route('/css/font-awesome-4.6.3/css/<fileName>')
def getFontAwesome(fileName):
    return send_from_directory(root_dir+'css/font-awesome-4.6.3/css/', fileName)

@app.route('/css/font-awesome-4.6.3/fonts/<fileName>')
def getFont(fileName):
    return send_from_directory(root_dir+'css/font-awesome-4.6.3/fonts/', fileName)

@app.route('/img/<fileName>')
def getImg(fileName):
    return send_from_directory(root_dir+'img',fileName)

@app.route('/js/<fileName>')
def getJs(fileName):
    return send_from_directory(root_dir+'js', fileName)

@app.route('/fonts/<fileName>')
def getFonts(fileName):
    return send_from_directory(root_dir+'fonts', fileName)

@app.route('/txt/<fileName>')
def getText(fileName):
    return send_from_directory(root_dir+'text', fileName)

@app.route('/demo')
def returnDemo():
    with open("visitors.txt",'a') as fw:
        fw.write("remote address: {}  time: {}\n".format(request.remote_addr,datetime.today()))
    with open(root_dir+"templates/chord_layout.html","r") as fr_html:
        html = "".join(fr_html.readlines())
    with open(root_dir+"demo/Data.txt","r") as fr:
        data = fr.readline()

    return data + html

@app.route('/my/img/my_logo.jpg')
def getMyLogo():
    with open("visitors.txt",'a') as fw:
        fw.write("remote address: {}  time: {}\n".format(request.remote_addr,datetime.today()))
    return send_from_directory(root_dir+'my/img','my_logo.jpg')

@app.route('/help')
def getHelp():
    return send_from_directory(root_dir+'templates','help.html')

@app.route('/help.html')
def getHelp1():
    return send_from_directory(root_dir+'templates','help.html')

@app.route('/getPic',methods=['POST','GET'])
def getPic():
    b64_string = request.form['svg']

    #b64_string += "=" * ((4 - len(request.form['png']) % 4) % 4)

    return Response(
        b64_string,
        mimetype="image/svg+xml",
        headers={"Content-disposition":
                 "attachment; filename=chart.svg"})

@app.route('/export',methods=['POST','GET'])
def export():
    string = request.form['file']

    return Response(
        string,
        mimetype="file/txt",
        headers={"Content-disposition":
                 "attachment; filename=export.monago"})

@app.route('/csv/<filename>')
def getDemoCSV(filename):
    return send_from_directory(root_dir+'csv', filename)

def parseInputGOsFromCSV(file,pVal):
    data = file.read()
    return parseInputGOs(data,pVal)

def parseInputIdsFromCSV(file):
    data = file.read()
    inputId = []
    for line in data.split("\n"):
        if line.strip("\n\r") != "":
            if line == "geneID":
                continue 
            inputId.append(line)
    return ','.join(inputId)



def loadGOHier():

    if len(GO_dict) == 0:
        with open(root_dir+"data/GO.go") as fr:
            for line in fr:
                go_inf = line.split(",",1)
                GO_dict.update({go_inf[0]:go_inf[1]})

def getGONameAndCatergory(GO_id):
    try:
        GO_inf = GO_dict[GO_id]
    except:
        logger.error("could not find GO id in GO_dict")
        return ""
    else:
        return GO_inf

def parseInputGOs(go_csvFormat,pVal):
    #GO_id,p-value,genes

    loadGOHier()

    goDictContainer = []
    lines = go_csvFormat.split("\n")

    for line in lines:
        if line.strip("\n\r") != "":
            cols = line.split(",",2) # do not split genes
            if cols[0] == "GO_id":
                continue
            genes = cols[2].split(";")
            count = len(genes)

            go_inf = getGONameAndCatergory(cols[0])

            if go_inf != "":
                go_cat,go_name = go_inf.split(",",1)
                if float(str(cols[1])) < pVal:
                    goDictContainer.append({"count": count,"genes":str(cols[2].strip("\r\"")), "GO_id": str(cols[0]), "GO_name": go_name, "cat": go_cat,
                        "pVal": str(cols[1])})

    return goDictContainer
def byteify(input, encoding='utf-8'):
    if isinstance(input, dict):
        return {byteify(key): byteify(value) for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode(encoding)
    else:
        return input

def parseInputMonagos(content_dict,pVal):
    goDictContainer = []
    num = len(content_dict["go_inf"])
    for index_monago in range(num):
        #content_dict["go_inf"][index_monago]["genes"] = ';'.join(map(str, content_dict["go_inf"][index_monago]["genes"]))
        pVal_each = float(str(content_dict["go_inf"][index_monago]["pVal"]))
        count = str(content_dict["go_inf"][index_monago]["count"])
        genes = str(content_dict["go_inf"][index_monago]["genes"])
        GO_id = str(content_dict["go_inf"][index_monago]["GO_id"])
        cat = str(content_dict["go_inf"][index_monago]["cat"])
        GO_name = str(content_dict["go_inf"][index_monago]["GO_name"])
        if pVal_each < pVal:
            goDictContainer.append({"count":count , "genes": genes, "GO_id": GO_id, "GO_name": GO_name,
                 "cat": cat,"pVal": pVal_each})
    return goDictContainer

@logTime
def getDataFromDavid(inputIds,idType,annotCat,pVal):
    '''
    send https request to David and get GO information
    
    Args:
        inputIds:a list of gene ids
        idType:the type of gene ids, such as AFFYMETRIX_3PRIME_IVT_ID
        annotcat: type of annotation wanted, such as GO
        pVal: p-value

    Return:
        A list of GO terms
    '''
    davidScrawler = DavidDataScrawler()
    davidScrawler.setParams(inputIds,idType,annotCat,pVal)

    try:
        go = davidScrawler.run()
    except Exception as e:
        logger.error(str(e))
        return [],False
    else:
        return go,True


@logTime
def processedData(go, similarity, clustComp):
    '''
    generate necessary data for visualization

    Args:
        a list of GO terms

    Return:
        matrix:a matrix M where M(i,j) represents the number of intersected genes betweeen GO[i] and GO[j]
        go_index_reord:an array representing the position change of GO terms after hieracical clustering
        go_hier:a list of GO that are ancesters of enriched GO terms.
        go_inf_reord:an array of enriched GO terms
        clusterHierData:an array storing hierarcical data use to generated hierarcical tree
    '''



    # redefine any go terms that are 'alternative ids' as their main id
    alt_ids = json.load(open(root_dir +"js/alt.js"))
    for index, id in enumerate(go):
        if id['GO_id'] in alt_ids.keys():
            id['GO_id'] = str(alt_ids[id['GO_id']])


    #for pre processing the data
    if similarity=='Resnik':
        from DataProcessResnikAve import DataProcess4
    elif similarity=='Genes':
        from DataProcessGenes import DataProcess4
    elif similarity=='Simrel':
        from DataProcessSimrelAve import DataProcess4


    dataProcess = DataProcess4()
    try:
        preProcessedData = dataProcess.dataProcess(go, clustComp)

    except KeyError:
        raise SemanticsError
    except Exception as e:
        logger.error(str(e))
    else:
        matrix = preProcessedData["matrix"]["matrix_count"]
        go_index_reord = preProcessedData["go_index_reord"]
        go_hier = preProcessedData["go_hier"]
        go_inf_reord = preProcessedData["go_inf"]
        clusterHierData = preProcessedData["clusterHierData"]
        simDict=preProcessedData['simDict']
        return matrix,go_index_reord,go_hier,go_inf_reord,clusterHierData, simDict


def loadConfig():
    pass

if __name__ == '__main__':
    loadConfig()
    loadGOHier()
    app.run(debug= (config["debug"]=="true"), host="0.0.0.0", port = 80, threaded = True)
