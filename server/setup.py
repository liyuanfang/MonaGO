import pip

def install(package):
    pip.main(['install', package])

if __name__ == '__main__':
    install('flask')
    install('numpy')
    install('requests')
    install('networkx')
    install('geotools')
    install('biopython')
    install('goatools')
    install('DavidDataScrawler')
