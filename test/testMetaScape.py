import unittest, sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'server'))
from MetaScapeScrawler import MetaScapeScrawler
# from DataProcess import DataProcess


class Test(unittest.TestCase):
	def setUp(self):
		self.URLWSGI = "http://metascape.org/gp_server"
		self.metaScapeScrawler = MetaScapeScrawler()
		# DataProcess = DataProcess()
		self.metaScapeScrawler.inputIds = ""

	def tearDown(self):
		pass

	def testGetSessionId(self):
		pass

	# def testGuessIdType(self):
	# 	self.metaScapeScrawler.inputIds = "Actc1,Actn2,Adk,Adnp,Adprhl1,Akap2,Ankrd1,Anp32e,Anxa6,App,Arpc2,Asb2,Atp1a1,Atp1b1,Atp2a2,Atp5o,Cald1,Cap2,Casq2,Cav3"
		
	# 	#with self.assertRaises(Exception):
	# 	result = self.metaScapeScrawler.guessIDType()
	# 	self.assertTrue(result == "Gene Symbol","guessIDType failed")

	def testSaveTextBoxToserver(self):
		self.metaScapeScrawler.sessionId = self.metaScapeScrawler.createNewSession()
		self.metaScapeScrawler.inputIds = "Actc1,Actn2,Adk,Adnp,Adprhl1,Akap2,Ankrd1,Anp32e,Anxa6,App,Arpc2,Asb2,Atp1a1,Atp1b1,Atp2a2,Atp5o,Cald1,Cap2,Casq2,Cav3"
		
		#with self.assertRaises(Exception):
		result = self.metaScapeScrawler.saveTextBoxToServer()
		self.assertTrue(result == "good","save textbox to server failed")

	# def testGuessSpeciesFromInput(self):
	# 	self.metaScapeScrawler.sessionId = self.metaScapeScrawler.createNewSession()
	# 	self.metaScapeScrawler.inputIds = "Actc1,Actn2,Adk,Adnp,Adprhl1,Akap2,Ankrd1,Anp32e,Anxa6,App,Arpc2,Asb2,Atp1a1,Atp1b1,Atp2a2,Atp5o,Cald1,Cap2,Casq2,Cav3"
	# 	result = self.metaScapeScrawler.guessSpeciesFromInput()
	# 	self.assertTrue("4896" in result ,"guess species from input failed")

	def testApplySpecies(self):
		self.testSaveTextBoxToserver()
		result = self.metaScapeScrawler.applySpecies(9606,-1)
		self.assertTrue(result == "{}" ,"apply species from input failed")

if __name__ == '__main__':
    unittest.main()