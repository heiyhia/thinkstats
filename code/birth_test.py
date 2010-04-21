import unittest
import birth

class Test(unittest.TestCase):

    def testRespondents(self):
        resp = birth.Respondents()
        self.assertEquals(len(resp.records), 7643)

    def testPregnancies(self):
        preg = birth.Pregnancies()
        self.assertEquals(len(preg.records), 13593)
        pass

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testCdf']
    unittest.main()
