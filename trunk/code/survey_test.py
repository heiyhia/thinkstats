import unittest
import survey

class Test(unittest.TestCase):

    def testRespondents(self):
        resp = survey.Respondents()
        resp.ReadRecords()
        self.assertEquals(len(resp.records), 7643)

    def testPregnancies(self):
        preg = survey.Pregnancies()
        preg.ReadRecords()
        self.assertEquals(len(preg.records), 13593)

if __name__ == "__main__":
    unittest.main()
