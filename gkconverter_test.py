import unittest
from gkconverter import GKConverter

class GKConverterTest(unittest.TestCase):
	def testSimpleConversion(self):
		(x, y) = GKConverter.convert_GK_to_lat_long(3477733, 5553274)

		self.assertAlmostEqual(50.11526435691097, x)
		self.assertAlmostEqual(8.687625204011725, y)

	def testSimpleGKtoLatLongConversion(self):
		right = 3477733
		height = 5553274

		(x, y) = GKConverter.convertGaussKruegerToLatitudeLongitude(right, height)

		self.assertAlmostEqual(50.1164273930041, x)
		self.assertAlmostEqual(8.6886330000005, y)

	def testHelmertTransformation(self):
		(x, y) = GKConverter.sevenParameterHelmertTransformation(50.1164273930041, 8.6886330000005)

		self.assertAlmostEqual(50.11526435691097, x)
		self.assertAlmostEqual(8.687625204011725, y)

if __name__ == '__main__':
	unittest.main()