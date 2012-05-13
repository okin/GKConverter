#-*- coding: utf-8 -*-
from math import cos, pi, sqrt, tan, sin, atan

floor = round


class GKConverter:
	@classmethod
	def is_valid_gk_coordinate(cls, x, y):
		return not (x > 1000000 and y > 1000000)

	@classmethod
	def convert_GK_to_lat_long(cls, right, height, use_wgs84=None):
		(x, y) = cls.convertGaussKruegerToLatitudeLongitude(right, height)

		return cls.sevenParameterHelmertTransformation(x, y, use_wgs84)

	@classmethod
	def convertGaussKruegerToLatitudeLongitude(cls, right, height):
		x = right
		y = height
		
		#Check for invalid Parameters
		if (not ((x > 1000000) and (y > 1000000))):
			raise ValueError("No valid Gauss-Kruger-Code.")
		
		#Variables to prepare the geovalues
		GKRight = x
		GKHeight = y
		e2 = 0.0067192188
		c = 6398786.849
		rho = 180 / pi
		bII = (GKHeight / 10000855.7646) * (GKHeight / 10000855.7646)
		print 'bII: %s' % (bII, )

		bf = 325632.08677 * (GKHeight / 10000855.7646) * ((((((0.00000562025 * bII + 0.00022976983) * bII - 0.00113566119) * bII + 0.00424914906) * bII - 0.00831729565) * bII + 1))
		print 'bf: %s' % (bf, )
		
		bf /= 3600 * rho
		print 'bf: %s' % (bf, )
		co = cos(bf)
		print 'co: %s' % (co, )
		g2 = e2 * (co * co)
		print 'g2: %s' % (g2, )
		g1 = c / sqrt(1 + g2)
		print 'g1: %s' % (g1, )
		t = tan(bf)
		print 't: %s' % (t, )
		fa = (GKRight - floor(GKRight / 1000000) * 1000000 - 500000) / g1
		print 'fa: %s' % (fa, )

		GeoDezRight = ((bf - fa * fa * t * (1 + g2) / 2 + fa * fa * fa * fa * t * (5 + 3 * t * t + 6 * g2 - 6 * g2 * t * t) / 24) * rho)
		dl = fa - fa * fa * fa * (1 + 2 * t * t + g2) / 6 + fa * fa * fa * fa * fa * (1 + 28 * t * t + 24 * t * t * t * t) / 120
		print 'dl: %s' % (dl, )
		GeoDezHeight = dl * rho / co + floor(GKRight / 1000000) * 3
		
		return (GeoDezRight, GeoDezHeight)

	@classmethod
	def sevenParameterHelmertTransformation(cls, right, height, use_wgs84=False):
		#Variables used in the transformation
		earthRadius = 6378137 #Earth is a sphere witht this radius
		aBessel = 6377397.155
		eeBessel = 0.0066743722296294277832
		ScaleFactor = 0.00000982
		RotXRad = -7.16069806998785E-06
		RotYRad = 3.56822869296619E-07
		RotZRad = 7.06858347057704E-06
		ShiftXMeters = 591.28
		ShiftYMeters = 81.35
		ShiftZMeters = 396.39
		LatitudeIt = 99999999
		
		if(use_wgs84):
			ee = 0.0066943799
		else:
			ee = 0.00669438002290
		
		GeoDezRight = (right / 180) * pi
		GeoDezHeight = (height / 180) * pi
		
		n = eeBessel * sin(GeoDezRight) * sin(GeoDezRight)
		n = 1 - n
		n = sqrt(n)
		n = aBessel / n
		
		CartesianXMeters = n * cos(GeoDezRight) * cos(GeoDezHeight)
		CartesianYMeters = n * cos(GeoDezRight) * sin(GeoDezHeight)
		CartesianZMeters = n * (1 - eeBessel) * sin(GeoDezRight)
		
		CartOutputXMeters = (1 + ScaleFactor) * CartesianXMeters + RotZRad * CartesianYMeters - RotYRad * CartesianZMeters + ShiftXMeters
		CartOutputYMeters = -1 * RotZRad * CartesianXMeters + (1 + ScaleFactor) * CartesianYMeters + RotXRad * CartesianZMeters + ShiftYMeters
		CartOutputZMeters = RotYRad * CartesianXMeters - RotXRad * CartesianYMeters + (1 + ScaleFactor) * CartesianZMeters + ShiftZMeters
		
		GeoDezHeight = atan(CartOutputYMeters / CartOutputXMeters)
		
		Latitude = (CartOutputXMeters * CartOutputXMeters) + (CartOutputYMeters * CartOutputYMeters)
		Latitude = sqrt(Latitude)
		Latitude = CartOutputZMeters / Latitude
		Latitude = atan(Latitude)
		
		not_accurate_enough = True

		while(not_accurate_enough):
			LatitudeIt = Latitude

			n = 1 - ee * sin(Latitude) * sin(Latitude)
			n = sqrt(n)
			n = earthRadius / n

			Latitude = CartOutputXMeters * CartOutputXMeters + CartOutputYMeters * CartOutputYMeters
			Latitude = sqrt(Latitude)
			Latitude = (CartOutputZMeters + ee * n * sin(LatitudeIt)) / Latitude
			Latitude = atan(Latitude)

			not_accurate_enough = (abs(Latitude - LatitudeIt) >= 0.000000000000001)
			
		GeoDezRight = (Latitude / pi) * 180;
		GeoDezHeight = (GeoDezHeight) / pi * 180;

		return (GeoDezRight, GeoDezHeight)
	"""
	#Convert the Value into classic notation
	sub convertToClassicNotation
	{
	  my $param = pop(@_);
	  my $SecondsAbs = ($param - floor($param - floor(($param - floor($param)) * 60) / 60) * 60 * 60);
	  my $ShortedSecondsRest = $SecondsAbs - floor($SecondsAbs);
	  $SecondsAbs = floor($SecondsAbs);
	  $ShortedSecondsRest = floor($ShortedSecondsRest * 100); # -> Zwei Stellen

	  #Building of the Output-String
	  return(floor($param) . '°' . floor(($param - floor($param)) * 60) . "\'" . $SecondsAbs . '.' . $ShortedSecondsRest . '"');
	}"""

