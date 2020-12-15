from skyfield.api import load
from skyfield.api import utc

import pandas as pd
import numpy as np
import math

def get_planetary_ephemeris():
	""" load planetary 'skyobject-files' using Skyfield API """
	# load Ephemeris
	# load library planet position file
	planets = load('de431t.bsp')

	_planets = {
		"sun": planets['SUN'],
		"earth": planets['EARTH'],
		"moon": planets['MOON'],
		"mercury": planets['MERCURY'],
		"venus": planets['VENUS BARYCENTER'],
		"mars": planets['MARS BARYCENTER'],
		"jupiter": planets['JUPITER BARYCENTER'],
		"saturn": planets['SATURN BARYCENTER'],
		"uranus": planets['URANUS BARYCENTER'],
		"neptune": planets['NEPTUNE BARYCENTER'],
		"pluto": planets['PLUTO BARYCENTER'],
	}

	ts = load.timescale(builtin=True)

	return _planets, ts

def generate_planets(_planets, ts, time_range):
	""" Generates a dataframe for the registered collection of dates
	
	Returns
		Dataframe with ephemeris of planets at location of 

	Usage
		from skyfield.api import load
		from skyfield.api import utc

		_planets, ts = get_planetary_ephemeris()
		hours_to_calculate = 24 * 365 * 10
		t_time_array = ts.utc(2010, 1, 1, range(0,hours_to_calculate), 0)
		df_planets = generate_planets(_planets, ts, t_time_array) # takes a while

	"""
	df = pd.DataFrame()
	df['time'] = np.array(time_range.tt)

	for planet_name, p in _planets.items():

		if planet_name == "earth":
			earth_offset = math.pi
			p = _planets['sun']
		else:
			earth_offset = 0

		df[planet_name] = _planets['earth'].at(time_range).observe(p).ecliptic_latlon()[1].radians + earth_offset # latitude
		
	return df

