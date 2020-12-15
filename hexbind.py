from skyfield.api import load
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from skyfield.api import utc
from scipy.optimize import brentq # machine learning

from datetime import timedelta, datetime
import pytz

# Custom helper functions
from definitions import *
from whereIsData import *


# convert binary to elemental representation
b_to_q = lambda x: [from_binary_to_element_symbols[x[i*2:i*2+2]] for i in [0,1,2]]
b_to_q_el = lambda x: [from_binary_to_element_ix[x[i*2:i*2+2]] for i in [0,1,2]]
b_to_ching = lambda x: [b for b in x]

# index elemental composition of the I Ching - binary becomes a 'triplet' with one of the four elements (air - fire - water - earth)
iching_ix = [b_to_q_el(str(x['binary']))for x in iching]

# binary position of the I Ching - binary becomes a string
iching_binary = [b_to_q(str(x['binary']))for x in iching]

# binary position of the I Ching - binary becomes an array of [1 and 0]
iching_binary_full = [b_to_ching(str(x['binary']))for x in iching]

def get_color_map():
	water = "#0E61B0"
	air = "#C29F17"
	earth = "#55A349"
	fire = "#C9280C"
	return [air, fire, water, earth]

def __test_hex_binary_to_element():
	print ([x for x in '10000'])

	print( from_binary_to_element_symbols['00'] )

	print ( b_to_q('000000'))
	print ( b_to_q_el('101010'))
	print ( b_to_ching('111000'))

def neutron_stream_pos(planet_position):
	""" returns mandala position (base 64) given planet position"""
	return ( (planet_position  + (2*line_width - 1*color_width - 2*tone_width) ) / (2*math.pi) * 64) % 64

def map_on_hexagram(df, include=[]):
	""" maps df planet positions onto position onto a hexagram and line """
	# convert dataframe to numpy array
	neutron_stream = df.to_numpy()

	hexagram_bin = np.floor(neutron_stream) # rounded up downwards

	# map bin number onto 'hexagram' (neutron stream is sequential order, hexagram is King Wen Sequence)
	strong = np.array(iching_map)
	flat = hexagram_bin.astype(int).flatten()

	previous_shape = neutron_stream.shape
	mapped = strong[flat]

	hexagram = mapped.reshape(previous_shape)
	hexagram_fraction = neutron_stream - hexagram_bin

	line = hexagram_fraction // (1/6) + 1 # count in which 6th this neutrino stream falls in
	line_fraction = (hexagram_fraction - (line - 1)*1/6 ) / (1/6)

	color = line_fraction // (1/6) + 1
	color_fraction = (line_fraction - (color -1) * 1/6) / (1/6)

	tone = color_fraction // (1/6) + 1

	
	return_info = [hexagram]
	
	if 'lines' in include:
		return_info += [line.astype(int)]
	
	return return_info

def get_elemental_ching_lines_map(df_planets):
	# input = ephemeris planets for a certain time period
	# map 'neutron' stream, aka influences of the probability field (the planets in the solar system physical space)
	# position to index on a wheel
	
	df_planets['earth'] = df_planets['sun'] - math.pi
	df_angles = neutron_stream_pos(df_planets)

	# index on a wheel to specific binary I-Ching Sequence - King Wen Version
	z = map_on_hexagram(df_angles, include=['lines'])

	hexagrams = z[0]
	lines = z[1]
	
	# many_2_b = np.array(iching_binary) # strong
	many_2 = np.array(iching_binary_full) # strong
	one_2  = hexagrams.astype(int).flatten() - 1 # flat

	# binary el
	#el_b = many_2_b[one_2]

	# normal el (0 -> 3)
	el = many_2[one_2]
	
	finish = el.reshape((df_angles.shape[0], df_angles.shape[1]*6))

	return finish.astype(int), lines

def get_elemental_ching_map(df_planets):
	# input = ephemeris planets for a certain time period
	# map 'neutron' stream, aka influences of the probability field (the planets in the solar system physical space)
	# position to index on a wheel
	
	df_planets['earth'] = df_planets['sun'] - math.pi
	df_angles = neutron_stream_pos(df_planets)

	# index on a wheel to specific binary I-Ching Sequence - King Wen Version
	z = map_on_hexagram(df_angles)
	
	# many_2_b = np.array(iching_binary) # strong
	many_2 = np.array(iching_binary_full) # strong
	one_2 = z.astype(int).flatten() - 1 # flat

	# binary el
	#el_b = many_2_b[one_2]

	# normal el (0 -> 3)
	el = many_2[one_2]
	
	finish = el.reshape((df_angles.shape[0], df_angles.shape[1]*6))

	return finish.astype(int)

def get_elemental_map(df_planets):
	# map 'neutron' stream, aka influences of the probability field (the planets in the solar system physical space)
	# position to index on a wheel
	
	df_planets['earth'] = df_planets['sun'] - math.pi
	df_angles = neutron_stream_pos(df_planets)

	# index on a wheel to specific binary I-Ching Sequence - King Wen Version
	z = map_on_hexagram(df_angles)
	
	many_2_b = np.array(iching_binary) # strong
	many_2 = np.array(iching_ix) # strong
	one_2 = z.astype(int).flatten() - 1 # flat

	# binary el
	el_b = many_2_b[one_2]

	# normal el (0 -> 3)
	el = many_2[one_2]
	
	finish = el.reshape((df_angles.shape[0], df_angles.shape[1]*3))

	return finish.astype(int)

def __test_neutron_stream_and_mapping(df_planets):
	# map 'neutron' stream, aka influences of the probability field (the planets in the solar system physical space)
	df_angles = neutron_stream_pos(df_planets.iloc[:, 1:6])
	z = map_on_hexagram(df_angles)
	print (z)

	many_2_b = np.array(iching_binary) # strong
	many_2 = np.array(iching_ix) # strong
	one_2 = z.astype(int).flatten() - 1 # flat

	print (many_2_b[63])

	# binary el
	el_b = many_2_b[one_2]

	# normal el (0 -> 3)
	el = many_2[one_2]
	print (el)

def get_crypto_planet_data(size, mapping_type='elemental'):
	""" Returns bitstamp data with size = # of ticks (in minutes for this dataset) 

	params:
		size: # of ticks (seconds in this case)
		mapping_type: 
			- 'elemental' (3-compound elements)
			- 'elemntal_ching' (I-Ching Plain Binary)
			- 'elemental_ching_lines' (I-Ching Plain Binary + Lines Dummified [1-6 lines -> 6 columns with 0 or 1])

	"""

	# get planetary positions 
	_planets, ts = get_planetary_ephemeris()

	color_map = get_color_map()

	df_c = pd.read_csv('bitstampUSD_1-min_data_2012-01-01_to_2020-09-14.csv', parse_dates=True)

	# make data timestamp
	df_c['date'] = pd.to_datetime(df_c['Timestamp'], unit='s')

	# cast down to hourly data
	groupkey = pd.to_datetime(df_c[-size:].date.dt.strftime('%Y-%m-%d %H'))
	df_hourly = df_c[-size:].groupby(groupkey).agg({'Close':'last','Volume_(BTC)':'sum'})
	df_hourly.head()

	first_date = df_hourly.iloc[0].name
	print ( first_date )

	# generate ephemerial elements
	h = first_date.hour
	hours_in_trading_code = len(df_hourly) # stock exchange count of number differences
	t_time_array = ts.utc(first_date.year, first_date.month, first_date.day, range(h,h+hours_in_trading_code), 0) # -3000 BC to 3000 BC, increments in hours


	# generate empheremis for time period
	df_crypto_planets = generate_planets(_planets, ts, t_time_array) # can take a while

	# selected desired planets for attribution
	r = ['earth','moon','mercury','venus','sun', 'mars', 'jupiter','saturn', 'uranus','neptune']
	r.reverse()

	# create elemental data map
	if mapping_type == 'elemental':
		data_tmp = get_elemental_map(df_crypto_planets.loc[:,r])
	elif mapping_type == 'elemental_ching':
		data_tmp = get_elemental_ching_map(df_crypto_planets.loc[:,r])
	elif mapping_type == 'element_ching_lines':
		data_tmp, lines = get_elemental_ching_lines_map(df_crypto_planets.loc[:,r])
		
#     return data_tmp, lines
		
	# plot data map
	fig, ax = plt.subplots(figsize=(5,5))
	sns.heatmap(data_tmp.transpose(), ax=ax, cmap=color_map, cbar=False)
	
	fig, ax = plt.subplots(figsize=(5,5))
	sns.heatmap(lines.transpose(), ax=ax, cmap=color_map, cbar=False)
	
	if mapping_type == 'elemental' or mapping_type == 'elemental_ching':
		# create the training dataset [Close, Solar System Time]
		df_solar = pd.DataFrame(data_tmp)
		df_solar.index = df_hourly.index
		
		df_dataset = pd.concat([df_hourly[['Close']], df_solar], axis=1)
		return df_dataset
	

	elif mapping_type == 'element_ching_lines':
		
		# create the training dataset [Close, Solar System Time]
		df_solar = pd.DataFrame(data_tmp)
		df_solar.index = df_hourly.index
		
		df_lines = pd.DataFrame(lines, columns = [str(x) for x in range(10)] ).astype(str)
		df_lines = pd.get_dummies(df_lines)
		df_lines.index = df_hourly.index
		
		df_dataset = pd.concat([df_hourly[['Close']], df_solar, df_lines], axis=1)
		return df_dataset
		



