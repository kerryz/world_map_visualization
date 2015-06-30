# Description
World map visualization using various map projections. A map projection is a systematic transformation of the latitudes and longitudes of locations on the surface of a sphere or an ellipsoid into locations on a plane. The map projections used are:

* Azimuthal equidistant projection
* Sinusoidal (also known as Mercator equal-area or Sanson-Flamsteed) projection
* Mercator projection

Color codes will be used to visualizes the prosperity levels of different regions, cities' population sizes, and you can choose to draw the geodesic path (shortest path between two lines on a sphere) between Beijing and Los Angeles through a command line interface. The areas (in pixels) will be calculated for Japan and Taiwan to compare how the different projection methods distort the projected areas.

# Requirements
* Python 2.7
* numpy
* matplotlib

# How to run

Have the files `map_projections.py` and `projection_formulae.py` in the same directory. Browse to that directory in a terminal and write the command:

	python map_projections.py [data_directory_filepath]

where `data_directory_filepath` is an optional argument that is the filepath to the root data directory.  
Defaults to `../data/`. 

# Data description

The data directory should have the following structure:

	data
	    |-- regions
	        |-- 000
	            |-- landparts.txt
	            |-- landshape.txt
	        |-- ...
	    |-- cities.txt
	    |-- landparts.txt
	    |-- landshape.txt
	    |-- regions.txt

* `landshape.txt`
	* Each row contains the coordinates of a shoreline.
* `landparts.txt`
	* Describes which of the coordinates in `landshape.txt` belong to the same land part.
	* Example: the first lines in this file is: 0,40,63, meaning that the coordinates in row 0 until 39, and 40 until 63 form their own land parts.
* `regions.txt`
	* Describes prosperity levels of 254 geographical regions.
	* Each line is of the form: `region_name|region_prosperity_level (1 highest, 7 lowest)`
* `regions` directory
	* Contains the landshape and landparts description of each region.
* `cities.txt`
	* Contains the attributes of 7322 cities. Each row contains the attributes of one city, which are separated by "|". The attributes are (in order):
		* City name
		* Population
		* Longitude
		* Latitude
