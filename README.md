# Prerequisites
* Python 2.7
* numpy
* matplotlib

# How to run

Have the files `map_projections.py` and `projection_formulae.py` in the same directory. Browse to that directory in a terminal and write the command:

	python map_projections.py [data_directory_filepath]

where `data_directory_filepath` is an optional argument that is the filepath to the root data directory.  
Defaults to `../data/`. 

This directory should have the following structure:

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
