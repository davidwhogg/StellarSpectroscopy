StellarSpectroscopy
===================

Quick summary:
-------------
**sqlcl.py** should be downloaded if it's not on the computer already. Otherwise, run **graphcode.py** (which will import **sqldata_notmain.py**) and it will produce the graph of two particular variables & a particular spectrum at a particular data point.

Explanation of code:
-------------------
* **filterstars.py** : This filters for the calibration star using a SQL query embedded within the python script.
* **graphcode.py** : This code graphs the data from sqldata_notmain and also generates a spectrum for a particular point.
* **redshiftcorr_separatepics.py** : This code corrects for the redshift for some spectra, using data from a fits file. It should be absorbed into **graphcode.py** eventually.
* **sqlcl.py** : This is a python module necessary for running SQL queries in Python.
* **sqldata_notmain**: This is imported as a module into **graphcode.py** and contains the query from **filterstars.py**, and essentially prepares the data into arrays.

*Note: I have made the codes that require the sqlcl module import it within the script, so that if you want to run the scripts, simply download sqlcl.py from here and put it in the same folder as the code to run (e.g. sqldata_notmain.py)

