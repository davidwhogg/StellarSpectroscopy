StellarSpectroscopy
===================

Quick summary:
-------------
* **sqlcl.py** should be downloaded if it's not on the computer already. 
* * Current work is in **eqw_estimate.py**: this function estimates the equivalent widths and continuum values of the h-b, h-gamma, h-d peaks for *n* stars with extinction values between *exta* and *extb*. It returns an array with the integrated fluxes, continuum values, equivalent widths and peak locations. Sample picture: http://postimg.org/image/nniirvr7r/. **NOTE: eqw_est_improved.py does the same thing as eqw_estimate.py but it chooses the region for flux integration by finding the nearest point of intersection of the spectrum with the continuum**

* All below programes will import **sqldata_notmain.py**, which is a copy of the query from **filterstars.py**.
Otherwise, run:
* * **graphcode.py** and it will produce the graph of two particular variables & a particular spectrum at a particular data point.
* * **histotriangle.py** which will graph histograms of the variables, using the stars from **filterstars.py** as well as a contour plot of two variables. Sample picture: http://postimg.org/image/s3pdr2dhx/
* * **LDA with SQLdata.py** which will take the same stars from **filterstars.py** and apply Linear Discriminant Analysis to them. Right now, the two classes have been arbitrarily assigned and the resulting analysis is meaningless. Sample picture: http://postimg.org/image/7d3d02as5/

Explanation of code:
-------------------
* **filterstars.py** : This filters for the calibration star using a SQL query embedded within the python script.
* **graphcode.py** : This code graphs the data from sqldata_notmain and also generates a spectrum for a particular point.
* **redshiftcorr_separatepics.py** : This code corrects for the redshift for some spectra, using data from a fits file. It should be absorbed into **graphcode.py** eventually.
* **sqlcl.py** : This is a python module necessary for running SQL queries in Python.
* **sqldata_notmain**: This is imported as a module into **graphcode.py** and contains the query from **filterstars.py**, and essentially prepares the data into arrays.

*Note: I have made the codes that require the sqlcl module import it within the script, so that if you want to run the scripts, simply download sqlcl.py from here and put it in the same folder as the code to run (e.g. sqldata_notmain.py)

