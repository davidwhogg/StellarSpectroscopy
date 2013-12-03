StellarSpectroscopy
===================

Data Compilation

* **eqw-newer.py**: this module calculates the continuum, flux and equivalent width for a dictionary of spectral line locations and outputs the data as a pickled file for use in **plotting.py** and **fit_coeffs_new.py** Sample picture: http://postimg.org/image/nniirvr7r/.
* 

Fitting
* **analysis.py** has a bunch of functions which allow for the individual inspection of the i-th spectrum. Currently it also plots a 3x3 figure of instances of failures of the above estimation module, e.g. http://postimg.org/image/6o10qb157/.
* **fit_coeffs_new** fits for the coefficients for an attentuation law involving brightness, HD-EW and extinction_g using a standardized wavelength grid, and generates a covariance matrix (currently investigating accuracy of this) and outputs the files for use in **plotting.py**.

Plotting
* **plotting.py** plots the coefficients obtained from **fit_coeffs_new.py**, with the possibility of labelling various spectral lines. Sample plot: http://postimg.org/image/gmulxqnov/
* * **histotriangle.py** which will graph histograms of the variables, using the stars from **filterstars.py** as well as a contour plot of two variables. Sample picture: http://postimg.org/image/s3pdr2dhx/

Currently not used:
* * **LDA with SQLdata.py** which will take the same stars from **filterstars.py** and apply Linear Discriminant Analysis to them. Right now, the two classes have been arbitrarily assigned and the resulting analysis is meaningless. Sample picture: http://postimg.org/image/7d3d02as5/
* **plotmanyspectra.py** will sort and plot the ratio of EW's to the median, broken up into 9 quantiles. The effectiveness of the current algorithm is disputed.
* **classes.py** tries to sort the data by extinction and HD-EW, currently facing same problem as **plotmanyspectra.py**

Notes:
* sqlcl.py, found in this repository, is required to run the SQL queries on SDSS in **eqw-newer.py**

Explanation of code:
-------------------
* **filterstars.py** : This filters for the calibration star using a SQL query embedded within the python script.
* **graphcode.py** : This code graphs the data from sqldata_notmain and also generates a spectrum for a particular point.
* **redshiftcorr_separatepics.py** : This code corrects for the redshift for some spectra, using data from a fits file. It should be absorbed into **graphcode.py** eventually.
* **sqlcl.py** : This is a python module necessary for running SQL queries in Python.
* **sqldata_notmain**: This is imported as a module into **graphcode.py** and contains the query from **filterstars.py**, and essentially prepares the data into arrays.


