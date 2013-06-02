def __main__():

#    import sys
#    directory = commands.getoutput("pwd")
#    sys.path.append(directory)
    import sqlcl

    alldata=  sqlcl.query("SELECT p.objID, \
    p.extinction_g, s.elodieTEff, s.elodieFeH, s.elodieObject, p.camcol, p.run, p.field, \
    p.obj, s.plate, s.fiberID, s.mjd, p.ra, p.dec, \
    p.psfMag_u, p.psfMag_g, p.psfMag_r, p.psfMag_i, p.psfMag_z \
    FROM PhotoObj AS p \
    JOIN SpecObj as s ON s.specobjID=p.specobjID \
    WHERE psfMag_r<19 \
    and psfMag_r>15 and type=6 \
    and dbo.fPhotoStatus('PRIMARY')>0 and dbo.fPhotoFlags('STATIONARY')>0 \
    and calibStatus_r=1 \
    and s.elodieTEff!=0 and s.elodieFeH!=0 and s.elodieLogG!=0 \
    and ((flags&dbo.fPhotoFlags('BLENDED')) \
    +(flags&dbo.fPhotoFlags('DEBLEND_TOO_MANY_PEAKS')) + \
    (flags&dbo.fPhotoFlags('SATURATED')) \
    +(flags&dbo.fPhotoFlags('BADSKY'))+ \
    (flags&dbo.fPhotoFlags('COSMIC_RAY')) \
    +(flags&dbo.fPhotoFlags('PEAKS_TOO_CLOSE'))+ \
    (flags&dbo.fPhotoFlags('NOTCHECKED_CENTER')) \
    +(flags&dbo.fPhotoFlags('SATUR_CENTER'))+ \
    (flags&dbo.fPhotoFlags('INTERP_CENTER')) \
    +(flags&dbo.fPhotoFlags('INTERP'))+ \
    (flags&dbo.fPhotoFlags('PSF_FLUX_INTERP')))=0 \
    and sqrt((power(psfMag_u-psfmag_g-0.82,2)+ \
    power(psfMag_g-psfMag_r-0.3,2)+ \
    power(psfMag_r-psfMag_i-0.09,2)+\
    power(psfMag_i-psfMag_z-0.02,2)))<0.20").read()

    # filters for calibration stars taken from
    # http://www.sdss3.org/dr9/algorithms/boss_std_ts.php

__main__()
