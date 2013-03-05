import sys
sys.path.append('/Users/admin/Desktop/Maser files')
import sqlcl

print sqlcl.query("SELECT top 10 objID \
FROM PhotoObj WHERE psfMag_r<19 \
and psfMag_r>15 and type=6 \
and  dbo.fPhotoStatus('PRIMARY')>0\
and dbo.fPhotoFlags('STATIONARY')>0\
and calibStatus_r=1\
and ((flags&dbo.fPhotoFlags('BLENDED'))+(flags&dbo.fPhotoFlags('DEBLEND_TOO_MANY_PEAKS'))+\
(flags&dbo.fPhotoFlags('SATURATED'))+(flags&dbo.fPhotoFlags('BADSKY'))+\
(flags&dbo.fPhotoFlags('COSMIC_RAY'))+(flags&dbo.fPhotoFlags('PEAKS_TOO_CLOSE'))+\
(flags&dbo.fPhotoFlags('NOTCHECKED_CENTER'))+(flags&dbo.fPhotoFlags('SATUR_CENTER'))+\
(flags&dbo.fPhotoFlags('INTERP_CENTER'))+(flags&dbo.fPhotoFlags('INTERP'))+\
(flags&dbo.fPhotoFlags('PSF_FLUX_INTERP')))=0\
and sqrt((power(psfMag_u-psfmag_g-0.82,2)+power(psfMag_g-psfMag_r-0.3,2)+\
power(psfMag_r-psfMag_i-0.09,2)+power(psfMag_i-psfMag_z-0.02,2)))<0.08").read()
