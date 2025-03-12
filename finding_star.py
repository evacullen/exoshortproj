from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

# Load the FITS file
fits_file = "exoshortproj/fits_files/wcs.fits"  # Use the same file uploaded to astrometry.net
hdul = fits.open(fits_file)
header = hdul[0].header
wcs = WCS(header)

# TOI 46432937 b RA/Dec (from ExoFOP)
ra_toi = "05:35:28.56"
dec_toi = "-14:35:49.89"
coord_toi = SkyCoord(ra_toi, dec_toi, unit=(u.hourangle, u.deg))

# Extract RA & Dec in degrees
ra_deg = coord_toi.ra.deg
dec_deg = coord_toi.dec.deg

# Convert RA/Dec to pixel coordinates
x_toi, y_toi = wcs.world_to_pixel_values(ra_deg, dec_deg)  # ✅ FIXED

print(f"TOI 46432937 b Pixel Coordinates: X = {x_toi:.2f}, Y = {y_toi:.2f}")
