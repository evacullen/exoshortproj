import numpy as np
from astropy.io import fits

# Load AXY file from Astrometry.net (detected stars)
axy_file = "exoshortproj/fits_files/axy.fits"
axy_data = fits.open(axy_file)[1].data  # Star data

# Extract X & Y positions
x_stars = axy_data['X']
y_stars = axy_data['Y']

# Print first 10 detected stars
for i in range(10):
    print(f"Star {i+1}: X = {x_stars[i]:.2f}, Y = {y_stars[i]:.2f}")
