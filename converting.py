from astropy.coordinates import SkyCoord
import astropy.units as u

# Star #1, #2, #10 in degrees (from your list):
stars_deg = [
    (83.820658, -14.696552),  # Star 1
    (83.827917, -14.667334),  # Star 2
    (83.855008, -14.654530),  # Star 10
]

for i, (ra_deg, dec_deg) in enumerate(stars_deg, start=1):
    sc = SkyCoord(ra_deg*u.deg, dec_deg*u.deg)
    print(f"Star {i} in sexagesimal = {sc.to_string('hmsdms')}")
