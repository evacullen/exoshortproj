import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus
from astroquery.astrometry_net import AstrometryNet
from astropy.time import Time
import os
from glob import glob

# UNCOMMENT BELOW AND PASTE IN TERMINAL. this code was developed in the Anaconda enviornment.
# pip install numpy
# pip install matplotlib
# pip install astropy
# pip install photutils
# pip install astroquery

# PLACE ALL FITS FILES DIRECTLY IN "fits_files" FOLDER. COPY THE PATH TO THIS FOLDER. 
# load your configuration
config = {
    "planetary_parameters": { # all parameters need to be updated for TOI. can someone get this info? 
        "Planet Name": "WASP-43 b", 
        "Target Star RA": "10:19:37.964856",
        "Target Star Dec": "-09:48:23.19516",
        "Published Mid-Transit Time (BJD-UTC)": 2457202.184881,
        "Orbital Period (days)": 0.81347406,
        },
    "user_info": {
        "Directory with FITS files": r"C:\Users\evacc\Downloads\exoshortproj", #IMPORTANT: replace with the path to the fits_files folder
        "Directory to Save Plots": r"C:\Users\evacc\Downloads\exoshortproj\output", # replace with the path to the output folder
        "Target Star X & Y Pixel": [1000, 408], #also needs to be updated for TOI.
        "Comparison Star(s) X & Y Pixel": [[530, 10], [900, 200]], #also also needs to be updated for TOI.
        "Plate Solution? (y/n)": "y"
    }
}

# photometry parameters
APERTURE_RADIUS = 5.0  # start with this, adjust based on FWHM, may need to be made smaller/larger depending on seeing conditions
ANNULUS_RADII = (8.0, 12.0)

def measure_flux_pixel(data, position):
    aperture = CircularAperture(position, r=APERTURE_RADIUS)
    annulus = CircularAnnulus(position, r_in=ANNULUS_RADII[0], r_out=ANNULUS_RADII[1])
    
    # background subtraction
    annulus_masks = annulus.to_mask(method='center')
    annulus_data = annulus_masks.multiply(data)
    mask = annulus_masks.data > 0
    bkg_median = np.median(annulus_data[mask])
    bkg = bkg_median * aperture.area
    
    # photometry
    phot_table = aperture_photometry(data, aperture)
    return phot_table['aperture_sum'][0] - bkg

def solve_astrometry(filename):
    ast = AstrometryNet()
    ast.api_key = 'uztwalfyyickpwxa'  # IMPORTANT replace with your API key from nova.astrometry.net (or use mine, it's fine)
    
    try:
        # try to see if there is already a valid WCS in the FITS header
        with fits.open(filename) as hdul:
            wcs = WCS(hdul[0].header)
            if wcs.is_celestial:
                return wcs
    except:
        pass

    # if no valid WCS, solve via Astrometry.net
    try:
        wcs_header = ast.solve_from_image(
            filename,
            publicly_visible="n",
            allow_commercial_use="n",
            solve_timeout=300,
            force_image_upload=True 
        )
        if wcs_header is not None:
            return WCS(wcs_header)
    except Exception as e:
        print(f"Astrometry failed for {filename}: {e}")
        return None
    
from astropy.coordinates import SkyCoord
from astropy import units as u

def verify_plate_solution(wcs, config):
    """Verify WCS solution matches expected target coordinates."""
    try:
        # convert target RA/Dec to pixels using WCS
        ra = config["planetary_parameters"]["Target Star RA"]
        dec = config["planetary_parameters"]["Target Star Dec"]
        
        # if RA/Dec are in sexagesimal format, you may need to parse them with SkyCoord (we need this for anything not in angular i.e. h:m:s...)
        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        calc_x, calc_y = wcs.all_world2pix(coord.ra.deg, coord.dec.deg, 0)
        
        calc_x, calc_y = wcs.all_world2pix(ra, dec, 0)
        
        # compare with user-provided coordinates
        expected_x, expected_y = config["user_info"]["Target Star X & Y Pixel"]
        offset = np.sqrt((calc_x - expected_x)**2 + (calc_y - expected_y)**2)
        
        # print warning if offset is above tolerance
        if offset > 3:
            print(f"WARNING: Plate solution offset: {offset:.2f} pixels")
    except Exception as e:
        print(f"Plate solution verification failed: {str(e)}")

def plot_lightcurve(times, fluxes, config):
    """Plot and save the time-series light curve."""
    plt.figure(figsize=(12, 6))
    
    # convert times to hours from first observation
    time_offset = np.min(times)
    hours = (times - time_offset) * 24
    
    plt.plot(hours, fluxes, 'bo', alpha=0.7)
    
    # formatting
    plt.title(f"{config['planetary_parameters']['Planet Name']} Light Curve")
    plt.xlabel('Hours from First Observation')
    plt.ylabel('Normalized Flux')
    plt.grid(True, alpha=0.3)
    
    # save plot
    output_path = os.path.join(config["user_info"]["Directory to Save Plots"], 
                               "lightcurve.png")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def plot_phased_lightcurve(times, fluxes, config):
    """Plot and save the phased light curve."""
    period = config["planetary_parameters"]["Orbital Period (days)"]
    t0 = config["planetary_parameters"]["Published Mid-Transit Time (BJD-UTC)"]
    
    phased_times = ((times - t0) / period) % 1
    phased_times[phased_times > 0.5] -= 1  # center transit at phase 0
    
    plt.figure(figsize=(10, 6))
    plt.plot(phased_times, fluxes, 'bo', alpha=0.7)
    plt.xlabel('Phase')
    plt.ylabel('Normalized Flux')
    plt.title('WASP-43 b Phased Light Curve')
    plt.grid(True)
    
    output_path = os.path.join(config["user_info"]["Directory to Save Plots"], 
                               "phased_lightcurve.png")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def process_images(config):
    """Main function that loops over FITS files, does photometry, solves astrometry if requested, 
    and produces light curves."""
    # create output directory if it doesn't exist
    os.makedirs(config["user_info"]["Directory to Save Plots"], exist_ok=True)

    # get FITS files
    fits_dir = config["user_info"]["Directory with FITS files"]
    files = glob(os.path.join(fits_dir, "*.fits"))
    
    times = []
    fluxes = []
    
    # photometry positions
    target_pos = config["user_info"]["Target Star X & Y Pixel"]
    comp_positions = config["user_info"]["Comparison Star(s) X & Y Pixel"]

    for file in files:
        with fits.open(file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            
        # solve astrometry & verify plate solution if required
        if config["user_info"]["Plate Solution? (y/n)"] == "y":
            wcs = solve_astrometry(file)
            if wcs:
                verify_plate_solution(wcs, config)

        # get observation time
        try:
            jd_time = Time(header['DATE-OBS'], format='isot').jd
        except KeyError:
            print(f"Missing DATE-OBS keyword in {file}; skipping this file.")
            continue

        # perform photometry
        target_flux = measure_flux_pixel(data, target_pos)
        comp_fluxes = [measure_flux_pixel(data, pos) for pos in comp_positions]
        
        if target_flux and all(comp_fluxes):
            norm_flux = target_flux / np.sum(comp_fluxes)
            times.append(jd_time)
            fluxes.append(norm_flux)

    # convert to arrays
    times = np.array(times)
    fluxes = np.array(fluxes)
    
    # generate and save plots
    plot_lightcurve(times, fluxes, config)
    plot_phased_lightcurve(times, fluxes, config)
    
    return times, fluxes

if __name__ == "__main__":
    times, fluxes = process_images(config)
    print("Times:", times)
    print("Fluxes:", fluxes)
