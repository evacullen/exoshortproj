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
    "planetary_parameters": {
        "Planet Name": "TIC 46432937 b",
        "Target Star RA": "05:35:28.56",
        "Target Star Dec": "-14:35:49.89",
        "Published Mid-Transit Time (BJD-UTC)": 10744.5338,
        "Orbital Period (days)": 1.44,
    },
    "user_info": {
        "Directory with FITS files": r"C:\Users\evacc\Downloads\exoshortproj\fits_files", #IMPORTANT: replace with the path to the fits_files folder
        "Directory to Save Plots": r"C:\Users\evacc\Downloads\exoshortproj\output", # replace with the path to the output folder
        "Target Star X & Y Pixel": [1026.80, 988.89], #updated.
        "Comparison Star(s) X & Y Pixel": [
            [1826.73, 1635.27],  #  Star 1
            [1576.36, 1506.55],  #  Star 2
            [1518.16, 1239.44]   #  Star 10
        ], #updated?
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
            scale_units='arcsecperpix',  # Add if you know the scale
            scale_type='ev',              # Estimate bounds
            scale_est=0.384,                # Example: 0.5"/pixel
            scale_err=20                   # 20% uncertainty
        )
        if wcs_header:
            return WCS(wcs_header)
    except Exception as e:
        print(f"Astrometry failed: {e}")
        return None
    
from astropy.coordinates import SkyCoord
from astropy import units as u

def verify_plate_solution(wcs, config):
    """Verify WCS solution matches expected target coordinates."""
    try:
        ra = config["planetary_parameters"]["Target Star RA"]
        dec = config["planetary_parameters"]["Target Star Dec"]
        
        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        calc_x, calc_y = wcs.all_world2pix(coord.ra.deg, coord.dec.deg, 0)
        
        expected_x, expected_y = config["user_info"]["Target Star X & Y Pixel"]
        offset = np.sqrt((calc_x - expected_x)**2 + (calc_y - expected_y)**2)
        
        if offset > 3:
            print(f"WARNING: Plate solution offset: {offset:.2f} pixels")
    except Exception as e:
        print(f"Plate solution verification failed: {str(e)}")

def plot_lightcurve(times, fluxes, config):
    """Plot and save the time-series light curve with binned averages."""
    plt.figure(figsize=(12, 6))
    
    # Convert times to hours from first observation
    time_offset = np.min(times)
    hours = (times - time_offset) * 24
    
    # Plot individual measurements
    plt.plot(hours, fluxes, 'bo', alpha=0.3, label='Individual measurements')
    
    # Add binned averages with error bars
    bin_width = 0.02  # hours - adjust this based on exposure time and cadence
    bins = np.arange(np.min(hours), np.max(hours) + bin_width, bin_width)
    
    bin_centers = []
    bin_means = []
    bin_errors = []
    
    for i in range(len(bins)-1):
        mask = (hours >= bins[i]) & (hours < bins[i+1])
        if np.sum(mask) > 0:
            bin_fluxes = fluxes[mask]
            bin_centers.append((bins[i] + bins[i+1])/2)
            bin_means.append(np.mean(bin_fluxes))
            bin_errors.append(np.std(bin_fluxes)/np.sqrt(len(bin_fluxes)))  # Standard error
            
    plt.errorbar(bin_centers, bin_means, yerr=bin_errors, 
                fmt='r-', ecolor='darkred', elinewidth=2, capsize=3,
                label=f'{bin_width} hr binned average')
    
    # Formatting
    plt.title(f"{config['planetary_parameters']['Planet Name']} Light Curve")
    plt.xlabel('Hours from First Observation')
    plt.ylabel('Normalized Flux')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Add expected transit duration shading
    transit_duration = 1.6  # Hours - adjust this based on transit duration
    plt.axvspan(-transit_duration/2, transit_duration/2, 
                color='gray', alpha=0.2, label='Expected transit window')
    
    # Save plot
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

        # ---> ADDED: Extract observation time <---
        try:
            jd_time = Time(header['DATE-OBS'], format='isot').jd
        except KeyError:
            print(f"Missing DATE-OBS keyword in {file}; skipping this file.")
            continue

        # Get WCS and update target position
        wcs = None
        current_target_pos = config["user_info"]["Target Star X & Y Pixel"]  # Default
        if config["user_info"]["Plate Solution? (y/n)"] == "y":
            wcs = solve_astrometry(file)
            if wcs and wcs.is_celestial:
                try:
                    ra = config["planetary_parameters"]["Target Star RA"]
                    dec = config["planetary_parameters"]["Target Star Dec"]
                    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
                    target_x, target_y = wcs.all_world2pix(coord.ra.deg, coord.dec.deg, 0)
                    current_target_pos = [target_x, target_y]
                    # Verify offset from user's expected position
                    expected_x, expected_y = config["user_info"]["Target Star X & Y Pixel"]
                    offset = np.hypot(target_x - expected_x, target_y - expected_y)
                    if offset > 10:  # Alert if large discrepancy
                        print(f"WARNING: {file} target offset {offset:.1f} pixels")
                except Exception as e:
                    print(f"WCS conversion failed: {e}")
                    current_target_pos = config["user_info"]["Target Star X & Y Pixel"]

        # perform photometry
        target_flux = measure_flux_pixel(data, current_target_pos)
        comp_fluxes = [measure_flux_pixel(data, pos) for pos in comp_positions]
        
        if target_flux and all(comp_fluxes):
            norm_flux = target_flux / np.sum(comp_fluxes)
            times.append(jd_time)  # <--- NOW jd_time IS DEFINED
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
