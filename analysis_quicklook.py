"""
ALMA data analysis (UK ALMA Workshop 2025)
Author: Daniel Walker

Run this code inside a CASA session.
Please ensure that you have the following packaged installed inside your CASA environment:
    spectral-cube
    pyspeckit
    astropy
    regions

To install, simply run "pip install spectral-cube pyspeckit astropy regions" inside your CASA session.
    
Basic analysis of ALMA data

In this short hands-on session we will look at some very basic analysis tools to:

*   Measure image statistics
*   Extract, plot, and fit spectra from cubes
*   Make and plot position-velocity plots
*   Make and plot moment maps

"""

import os
import shutil
import pyspeckit
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
from regions import CirclePixelRegion, PixCoord

# Create output directory for plots if it doesn't exist
output_dir = 'plots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Update the filename and mask_file variables to match your data
filename = 'PN_Hb_5.spw_0.image.fits'
mask_file = 'region.crtf'

"""
Plot a single channel in the cube and calculate some statistics of this channel and of the full cube
Feel free to change this channel number and to explore different stats that you get from 
imstat: https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.information.imstat.html
"""

channel = 522
with fits.open(filename) as hdul:
    data = hdul[0].data[channel, :, :]

plt.figure()
plt.imshow(data, origin='lower', cmap='inferno')
plt.colorbar(label='Intensity (Jy/beam)')
plt.title(f'Channel {channel} of {os.path.basename(filename)}')
plt.savefig(os.path.join(output_dir, 'single_channel.png'), bbox_inches='tight', dpi=300)
plt.close()

print('\n')

stats_chan = imstat(filename, chans=str(channel))
min_val_chan = stats_chan['min'][0]
max_val_chan = stats_chan['max'][0]
rms_val_chan = stats_chan['rms'][0]
print(f"Channel {channel}: Min: {min_val_chan * 1e3:.2f} mJy/beam, Max: {max_val_chan * 1e3:.2f} mJy/beam, RMS: {rms_val_chan * 1e3:.2f} mJy/beam")

stats_full_cube = imstat(filename)
min_val_cube = stats_full_cube['min'][0]
max_val_cube = stats_full_cube['max'][0]
rms_val_cube = stats_full_cube['rms'][0]
print(f"Full Cube: Min: {min_val_cube * 1e3:.2f} mJy/beam, Max: {max_val_cube * 1e3:.2f} mJy/beam, RMS: {rms_val_cube * 1e3:.2f} mJy/beam")

"""
Extract a mean spectrum with [spectral-cube](https://spectral-cube.readthedocs.io/en/latest/), 
then plot it with [pyspeckit](https://pyspeckit.readthedocs.io/en/latest/index.html)
"""

cube = SpectralCube.read(filename)
cube.allow_huge_operations=True
mean_spectrum = cube.mean(axis=(1, 2))

sp = pyspeckit.Spectrum(data=mean_spectrum.value, xarr=mean_spectrum.spectral_axis.value)
fig = plt.figure(figsize=(8, 7))
sp.plotter(figure=fig)
sp.plotter.axis.set_xlabel('Frequency (Hz)')
sp.plotter.axis.set_ylabel('Intensity (Jy/beam)')
plt.savefig(os.path.join(output_dir, 'mean_spectrum.png'), bbox_inches='tight', dpi=300)
plt.close()

"""
Repeat for a circular region (roughly) centred on the source
"""

ny, nx = cube.shape[1:]
centre = PixCoord(x=nx/2, y=ny/2)
region = CirclePixelRegion(center=centre, radius=10)
channel_55 = cube[54].value

subcube = cube.subcube_from_regions([region])
mean_spectrum = subcube.mean(axis=(1, 2))

sp = pyspeckit.Spectrum(data=mean_spectrum.value, xarr=mean_spectrum.spectral_axis.value)
fig = plt.figure(figsize=(8, 7))
sp.plotter(figure=fig)
sp.plotter.axis.set_xlabel('Frequency (Hz)')
sp.plotter.axis.set_ylabel('Intensity (Jy/beam)')
plt.savefig(os.path.join(output_dir, 'circular_region_spectrum.png'), bbox_inches='tight', dpi=300)
plt.close()

"""
Fit multiple Gaussian components to the spectrum
* The guesses are initialised with amplitude, centroid, and FWHM
* Guesses can be quite crude and the fitter will usually do a good job
* Guessed don't have to be hard-coded ([see example](https://pyspeckit.readthedocs.io/en/latest/example_fromscratch.html) using moments to initialise guesses)

**The fit shown in this example is not necesarrily a good fit. This is just illustrative.**
"""

guesses = [0.005, 2.2635e11, 5e9,
           0.05, 2.267e11, 5e9,
           0.1, 2.269e11, 5e9]

sp = pyspeckit.Spectrum(data=mean_spectrum.value, xarr=mean_spectrum.spectral_axis.value)
fig = plt.figure(figsize=(8, 7))
sp.plotter(figure=fig)
sp.plotter.axis.set_xlabel('Frequency (Hz)')
sp.plotter.axis.set_ylabel('Intensity (Jy/beam)')
sp.specfit(fittype='gaussian', guesses=guesses)
plt.savefig(os.path.join(output_dir, 'gaussian_fit.png'), bbox_inches='tight', dpi=300)
plt.close()

"""
Create a position-velocity plot across the source
* This example is using CASA's [impv](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.impv.html) task
     * Also see [pvextractor](https://pvextractor.readthedocs.io/en/latest/#) for a standalone Python package
* The slice is defined by the start and end pixel coordinates
     * Feel free to play with these, and the position angle, to see how the PV plot changes
"""

channel = 522
start = [148,122]
end = [175,175]
pv_output = filename.replace('.fits', '.pv')

# Create the PV diagram using CASA
impv(imagename=filename,
     outfile=pv_output,
     mode='coords',
     start=start,
     end=end,
     overwrite=True)

pv_fits = pv_output + '.fits'
exportfits(imagename=pv_output,
           fitsimage=pv_fits,
           overwrite=True)

# Fix the WCS header before creating the plot
# Sometimes it complains about the time keywords and demands them to be lowercase, for some reason ...
def fix_wcs_header(header):
    time_keys = ['TIMESYS']
    for key in time_keys:
        if key in header and isinstance(header[key], str):
            header[key] = header[key].lower()
    return header

# Create the plot with the fixed header
fig = plt.figure(figsize=(20, 8))

with fits.open(filename) as hdul:
    channel_data = hdul[0].data[channel, :, :]
    header = hdul[0].header

ax1 = fig.add_subplot(1, 2, 1)
im1 = ax1.imshow(channel_data, origin='lower', cmap='inferno', aspect='auto')
cbar1 = plt.colorbar(im1, ax=ax1)
ax1.set_title('Channel ' + str(channel))
ax1.plot([start[0], end[0]], [start[1], end[1]], color='white', linestyle='--')

with fits.open(pv_fits) as hdul:
    pv_data = hdul[0].data
    pv_header = hdul[0].header
    pv_header = fix_wcs_header(pv_header)
    wcs_pv = WCS(pv_header)

ax2 = fig.add_subplot(1, 2, 2, projection=wcs_pv)
im2 = ax2.imshow(pv_data, origin='lower', cmap='inferno', aspect='auto')
cbar2 = plt.colorbar(im2, ax=ax2, label='Intensity (Jy/beam)')
ax2.set_title('Position-Velocity Diagram')
ax2.coords[0].set_axislabel('Position (arcsec)')
ax2.coords[1].set_axislabel('Velocity (km/s)')

plt.savefig(os.path.join(output_dir, 'position_velocity.png'), bbox_inches='tight', dpi=300)
plt.close()

"""
Create and plot moment maps
* Use CASA's [immoments](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.immoments.html) task to create moment maps
* In this example we have moment 0 (integrated intensity), moment 1 (intensity weighted coordinate / velocity field), and moment 8 (peak intensity).
     * Feel free to try different moments. You can find the explanations at the above link for the task documentation.
* Notice that the chans, includepix, and region parameters are commented out
     * I'd encourage you to re-run the below blocks, each time un-commenting one of these lines. How does this change the resulting plots, and why?
* Note that you can also use other tools to create moment maps, such as [spectral-cube](https://spectral-cube.readthedocs.io/en/latest/)
"""

moment_files = [filename.replace('.fits', '.moment.integrated'),
                filename.replace('.fits', '.moment.maximum'),
                filename.replace('.fits', '.moment.weighted_coord')]

for fn in moment_files:
    if os.path.exists(fn):
        shutil.rmtree(fn)
        print(f"Removed existing directory: {fn}")

immoments(imagename=filename,
          moments=[0, 1, 8],
          chans='420~630',
          includepix=[0.03, 100],
          region=mask_file,
          outfile=filename.replace('.fits', '.moment'))

# CASA alyways outputs images in CASA image format
# Use task exportfits to make them in FITS files
for mom in moment_files:
    exportfits(imagename=mom,
               fitsimage=mom + '.fits',
               dropdeg=True,
               overwrite=True)

# Plot the moments
fig, axes = plt.subplots(1, 3, figsize=(12, 18))

titles = ['Moment 0 (Integrated intensity)',
          'Moment 8 (Maximum intensity)',
          'Moment 1 (Weighted velocity)']

colourmaps = ['inferno', 'inferno', 'seismic']

for i, (mom_fits, ax, cmap) in enumerate(zip(moment_files, axes, colourmaps)):
    with fits.open(mom_fits + '.fits') as hdul:
        data = hdul[0].data
        im = ax.imshow(data, origin='lower', cmap=cmap)
        ax.set_title(titles[i])
        fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'moment_maps.png'), bbox_inches='tight', dpi=300)
plt.close()