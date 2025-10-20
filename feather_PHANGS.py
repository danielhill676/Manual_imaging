import os, shutil, glob

lowres = 'NGC3351.fits'
highres = 'NGC3351_12m_co21.image'
regridname = lowres.replace('.fits','regrid.image')
highresnostokes = highres.replace('.image','_nostokes.image')
feathername = lowres.replace('.fits','_feather.image')
jybeamname = lowres.replace('.fits','_jyperbeam.image')

print("Current working directory:", os.getcwd())
print("Files in working directory (first 200 entries):")
for i,f in enumerate(sorted(os.listdir('.'))):
    if i<200:
        print(" ", f)

# helper to check existence covering both files and CASA image directories
def exists(path):
    return os.path.exists(path)

def is_casa_image(path):
    return path.endswith('.image') and os.path.isdir(path)

# check inputs exist
for fname in [lowres, highres]:
    if exists(fname):
        if is_casa_image(fname):
            print(f"Found CASA image directory: {fname} (dir)")
        else:
            print(f"Found file: {fname}")
    else:
        # try fuzzy match using glob (case-insensitive-ish)
        matches = glob.glob('*' + os.path.basename(fname) + '*')
        print(f"Did not find exact entry for {fname}. Glob matches: {matches}")

# remove targets if present (use rmtree for .image directories)
for fname in [regridname, highresnostokes, feathername, jybeamname]:
    if exists(fname):
        try:
            if is_casa_image(fname):
                print(f"Removing CASA image directory: {fname} (this may take a moment)...")
                shutil.rmtree(fname)
            else:
                print(f"Removing file: {fname}")
                os.remove(fname)
        except Exception as e:
            print(f"Could not remove {fname}: {e}")
    else:
        print(f"Not present (skipping): {fname}")


# Remove the degenerate Stokes axis
imsubimage(
    imagename=highres,                     # your high-res cube
    outfile=highresnostokes,
    chans='',
    stokes='I',                                # pick Stokes I plane
    dropdeg=True,                              # this removes the 3rd axis (Stokes)
    overwrite=True
)

import math
nu = 2.30538e11
c = 299792458.0
k = 1.380649e-23
# lowres beam in arcsec (replace if different)
bmaj = imhead(lowres, mode='get', hdkey='bmaj')['value']
print('BMAJ (arcsec)=', bmaj)
bmin = imhead(lowres, mode='get', hdkey='bmin')['value']
print('BMIN (arcsec)=', bmin)
rad = math.pi/(180.0*3600.0)
theta_maj = bmaj*rad
theta_min = bmin*rad
omega = math.pi/(4.0*math.log(2.0)) * theta_maj * theta_min
factor = 2.0*k/((c/nu)**2) * 1e26 * omega
print('Jy/beam per K =', factor)

immath(imagename=lowres,
       expr='IM0 * {factor}'.format(factor=factor),
       outfile=jybeamname
)
imhead(jybeamname, mode='put', hdkey='bunit', hdvalue='Jy/beam')

imregrid(imagename=jybeamname,   # or 'low_orig.image' if reframe not needed
         template=highres,
         output=regridname,
         axes=[0,1,2],                # force spatial+spectral axes to be matched
         interpolation='linear',
         overwrite=True)

imhead(highresnostokes, mode='put', hdkey='crpix3', hdvalue=1.0)
imhead(regridname, mode='put', hdkey='crpix3', hdvalue=1.0)


feather(
    imagename=feathername,
    highres=highresnostokes,
    lowres=regridname
)