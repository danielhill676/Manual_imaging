import os, shutil, glob

lowres = 'NGC4254.fits'
highres = 'NGC4254_12m_co21_pbcorr_trimmed_k.fits'
highres_casa = highres.replace('.fits','.image')
regridname = lowres.replace('.fits','regrid.image')
highresnostokes = highres.replace('.fits','_nostokes.image')
feathername = lowres.replace('.fits','_feather.image')
featherfits = feathername.replace('.image','.fits')
jybeamname_low = lowres.replace('.fits','_jyperbeam.image')
jybeamname_high = highres.replace('.fits','_jyperbeam.image')

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
for fname in [regridname, highresnostokes, feathername, featherfits, highres_casa]:
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


importfits(fitsimage=highres,imagename=highres_casa)

# Remove the degenerate Stokes axis
imsubimage(
    imagename=highres_casa,                     # your high-res cube
    outfile=highresnostokes,
    chans='',
    stokes='I',                                # pick Stokes I plane
    dropdeg=True,                              # this removes the 3rd axis (Stokes)
    overwrite=True
)
######
import math
nu = 2.30538e11
c = 299792458.0
k = 1.380649e-23
# lowres beam in arcsec (replace if different)
bmaj_low = imhead(lowres, mode='get', hdkey='bmaj')['value']
print('BMAJ (arcsec)=', bmaj)
bmin_low = imhead(lowres, mode='get', hdkey='bmin')['value']
print('lowres BMIN (arcsec)=', bmin_low)
rad = math.pi/(180.0*3600.0)
theta_maj_low = bmaj_low*rad
theta_min_low = bmin_low*rad
omega_low = math.pi/(4.0*math.log(2.0)) * theta_maj_low * theta_min_low
factor_low = 2.0*k/((c/nu)**2) * 1e26 * omega_low
print('lowres Jy/beam per K =', factor_low)

immath(imagename=lowres,
       expr='IM0 * {factor}'.format(factor=factor),
       outfile=jybeamname_low
)
imhead(jybeamname_low, mode='put', hdkey='bunit', hdvalue='Jy/beam')
#####

######
bmaj_high = imhead(highresnostokes, mode='get', hdkey='bmaj')['value']
print('BMAJ (arcsec)=', bmaj_high)
bmin_high = imhead(highresnostokes, mode='get', hdkey='bmin')['value']
print('highres BMIN (arcsec)=', bmin_high)
theta_maj_high = bmaj_high*rad
theta_min_high = bmin_high*rad
omega_high = math.pi/(4.0*math.log(2.0)) * theta_maj_high * theta_min_high
factor_high = 2.0*k/((c/nu)**2) * 1e26 * omega_high
print('highres Jy/beam per K =', factor_high)

immath(imagename=highresnostokes,
       expr='IM0 * {factor}'.format(factor=factor_high),
       outfile=jybeamname_high
)
imhead(jybeamname_high, mode='put', hdkey='bunit', hdvalue='Jy/beam')
#####



imregrid(imagename=jybeamname_low,
         template=jybeamname_high,
         output=regridname,
         axes=[0,1,2],                # force spatial+spectral axes to be matched
         interpolation='linear',
         overwrite=True)

imhead(jybeamname_high, mode='put', hdkey='crpix3', hdvalue=1.0)
imhead(regridname, mode='put', hdkey='crpix3', hdvalue=1.0)


feather(
    imagename=feathername,
    highres=jybeamname_high,
    lowres=regridname
)

exportfits(imagename=feathername,
           fitsimage=featherfits,
           overwrite=True
)