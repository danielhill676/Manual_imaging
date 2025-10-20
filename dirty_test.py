from casatools import table
# input vis name here ###########
vis = 'NGC3351_12m_co21.ms'
#################################
name = vis + '.dirty_test'
column = 'data'  # or 'corrected' if you have CORRECTED_DATA column

tb = table()
tb.open(vis + '/FIELD')
phase_dir = tb.getcol('PHASE_DIR')
tb.close()

# Flatten safely: each field can have multiple reference directions
is_mosaic = len({tuple(phase_dir[:, i, 0]) for i in range(phase_dir.shape[1])}) > 1
gridder = 'mosaic' if is_mosaic else 'standard'
print('Gridder selected:', gridder)

# --- Quick dirty image test ---
tclean(
    vis=vis,
    imagename=name,
    specmode='cube',
    restfreq='230.538GHz',
    outframe='LSRK',
    nchan=1,
    cell='0.5arcsec',
    imsize=[256, 256],
    weighting='natural',
    gridder=gridder,
    niter=0,
    datacolumn=column,
    calcpsf=True,
    calcres=True,
    restoration=False,
    pbcor=False,
    interactive=False
)

# if continuum

tclean(
    vis=vis,  # your continuum MS
    imagename=name,
    field='',               # all fields
    spw='',                 # all SPWs
    specmode='mfs',         # continuum
    niter=0,                # no cleaning, just check if data reads
    imsize=[64,64],         # small image to save time
    cell='1.0arcsec',       # coarse cell
    weighting='natural',    # maximize sensitivity
    interactive=False       # no GUI
)
