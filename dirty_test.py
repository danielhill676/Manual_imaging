from casatools import table
# input vis name here ###########
vis = 'NGC1365_12m_co21.ms'
#################################
name = vis + '._test'
column = 'data'  # or 'corrected' if you have CORRECTED_DATA column
interactive = False
niter = 0  # dirty image
cell = '0.076arcsec'
imsize = [600,460]
field = 'NGC1365'

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
    cell=cell,
    imsize=imsize,
    weighting='natural',
    gridder=gridder,
    niter=niter,
    datacolumn=column,
    calcpsf=True,
    calcres=True,
    restoration=False, # doesn't make restored image for dirty
    pbcor=False,
    interactive=interactive
)

# if continuum

tclean(
    vis=vis,  # your continuum MS
    imagename=name,
    field='',               # all fields
    spw='',                 # all SPWs
    specmode='mfs',         # continuum
    niter=niter,                # no cleaning, just check if data reads
    imsize=imsize,         # small image to save time
    cell=cell,       # coarse cell
    weighting='natural',    # maximize sensitivity
    interactive=interactive       # no GUI
)
