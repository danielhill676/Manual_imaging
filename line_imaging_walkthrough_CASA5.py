import os

# File and target configuration
vis_file = 'uid___A002_Xb945f7_X1b14.ms.split.cal'
# Define target
target_name = 'NGC7582'
target_spw = '0'

# Define paths
split_vis = vis_file + '.target'
contsub_vis = vis_file + '.target.contsub'
image_basename = 'ngc7582'

# Continuum channels (line-free)
CONT_CHANNELS = ('0: 229.029325914~229.138712154GHz; 229.572350465~230.236481211GHz',)

NITER = 100000
THRESHOLD = '2.1mJy'
ROBUST = 0.5
INTERACTIVE = True
MASKTYPE = 'auto-multithresh'

LINE_CHUNKS = [{'start': 186, 'width': 1, 'nchan': 113}]
IMSIZE = 80

# CASA 5.7: synthesisutils via casac
su = casac.synthesisutils()
size = su.getOptimumSize(IMSIZE)

# Imaging parameters
tclean_params = {
    'imsize': [size, size],
    'cell': ['0.5arcsec'],
    'phasecenter': 'ICRS 23:18:23.60 -42.22.14.00000',
    'gridder': 'mosaic',
    'deconvolver': 'multiscale',
    'robust': ROBUST,
    'pbcor': True,
    'niter': NITER,
    'usemask': MASKTYPE,
    'interactive': INTERACTIVE,
    'specmode': 'cube',
    'spw': '0',
    'threshold': THRESHOLD,
    'weighting': 'briggs',
    'restoringbeam': 'common',
    'minbeamfrac': 0.3,
    'noisethreshold': 5.0
}



# --- Step-by-step walkthrough ---

## 1. Split out the target data and science spectral windows
"""
You should have already done this for yesterday's continuum imaging exercise, so you can skip this step.
I just left it in for completeness.
"""

if not os.path.isdir(split_vis):
    print(f"Splitting {vis_file} to {split_vis}")
    
    tb.open(vis_file)
    colnames = tb.colnames()
    tb.close()
    column = 'corrected' if 'CORRECTED_DATA' in colnames else 'data'
    print(f"Using {column.upper()} column for split")

    split(vis=vis_file,
      outputvis=split_vis,
      field=target_name,
      spw=target_spw,
      datacolumn=column,  # will copy the values
      keepflags=True)          # optional, keep flags

## 2. Subtract continuum from the split MS

if not os.path.exists(contsub_vis):
    print(f"Performing continuum subtraction")
    
    uvcontsub(vis=split_vis,
            outputvis=contsub_vis,
            fitspec=CONT_CHANNELS,
            fitorder=0,
            datacolumn='data')

## 3. Make dirty cube of full spectral window

dirty_line_name = f"{image_basename}.spw0.dirty"

if not os.path.isdir(f"{dirty_line_name}.image"):
    print(f"Creating dirty line cube: {dirty_line_name}")
    
    # Create a modified parameter dict for dirty cube
    dirty_params = tclean_params.copy()
    dirty_params.update({
        'niter': 0,
        'usemask': None,
        'threshold': None,
        'interactive': False
    })
    
    tb.open(contsub_vis)
    colnames = tb.colnames()
    tb.close()
    column = 'corrected' if 'CORRECTED_DATA' in colnames else 'data'
    print(f"Using {column.upper()} column for tclean")
    
    tclean(vis=contsub_vis,
           imagename=dirty_line_name,
           selectdata=True,
           datacolumn=column,
           **dirty_params)

## 4. Make clean line images for each chunk
for chunk_idx, chunk in enumerate(LINE_CHUNKS):
    line_name = image_basename + '.spw0.chunk' + str(chunk_idx)
    os.system('rm -rf ' + line_name + '.image')  # remove if exists
    if not os.path.isdir(line_name + '.image'):
        print('Creating line cube for chunk {}'.format(chunk_idx))
        
        # Create a modified parameter dict for this chunk
        chunk_params = tclean_params.copy()
        chunk_params.update({
            'start': chunk['start'],
            'width': chunk['width'],
            'nchan': chunk['nchan']
        })
        
        tb.open(contsub_vis)
        colnames = tb.colnames()
        tb.close()
        
        if 'CORRECTED_DATA' in colnames:
            column = 'corrected'
        else:
            column = 'data'
        
        print('Using {} column for tclean'.format(column.upper()))
        
        tclean(vis=contsub_vis,
               imagename=line_name,
               selectdata=True,
               datacolumn=column,
               **chunk_params)