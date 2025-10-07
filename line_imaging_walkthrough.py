import os

# File and target configuration
vis_file = 'uid___A002_Xb945f7_X1b14.ms.split.cal'
target_name = 'NGC7582'
target_spw = '0'


# Define paths
split_vis = f"{vis_file}.target"
contsub_vis = f"{vis_file}.target.contsub"
image_basename = 'ngc7582'

"""
Continuum channels (line-free)
You will need to manually set this to specific the continuum channels for each spectral window
Note that SPWs 3 and 4 are line-free and can be used as-is. Only SPWs 1 and 2 channels need to be specified.
[You should have already done this for yesterday's continuum imaging exercise, so you can reuse the same values here]
"""
# Continuum channels (line-free) NOTE: Usually after splitting it just contains target spw which is then spw 0
CONT_CHANNELS = ('0: 229.029325914~229.138712154GHz; 229.572350465~230.236481211GHz' #, 1: 230.030519105~230.514943969GHz ; 230.675116706~231.23376796GHz ; 231.507233609~231.538486826GHz'
                 )
NITER = 100000  ### <- Update this for full clean ### 0 for dirty, many for full clean
THRESHOLD = '2.1mJy'  ### <- Update this for full clean ### 0 for dirty, then use to get for full clean
ROBUST = 0.5  # Feel free to play with the robust parameter to see how it affects the image
INTERACTIVE = False  # Set to True if you want to run the imaging interactively
MASKTYPE = 'auto-multithresh'

"""
To image specific lines, you will need to specify the start channel, channel width, and number of channels for each chunk.
You can also use units of frequency or velocity instead of channel numbers (see tclean documentation).
The below numbers are just placeholders and will need to be updated based on the spectral line(s) in the data.
"""
LINE_CHUNKS = [
    {'start': 186, 'width': 1, 'nchan': 113}
]

IMSIZE = 1152

from casatools import synthesisutils
su = synthesisutils()
size = su.getOptimumSize(IMSIZE)

# Imaging parameters
tclean_params = {
    'imsize': [size, size],
    'cell': ['0.04arcsec'],
    'phasecenter': 'ICRS 23:18:23.60 -42.22.14.00000', 
    'gridder':'mosaic',
    'deconvolver': 'multiscale',
    'robust': ROBUST,
    'pbcor': True,
    'niter': NITER,
    'usemask': MASKTYPE,
    'interactive': INTERACTIVE,
    'specmode': 'cube',
    'spw': '0',
    'threshold': THRESHOLD,
    'weighting': 'briggsbwtaper',
    'restoringbeam': 'common',
    'minbeamfrac': 0.3, # default is 0.3 reduce when automasking is not working well
    'noisethreshold': 5.0 # default is 5.0 reduce when automasking is not working well
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
    line_name = f"{image_basename}.spw0.chunk{chunk_idx}"
    
    if not os.path.isdir(f"{line_name}.image"):
        print(f"Creating line cube for chunk {chunk_idx}")
        
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
        column = 'corrected' if 'CORRECTED_DATA' in colnames else 'data'
        print(f"Using {column.upper()} column for tclean")
        
        tclean(vis=contsub_vis,
               imagename=line_name,
               selectdata=True,
               datacolumn=column,
               **chunk_params)