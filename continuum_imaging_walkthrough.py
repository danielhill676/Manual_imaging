import os

# File and target configuration
vis_file = 'uid___A002_X1003af4_Xa540.ms'  # Change name if necessary
target_name = 'PN_Hb_5'
target_spw = '25,27,29,31'

# Define paths
split_vis = f"{vis_file}.split.cal.target"
image_basename = target_name

# Continuum channels (line-free)
CONT_CHANNELS = ('0:226.2010990572~226.2665287447GHz;226.3939701509~226.5687748384GHz;226.9725834322~227.1312748384GHz,'
                 '1:230.0764655102~230.1931647289GHz;230.2170905102~230.3958014477GHz;230.6809576977~230.8245123852GHz;230.8484381664~231.0046881664GHz,'
                 '2:240.9498861665~240.9567221040GHz;241.1676596040~241.3034017915GHz;241.3991049165~241.4098471040GHz;241.5182455415~241.6334799165GHz;241.7633627290~241.8004721040GHz;241.9303549165~241.9967611665GHz;242.0905111665~242.1666830415GHz;242.3834799165~242.4996908540GHz;242.5993002290~242.7106283540GHz,'
                 '3:243.9977299165~244.0231205415GHz;244.1656986665~244.2057377290GHz;244.2233158540~244.2467533540GHz;244.3033939790~244.3658939790GHz;244.4069096040~244.4615971040GHz;244.5074955415~244.5846439790GHz;244.6246830415~244.6715580415GHz;244.6949955415~244.7496830415GHz;244.7819096040~244.7897221040GHz;244.7994877290~244.8043705415GHz;244.8307377290~244.8385502290GHz;244.8492924165~244.8551517915GHz;244.8990971040~244.9108158540GHz;244.9879642915~245.0768314790GHz;245.1110111665~245.1530033540GHz;245.2487064790~245.2731205415GHz;245.3453861665~245.3893314790GHz;245.4274174165~245.5621830415GHz;245.6285892915~245.6823002290GHz;245.7545658540~245.8151127290GHz;245.8561283540~245.8619877290GHz'
                 )
NITER = 1000000  ### <- Update this for full clean ###
THRESHOLD = '1.17mJy'  ### <- Update this for full clean ###
ROBUST = 0.5  # Feel free to play with the robust parameter to see how it affects the image
"""
If you're on a Mac with an M-chip, you probably want to set this to False
"""
INTERACTIVE = False  # Set to False if you want to run the imaging non-interactively
MASKTYPE = 'auto-multithresh' # Set to 'auto-multithresh' if you want to use the auto-masking feature 

# Imaging parameters for continuum
tclean_params = {
    'imsize': [320, 300],
    'cell': ['0.22arcsec'],
    'phasecenter': 'ICRS 17:47:56.2008 -029.59.39.588',
    'gridder': 'mosaic',
    'deconvolver': 'hogbom',
    'robust': ROBUST,
    'pbcor': True,
    'niter': NITER,
    'usemask': MASKTYPE,
    'interactive': INTERACTIVE,
    'specmode': 'mfs',
    'spw': CONT_CHANNELS,
    'threshold': THRESHOLD,
    'weighting': 'briggs'
}

# --- Step-by-step walkthrough ---

## 1. Split out the target data and science spectral windows

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
          datacolumn=column)

## 2. Make continuum images

# 2.1 Make dirty continuum image

continuum_name = f"{image_basename}.continuum"
dirty_continuum_name = f"{continuum_name}.dirty"

if not os.path.isdir(f"{dirty_continuum_name}.image"):
    print(f"Creating dirty continuum image: {dirty_continuum_name}")
    
    tb.open(split_vis)
    colnames = tb.colnames()
    tb.close()
    column = 'corrected' if 'CORRECTED_DATA' in colnames else 'data'
    print(f"Using {column.upper()} column for tclean")
    
    # Create a modified parameter dict for dirty image to override clean parameters
    dirty_params = tclean_params.copy()
    dirty_params.update({
        'niter': 0,
        'usemask': None,
        'threshold': None,
        'interactive': False
    })
    
    tclean(vis=split_vis,
           imagename=dirty_continuum_name,
           selectdata=True,
           datacolumn=column,
           **dirty_params)

# 2.2 Make clean continuum image

if not os.path.isdir(f"{continuum_name}.image"):
    print(f"Creating clean continuum image: {continuum_name}")
    
    tb.open(split_vis)
    colnames = tb.colnames()
    tb.close()
    column = 'corrected' if 'CORRECTED_DATA' in colnames else 'data'
    print(f"Using {column.upper()} column for tclean")

    tclean(vis=split_vis,
           imagename=continuum_name,
           selectdata=True,
           datacolumn=column,
           **tclean_params)