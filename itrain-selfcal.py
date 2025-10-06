# I-TRAIN with the European ARC Network
# Improving image fidelity through self-calibration
# Richards, Moravec, Perez-Sanchez, Toribio
# 2021/05/24
#==============================================================================================

#-- This script was prepared for a tutorial on improving image fidelity through self-calibration
#-- of ALMA data.

#-- The dataset used in this script corresponds to Science Verification data of VY CMA in Band 7
#-- that was preprocessed for the purpose
#-- Link to the original dataset VY CMA from SV data: 
#-- https://almascience.nrao.edu/almadata/sciver/VYCMaBand7/
#-- Preparation steps applied the original dataset can be found commented below.

#-- The script is organized in STEPS. Search for the variable "step_title" that holds all the steps
#-- executed by the script. To execute the full script, use the CASA command 'execfile'. If you 
#-- would only like to run certain steps, define the variable 'mysteps' and then execute the script
#-- with 'execfile', eg, to execute steps 0 & 1 only, run in the CASA terminal:
#-- mysteps=[0,1]
#-- execfile('itrain-selfcal.py')

#-- In this script self-calibration of continuum is performed in several cycles, and it is applied 
#-- to all channels/spw.


"""
# PREPARATION ALREADY DONE ON THE ORIGINAL SV DATASET
# For future reference, we include here the steps that have been applied to the orignal 
# Science Verification dataset of VY CMA to obtain the measurement set used in this tutorial.
# Starting with the Science Verification data of VY CMA
# Select just one of the 3 EBs, with Tsys, WVR, bandpass and phase-ref
# solutions applied.  Shift to constant LSR velocity and average target

basename=['uid___A002_X6d0f96_X1de2']
refantenna='DV15'

vislist=[]
for name in basename:
    for s in range(2):  
        os.system('rm -rf '+name+'spw'+str(s)+'_VYCMa_325.ms')
        mstransform(vis=name+'.ms.split',
        spw=str(s),
        outputvis=name+'spw'+str(s)+'_VYCMa_325.ms',datacolumn='corrected',
        field='V*',
        regridms=True,
        width=2,
        outframe='lsrk',
        timeaverage=True,
        timebin='12.1s')
        vislist.append(name+'spw'+str(s)+'_VYCMa_325.ms')

concat(vis=vislist,
       concatvis='X1de2_VYCMa_325.ms')

#Flag antenna DA50
flagdata(vis=vis,
        mode='manual',
        antenna='DA50',
        flagbackup=True)


"""

#============================================================================
# Python packages
#============================================================================

import os
import sys
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np


#============================================================================
# PARAMETER DEFINITION
#============================================================================ 


#-- Measurement set
visname = '7582_selfcal.ms'
plot_prefix = visname+ '_spw'
vis = visname
field = 'NGC7582'

#-- reference antenna
refantenna ='DV14'

#-- Continuum channels selection
#-- This selection was made after calibration on the maser and imaging but 
#-- you could make a selection from the visibility spectrum before self-cal                                         
#contchans='4:0~3,14:0~0,15:0~0,16:0~0,25:0~276;613~959,27:0~239,28:0~0,29:0~239,30:0~0,31:0~239,32:0~0'     

contchans='0:166~194;304~475,1:50~172;216~356;428~436'
cell='0.018arcsec'
imsize=2304

#===========================================================================
# FUNCTIONS
#===========================================================================

# Useful functions for our purposes are defined here
 
def get_im_stats(im_name):
    # Calculate image statistics
    noise=imstat(imagename=im_name,region='ellipse[[1142pix,632pix],[253pix,708pix],0deg]')['rms'][0] #region representative of the image RMS, large enough and without source signal
    peak=imstat(imagename=im_name,region='ellipse[[1198pix,1313pix],[422pix,366pix],0deg]')['max'][0] #region including our target
    print('rms {0:.3f}, peak {1:.3f}, snr {2:.0f}'.format(noise, peak, peak/noise))


def plot_gaincal_table(caltable):
    # make plots for antenna triplets that will be saved in png files
    for ant in range(3):
        plotms(vis=caltable,xaxis='time',yaxis='phase',iteraxis='antenna',gridrows=4, gridcols=2,coloraxis='spw', 
            xaxisfont = 7, yaxisfont = 7, highres =  True, 
            antenna=str(ant*8)+'~'+str(ant*8+7), 
            showgui = False, plotfile=caltable+'_ant'+str(ant*8)+'-'+str(ant*8+7)+'.png')

def plot_gaincal_snr_dist(path, visname, selfcal_cycle, solints):
    # Make plot of SNR of gaintables for a list of solution intervals
    # path: path to caltables; plot will be saved in that path too
    # solints: list of solution intervals (list of strings)
    print("Percentile of score for snr=6, i.e. percentage of solutions with snr <=6:")
    for ss in solints:
        tb.open(path+'/'+visname+'.'+selfcal_cycle+'.solint_'+ss+'.tb')
        snr = tb.getcol( 'SNR' ).ravel()
        tb.close()
      
        if(matplotlib.__version__>="2.0"):
            plt.hist( snr, bins=50, density=True, histtype='step', label=ss )
        else:
            plt.hist( snr, bins=50, normed=True, histtype='step', label=ss )

        print( 'P(<=6) = {0}  ({1})'.format( stats.percentileofscore( snr, 6 ), ss ) )

    plt.legend( loc='upper right' )
    plt.xlabel( 'SNR' )
    plt.savefig(path+'/'+selfcal_cycle+'_SNR_hist_solint_all.png')



#===========================================================================
# STEPS
#===========================================================================

# List of steps executed by this script
thesteps=[]
step_title = {0:  'List the data set and plot antennas and visibility spectrum',
              1:  'Make dirty image of continuum',
	      ### INITIAL MODEL 
              2:  'Make an initial, conservative cleaning',
              3:  'Check and save model',
              ### FIRST ROUND OF SELF-CALIBRATION - PHASE
              4:  'Calculate gain solution table - phase-only, solution interval = scan-length',
              5:  'Explore different solution intervals', 
              6:  '[ADVANCED] Calculate SNR of the different solution intervals', 
              7:  'Apply calibration table',
              8:  'Make second, conservative cleaning and save model',
              ### SECOND ROUND OF SELF-CALIBRATION - PHASE
              9:  'Explore different solution intervals',
              10: '[ADVANCED] Calculate SNR of the different solution intervals',
              11: 'Calculate gain solution table - phase-only, solution interval = 60s applying round 1 table on-the-fly',
              12: 'Apply calibration tables',
              13: 'Make image of continuum and save model',
              ### THIRD ROUND OF SELF-CALIBRATION - AMPLITUDE & PHASE
              14: 'Calculate gain solution table - amplitude and phase, long solution interval',
              15: 'Apply calibration tables',
              16: 'Make image of continuum and save model',
              ### FOURTH ROUND OF SELF-CALIBRATION - AMPLITUDE & PHASE
              17: 'Calculate gain solution table - amplitude and phase, short solution interval',
              18: 'Apply calibration table',
              ### FINAL CONTINUUM IMAGE
              19: 'Make image of continuum and save model', 
	     }
 
# The Python variable 'mysteps' will control which steps are executed when you start the script using
#   execfile('[This script].py')
# e.g. to execute only steps 2, 3, and 4 of the script define first 'mysteps' and then execute the script:
#   mysteps = [2,3,4]
#   execfile('script.py')
# Setting mysteps = [] will make it execute all steps.

try:
  print('List of steps to be executed ...', mysteps)
  thesteps = mysteps
except:
  print('global variable mysteps not set.')
if (thesteps==[]):
  thesteps = list(range(0,len(step_title)))
  print('Executing all steps: ', thesteps)

#---------
mystep = 0 
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  #print listobs
  os.system('rm '+visname+'_listobs.txt')
  listobs(vis=vis, listfile=visname+'_listobs.txt', verbose=True)

  #plot antenna locations
  os.system('rm '+visname+'*plotants*png')
  plotants(vis, figfile=visname+'_plotants.png')
  plotants(vis, logpos=True, figfile=visname+'_plotants_log.png')

  # Plots of amplitude vs. frequency 
  # --You will notice the very bright maser line in spw 0:981
  # --even for continuum it is worth plotting this to reveal bad data
  for s in ['0','1']:
      os.system('rm -rf '+visname+'_spw'+s+'_vis-spectrum.png')
      plotms(vis = vis,
             xaxis='frequency', yaxis='amp',
             selectdata=True, spw=s,
             avgtime='1e8',avgscan=True,avgbaseline=True,
             coloraxis='baseline',
             showgui = False,
             highres = True,
             plotfile=visname+'_spw'+s+'_vis-spectrum.png')
  
 

#---------
mystep = 1 
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])


  ## Make a first dirty imaging of the continuum to get a sense of the structure of the object
  imagename = visname + '_cont.dirty'
  os.system('rm -rf '+imagename+'.*')
  tclean(vis = vis,
        imagename = imagename,
        field = field,
        spw=contchans,
        specmode='mfs',
        cell=cell,
        imsize=imsize,
        niter=0,
        interactive=False)

  # view image
  # imview(imagename+'.image')



### INITIAL MODEL
#---------
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  #Delete existing models and calibration tables to make sure we have a neat start:
  #delmod(vis=vis, scr=True)
  #clearcal(vis=vis)

  # make an initial, conservative clean 
  imagename = visname + '_cont0.init.clean'
  os.system('rm -rf '+imagename+'.*')
  tclean(vis = vis,
        imagename = imagename,
        field = field,
        spw=contchans,
        specmode='mfs',
        cell=cell,
        imsize=imsize,
        niter = 200,
        interactive=False, usemask='user', mask='7582_cont_cleanmask.mask')
  #Get image statistics for comparison 
  get_im_stats(imagename+'.image')


#---------
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  #model name
  modelname=visname+'_cont0.init.clean.model'

  # check that model has saved
  plotms(vis=vis, xaxis='UVwave', yaxis='amp', ydatacolumn='model',showgui=False,plotfile=modelname+'.png')

  # force model to save
  ft(vis=vis,model=modelname,usescratch=True)

  # check that model has saved after ft
  plotms(vis=vis, xaxis='UVwave', yaxis='amp', ydatacolumn='model',showgui=False,plotfile=modelname+'_ft.png')



### FIRST ROUND OF SELF-CALIBRATION - PHASE
#---------
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  # calculate gain table for solint='inf' = scan length
  solint='inf'
  caltable=visname+'_cont.ph1.solint_'+solint+'.tb'
  os.system('rm -rf '+caltable)
  gaincal(vis = vis,
          field= field,
          refant=refantenna,
          caltable=caltable,
          spw=contchans,
          calmode='p',
          solint=solint,
          gaintype='G',
          minsnr=3)


  # view gain table for all antennas
  plotms(vis=caltable,xaxis='time',yaxis='phase',iteraxis='antenna',gridrows=3, gridcols=3,coloraxis='spw')


#---------
mystep = 5
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print('Step ', mystep, step_title[mystep])

    ### COMMANDS FOR EXPLORING THE SOLUTION INTERVAL
    # # We recommend investigating the phases with various solution intervals to get a sense of the interval on which the solutions vary

    # The following loop calculates gaincal solutions for a list of intervals and makes corresponding plots
    # The output is saved in a separate folder  
    selfcal_cycle = 'ph1_checks'
    solint_all = ['int', '20s', '40s', '60s', '80s', '160s', 'inf'] 
    for solint in solint_all:
        print('Solint:', solint)
        caltable = visname+'.'+selfcal_cycle+'.solint_'+solint+'.tb'
        gaincal(vis=vis,caltable=caltable,solint=solint,refant=refantenna,spw=contchans,calmode='p',gaintype='G',minsnr=3)

        # make plots for antenna triplets that will be saved in png files
        plot_gaincal_table(caltable)
    
    os.system('rm -r '+selfcal_cycle)
    if not os.path.exists(selfcal_cycle):
        os.makedirs(selfcal_cycle)
    os.system('mv *.'+selfcal_cycle+'*tb* '+selfcal_cycle+'/')
    print("Check output of this step in folder: "+selfcal_cycle)

#---------
mystep = 6
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print('Step ', mystep, step_title[mystep])

    # [ADVANCED]
    # Calculate the distribution of SNR of the gaincal tables for different solution intervals
    # using the tables generated in the previous step; plots are moved to the same output folder as in the step above
    selfcal_cycle = 'ph1_checks'
    solint_all = ['int', '20s', '40s', '60s', '80s', '160s', 'inf']
    
    try:
        path = selfcal_cycle + '/'
        plot_gaincal_snr_dist(path, visname, selfcal_cycle, solint_all)
        print("Check output of this step in folder: "+selfcal_cycle)
    except:
        print("You need to run the previous step for the same list of solution intervals")


#---------
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])
 
  # apply the solutions to the MS  
  caltable=visname+'_cont.ph1.solint_inf.tb'
  applycal(vis = vis,
           field= field,
           spw='0,1',
           spwmap=[0,1],
           gaintable=caltable,
           calwt = False,
           applymode='calonly',
           flagbackup = False)



#---------
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])


  # make a second, conservative clean
  imagename = visname + '_cont.ph1.clean'
  os.system('rm -rf '+imagename+'.*')
  tclean(vis = vis,
        imagename = imagename,
        field = field,
        spw=contchans,
        specmode='mfs',
        cell=cell,
        imsize=imsize,
        niter=200,
        interactive=False, usemask='user', mask='7582_cont_cleanmask.mask')
  #Get image statistics for comparison 
  get_im_stats(imagename+'.image')


  #force model to save
  modelname=imagename+'.model'
  ft(vis=vis,model=modelname,usescratch=True)



### SECOND ROUND OF SELF-CALIBRATION - PHASE
#---------
mystep = 9
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print('Step ', mystep, step_title[mystep])

    ### COMMANDS FOR EXPLORING THE SOLUTION INTERVAL
    # # We recommend investigating the phases with various solution intervals to get a sense of the interval on which the solutions vary

    # The following loop calculates gaincal solutions for a list of intervals and makes corresponding plots
    # The output is saved in a separate folder  
    selfcal_cycle = 'ph2_checks'
    solint_all = ['int', '20s', '40s', '60s', '80s', '160s', 'inf']
    for solint in solint_all:
        print('Solint:', solint)
        solint_1='inf'
        caltable = visname+'.'+selfcal_cycle+'.solint_'+solint+'.tb'
        gaincal(vis=vis,caltable=caltable,solint=solint,refant=refantenna,spw=contchans,
            gaintable = [visname + '_cont.ph1.solint_'+solint_1+'.tb'], spwmap=[0,1], calmode='p',gaintype='G',minsnr=3)

        # make plots for antenna triplets that will be saved in png files
        plot_gaincal_table(caltable)

    os.system('rm -r '+selfcal_cycle)
    if not os.path.exists(selfcal_cycle):
        os.makedirs(selfcal_cycle)
    os.system('mv *.'+selfcal_cycle+'*tb* '+selfcal_cycle+'/')
    print("Check output of this step in folder: "+selfcal_cycle)



#---------
mystep = 10
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print('Step ', mystep, step_title[mystep])

    # [ADVANCED]
    # Calculate the distribution of SNR of the gaincal tables for different solution intervals
    # using the tables generated in the previous step; plots are moved to the same output folder as in the step above
    selfcal_cycle = 'ph2_checks'
    solint_all = ['int', '20s', '40s', '60s', '80s', '160s', 'inf']
   
    try:
        path = selfcal_cycle + '/'
        plot_gaincal_snr_dist(path, visname, selfcal_cycle, solint_all)
        print("Check output of this step in folder: "+selfcal_cycle)
  
    except:
        print("You need to run the previous step for the same list of solution intervals")



#---------
mystep = 11
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  # calculate gain table for solint='60s' while applying the table from round 1 'on the fly'
  solint_1='inf'
  solint='60s'
  caltable = visname + '_cont.ph2.solint_'+solint+'.tb'
  os.system('rm -rf '+caltable)
  gaincal(vis = vis,
          field= field,
          refant=refantenna,
          caltable=caltable,
          spw=contchans,
          gaintable = [visname + '_cont.ph1.solint_'+solint_1+'.tb'],
          spwmap=[0,1],
          calmode='p',
          solint=solint,
          gaintype='G',
          minsnr=3)


  # view gain table for all antennas
  plotms(vis=caltable,xaxis='time',yaxis='phase',iteraxis='antenna',gridrows=3, gridcols=3,coloraxis='spw')



#---------
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])
 
  # apply the cumulative solutions to the MS 
  solint_1='inf'
  solint_2='60s' 
  applycal(vis = vis,
           field= field,
           spw='0,1',
           gaintable=[visname+'_cont.ph1.solint_'+solint_1+'.tb', visname+'_cont.ph2.solint_'+solint_2+'.tb'],
           spwmap=[[0,1],[0,1]],
           calwt = False,
           applymode='calonly',
           flagbackup = False)



#---------
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  # If you would like to compare what the second round of phase self-cal accomplished compared to the first
  # make a third, conservative clean
  imagename = visname + '_cont.ph2.clean'
  os.system('rm -rf '+imagename+'.*')
  tclean(vis = vis,
      imagename=imagename,
      field=field,
      spw=contchans,
      specmode='mfs',
      cell=cell,
      imsize=imsize,
      niter=200,
      interactive=False, usemask='user', mask='7582_cont_cleanmask.mask')
  # get image statistics for comparison
  get_im_stats(imagename+'.image')


  #force model to save
  modelname=imagename+'.model'
  ft(vis=vis,model=modelname,usescratch=True)


### THIRD ROUND OF SELF-CALIBRATION - AMPLITUDE & PHASE
#---------
mystep = 14
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  # long solution interval, applying previous solutions on the fly
  solint='120s'
  caltable = visname + '_cont.ap1.solint_'+solint+'.tb'
  # apply the cumulative solutions to the MS 
  solint_1='inf'
  solint_2='60s'
  os.system('rm -rf '+caltable)
  gaincal(vis = vis,
          field= field,
          refant=refantenna,
          caltable=caltable,
          gaintable = [visname + '_cont.ph1.solint_'+solint_1+'.tb', visname + '_cont.ph2.solint_'+solint_2+'.tb' ], 
          spwmap=[[0,1],[0,1]],
          spw=contchans,
          calmode='ap',
          solint=solint,
          gaintype='G',
          minsnr=3)


  # view gain table for all antennas
  plotms(vis=caltable,xaxis='time',yaxis='phase',iteraxis='antenna',gridrows=3, gridcols=3,coloraxis='spw')





#---------
mystep = 15
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  # apply the cumulative solutions to the MS 
  solint_1='inf'
  solint_2='60s'
  solint_3='120s'
  applycal(vis = vis,
           field= field,
           spw='0,1',
           spwmap=[[0,1],[0,1],[0,1]],
           gaintable=[visname+'_cont.ph1.solint_'+solint_1+'.tb', visname+'_cont.ph2.solint_'+solint_2+'.tb', visname+'_cont.ap1.solint_'+solint_3+'.tb'],
           calwt = False,
           applymode='calonly',
           flagbackup = False)



#---------
mystep = 16
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  # make yet another, conservative clean
  imagename = visname + '_cont.ap1.clean'
  os.system('rm -rf '+imagename+'.*')
  tclean(vis = vis,
      imagename=imagename,
      field=field,
      spw=contchans,
      specmode='mfs',
      cell=cell,
      imsize=imsize,
      niter=300,
      interactive=False, usemask='user', mask='7582_cont_cleanmask.mask')
  # get image statistics for comparison
  get_im_stats(imagename+'.image')


  #force model to save
  modelname=imagename+'.model'
  ft(vis=vis,model=modelname,usescratch=True)



### FOURTH ROUND OF SELF-CALIBRATION - AMPLITUDE & PHASE
#---------
mystep = 17
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  # calculate gain table for solint='60s' while applying previous solutions 'on the fly'
  solint='60s'
  caltable = visname + '_cont.ap2.solint_'+solint+'.tb'
  solint_1='inf'
  solint_2='60s'
  solint_3='120s'
  os.system('rm -rf '+caltable)
  gaincal(vis = vis,
          field= field,
          refant=refantenna,
          caltable=caltable,
          spwmap=[[0,1],[0,1],[0,1]],
          gaintable=[visname+'_cont.ph1.solint_'+solint_1+'.tb', visname+'_cont.ph2.solint_'+solint_2+'.tb', visname+'_cont.ap1.solint_'+solint_3+'.tb'],
          spw=contchans,
          calmode='ap',
          solint=solint,
          gaintype='T',    # notice gaintype option
          minsnr=3)


  # view gain table for all antennas
  plotms(vis=caltable,xaxis='time',yaxis='phase',iteraxis='antenna',gridrows=3, gridcols=3,coloraxis='spw')



#---------
mystep = 18
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  # apply the cumulative solutions to the MS
  solint_1='inf'
  solint_2='60s'
  solint_3='120s' 
  solint_4='60s'
  applycal(vis = vis,
           field= field,
           spw='0,1',
           spwmap=[[0,1],[0,1],[0,1],[0,1]],
           gaintable=[visname+'_cont.ph1.solint_'+solint_1+'.tb', visname+'_cont.ph2.solint_'+solint_2+'.tb', 
           visname+'_cont.ap1.solint_'+solint_3+'.tb', visname+'_cont.ap2.solint_'+solint_4+'.tb'],
           calwt = False,
           flagbackup = False, applymode='calfalg')



### FINAL CONTINUUM IMAGE
#---------
mystep = 19
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  imagename = visname + '_cont.ap2.clean'
  os.system('rm -rf '+imagename+'.*')
  tclean(vis = vis,
      imagename=imagename,
      field=field,
      spw=contchans,
      specmode='mfs',
      cell=cell,
      imsize=2304,
      niter=300,
      deconvolver = 'multiscale', 
      scales=[0,4,8,12],
      interactive=False, usemask='user', mask='7582_cont_cleanmask.mask')
  # get image statistics for comparison
  get_im_stats(imagename+'.image')


  #force model to save
  modelname=imagename+'.model'
  ft(vis=vis,model=modelname,usescratch=True)





#---------------------------------------------------------------------------------
#----- End of script.
#----------------------------------------------------------------------------------

