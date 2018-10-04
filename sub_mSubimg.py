#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#delete pixels in red and nb ccds that don't overlap and save new cdds. Uses mSubimage from Montage
#The create Halpha - r images


from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import glob
import shutil
import sys
import subprocess as sp
import math

import make_lists




args = make_lists.get_vphas_num()
ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'

#returns blocks in the order [u, g, r, r2, i, NB]
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)


#list pairs of nb and (red block) r dirs
pair_list = [ [a_block[2], a_block[5]] , [b_block[2], b_block[5]] ]
      
      
  

for ind,pair in enumerate(pair_list):


	if ind==0:
		block = 'a'
	elif ind==1:
		block = 'b'
       
       
       #make some new directories to store files in
        _, red_name = pair[0].rsplit('/', 1)
        _, nb_name = pair[1].rsplit('/', 1)
        
                
        red_findir = ex_path+'/trimmed_red/' + red_name
        nb_findir = ex_path +'/trimmed_nb/' + nb_name
                
        if not os.path.exists(red_findir):
                os.makedirs(red_findir)
                print "Created final directory: %s" %red_findir

        if not os.path.exists(nb_findir):
                os.makedirs(nb_findir)
                print "Created final directory: %s" %nb_findir

        sub_findir = ex_path+'/Halpha_min_r/no_bin/'+nb_name+red_name
        if not os.path.exists(sub_findir):
                os.makedirs(sub_findir)
                print "Created final directory: %s" %sub_findir
	print



	#loop over the CCDS
        for ccdnum in range(1,33):
        	print 'CCD', ccdnum
        
        	"""
        	#Use confidence corrected files
        	redpath = pair[0] + '/confcorr/'+block+'_ccd'+str(ccdnum)+'.fits'
        	if not os.path.exists(redpath):
        		print 'Path to confidence corrected red ccd could not be found'
        		print redpath
        	
        	        		
              	nbpath = pair[1] + '/confcorr/'+block+'_ccd'+str(ccdnum)+'.fits'
        	if not os.path.exists(nbpath):
        		print 'Path to NB ccd could not be found'
        		print nbpath
        		print 
        	"""
        	

		
		#Use original image files
        	redpath = glob.glob( pair[0] + '/single/*_ccds/*_ccd'+str(ccdnum)+'.fit' )
        	redpath = redpath[0]
        	if not os.path.exists(redpath):
        		print 'Path to r single file could not be found'
        		print redpath
        		raw_input('Paused')

        	
        	nbpath = glob.glob( pair[1] + '/single/*_ccds/*_ccd'+str(ccdnum)+'.fit')
        	nbpath = nbpath[0]
        	if not os.path.exists(nbpath):
        		print 'Path to NB single file could not be found'
        		print nbpath
        		raw_input('Paused')
        	
        
        
        	print 'Input files:'
        	print redpath
        	print nbpath
        	print
        	
        	
        	
        	sub_name = sub_findir+'/'+nb_name+red_name+'_ccd'+str(ccdnum)+'.fits'
                #break if ccd already exists
                if os.path.exists(sub_name):
                        print "Halpha-r ccd already exists"
                        continue
        
        
        
                #make new image files: trim the CCDs to they're the same size                 
                newRed = red_findir+'/'+red_name+'_ccd'+str(ccdnum)+'.fit'
                newNB = nb_findir+'/'+nb_name+'_ccd'+str(ccdnum)+'.fit'
                #break if ccd already exists
                if os.path.exists(newRed) or os.path.exists(newNB) :
                        print "Trimmed red ccd already exists"
                        
                        
                else:                     

                	nb = fits.open(nbpath)
                	red = fits.open(redpath)
                
                	"""
                	# 0 for confidence corrected files, 1 for single
                	red_header = red[0].header
                	nb_header = nb[0].header
                	red_img = red[0].data
                	nb_img = nb[0].data
                	"""
                
                
                	red_header = red[1].header
                	nb_header = nb[1].header
                	red_img = red[1].data
                	nb_img = nb[1].data
                
                
                        red_wcs = wcs.WCS(red_header)
                	nb_wcs = wcs.WCS(nb_header)


                
                	#check for nb wcs shifts relative to red in raw images and shift nb ccd into wcs of the red
                	cx = red_img.shape[0]/2
                	cy = red_img.shape[1]/2
                	cra, cdec = red_wcs.wcs_pix2world(cx, cy, 1) #red as reference
                	refx, refy = nb_wcs.wcs_world2pix(cra, cdec, 1) #pixels in nb with same wcs as cra and cdec in red

                	#shiftx = change of centre col, shifty = change of centre row
                	shiftx = refx-cx
                	if shiftx <1 :
                	        shiftx = int(math.floor(shiftx))
                	else:
                	        shiftx = int(math.ceil(shiftx))
                        
                	shifty = refy-cy
                	if shifty <1:
                	        shifty = int(math.floor(shifty))
                	else:
               	        	shifty = int(math.ceil(shifty))
                        
                        
                	if shiftx<0:
                        	red_x_start=str(-shiftx) 
                        	red_x_size = str(red_img.shape[1]+shiftx)
                        	nb_x_start = '1'
                        	nb_x_size =str(nb_img.shape[1]+shiftx)
                	elif shiftx>=0:
                        	red_x_start='1'
                        	red_x_size = str(red_img.shape[1]-shiftx)
                        	nb_x_start = str(shiftx)
                        	nb_x_size =str(nb_img.shape[1]-shiftx)
                        

                	if shifty>=0:
                        	red_y_start = '1'
                        	red_y_size= str(red_img.shape[0]-shifty)
                        	nb_y_start = str(shifty)
                        	nb_y_size = str(nb_img.shape[0]-shifty)  
                	elif shifty<0:
                        	red_y_start = str(-shifty)
                        	red_y_size= str(red_img.shape[0]+shifty)
                        	nb_y_start = '1'
                        	nb_y_size = str(nb_img.shape[0]+shifty)      
                                         
                        
               
                
               		print 'Calling mSubimage'
               		#sp.call(["/usr/local/Montage_v5.0/mSubimage", "-p", redpath, newRed, red_x_start, red_y_start, red_x_size, red_y_size])
               		#sp.call(["mSubimage", "-p", nbpath, newNB, nb_x_start, nb_y_start, nb_x_size, nb_y_size])
                
               		#hdu 1 for the single images
               		sp.call(["/usr/local/Montage_v5.0/bin/mSubimage", "-p", "-h", "1", redpath, newRed, red_x_start, red_y_start, red_x_size, red_y_size])
               		sp.call(["/usr/local/Montage_v5.0/bin/mSubimage", "-p", "-h", "1", nbpath, newNB, nb_x_start, nb_y_start, nb_x_size, nb_y_size])
               		
               		
               		#close files
               		red.close()
                	nb.close()
                	
                	
                	
                               
                
                if not os.path.exists(newRed):
                	print 'Something went wrong'
                	print 'No trimmed r file'
                	sys.exit()               
                	
                if not os.path.exists(newNB):
                	print 'Something went wrong'
                	print 'No trimmed NB file'
                	sys.exit()  
                               
                               

               	#Reload trimmed mosaics and divide
               	trimmed_red = fits.open(newRed)
               	trimmed_nb = fits.open(newNB)
               	trimmed_red_img = trimmed_red[0].data
               	trimmed_nb_img = trimmed_nb[0].data
               	
               	trimmed_nb_header = trimmed_nb[0].header
                
                                
                                
                #apply scaling factor to Halpha images
                if ccdnum in range(0,9): #A
                        scaling_factor = 2.7495
                elif ccdnum in range(9,17): #D
                        scaling_factor = 2.7265
                elif ccdnum in range(17, 25):#B
                        scaling_factor = 2.6235
                elif ccdnum in range(25, 33): #C
                        scaling_factor = 2.7195


                        
                        
                scale_array = np.full(trimmed_nb_img.shape, scaling_factor)
                trimmed_nb_img = np.multiply(trimmed_nb_img, scale_array)
                
                #Create Halpha - r image
                sub_img = np.subtract(trimmed_nb_img, trimmed_red_img)
                sub = fits.PrimaryHDU(sub_img, header=trimmed_nb_header)
                sub.writeto(sub_name, clobber=True)


 
                #close files
                trimmed_red.close()
                trimmed_nb.close()



        
