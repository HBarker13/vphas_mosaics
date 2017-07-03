#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python


from astropy.io import fits
import numpy as np
import sys, os
import shutil
import glob
import subprocess as sp

import make_lists

args = make_lists.get_vphas_num()
ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'



#loop over directories of each filter
filternames = ['u', 'g', 'r', 'r2', 'i', 'NB']
#filternames = ['NB']
for filtername in filternames:

	print 'Filter', filtername
	dirnames = glob.glob( ex_path + '/' + filtername + '*')
	
	for dirname in dirnames:
	
		block = dirname[-1]
	
		single_ccds = glob.glob( dirname + '/single/*ccds/*.fit')
		
		conf_ccds = glob.glob( dirname + '/calib/*ccds/*.fit')
		
		
		#create a directory to put the condifence corrected images in
		confcorr_path = dirname + '/confcorr'
		if not os.path.exists( confcorr_path):
			os.makedirs(confcorr_path)
			


        	#match the image and confidence files
        	for ccdnum in range(1,33):
        	
        	
        		#skip if the confidence corrected file already exists
			confcorr_fpath = confcorr_path + '/'+block+'_ccd'+str(ccdnum)+'.fits'
			#if os.path.exists(confcorr_fpath):
			#	continue
        	
        		single_file = [ f for f in single_ccds if '_ccd'+str(ccdnum)+'.fit' in f]
        		if len(single_file)!=1:
        			print "Can't find single file"
        			print single_file
        			raw_input('Paused')
        		single_file = single_file[0]
        			
        	
        		single = fits.open(single_file)
        		conf_name = single[1].header['CIR_CPM']
        		conf_name, _ = conf_name.split('.fit', 1) #strip the .fit[1]

        		
        		conf_ccd = dirname + '/calib/'+conf_name+'_ccds/'+conf_name+'_ccd'+str(ccdnum)+'.fit'
   			if not os.path.exists( conf_ccd):
        			print "Can't find conf file"
        			print single_file
        			print conf_ccd
        			print conf_name
        			raw_input('Paused')
        			
        		#open the confidence file
        		conf = fits.open(conf_ccd)
        		conf_array = conf[1].data
        		
        		
        
                    	#set confidence level and remove pixels below
                    	stdev = np.std(conf_array)
                    
                    	median = 100 #conf maps are normalised to median 100
                    	#confLevel = 95  #limit used in Drew2014 vphas intro paper
                    	confLevel = int(median-(3*stdev)) #use 3 * stdev. Rounds down
                    	
                    	#print 'CCD', ccdnum
                    	#print "Standard deviation: %f" %stdev
                    	#print "Confidence cut-off: %i " %confLevel
                    	
                    
                    
                    	#remove image pixels where the confidence is below the chosen limit
                    	im_array = single[1].data
                    	#newData = np.where(conf_array>confLevel, im_array, float('NaN'))
                    	#remove negative image pixels
                    	#newData = np.where(newData>=0, newData, float('NaN'))
                    
                    
                    	newData = im_array
                    
                    	#save confidence corrected image                
                    	newHeader = single[1].header
                    	newHeader.append(('COMMENT', 'Low confidence pixels (median-3standard deviations) set to nan'))
                    	confcorr = fits.PrimaryHDU(newData, newHeader)
                    	confcorr.writeto( confcorr_fpath, clobber=True)

			single.close()
			conf.close()
        	
        	
        	
        	
        	
   
