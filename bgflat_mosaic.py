#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#Make (background checked) mosaics from confidence corrected ccds of the same filter using my data tree

#can the bin = y/n be made more concise?

import subprocess as sp
import os
import glob
import shutil
from datetime import datetime
import sys
import argparse

import make_lists


# input form:  bgflat_mosaic.py -v (pointing_number) -b (y/n) -s (#)
#where -b is the choice whether to use binned images as input (yes or no) and -s is the shrink factor of the binned images.
# The -s flag is not needed if binned data isn't used
args = make_lists.get_args()



startTime = datetime.now()
print "Start time: %s" %startTime

ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)



filternames = ['u', 'g', 'r', 'r2', 'i', 'NB', 'Halpha_div_r']
filtername = None
while filtername not in filternames:
	filtername = raw_input('Choose a filter: (u, g, r, r2, i, NB, Halpha_div_r) ' )



#Make finDir to save mosaics to
fin_dir_path = os.getcwd() + '/vphas_' + args.vphas_num + '_fin'
if not os.path.exists(fin_dir_path):
	os.makedirs(fin_dir_path)


#Make working directory
workdir_path = os.getcwd() + '/'+args.vphas_num + '_bgflat_' + filtername +'workdir'
if not os.path.exists(workdir_path):
	os.makedirs(workdir_path)
	print "Working directory created"
	

#find the directories containing files to mosaic 
dirnames = glob.glob( ex_path + '/'+filtername+'*' )

#break if there is only one directory: the background smoothing needs
#at least 2 directories to work
if len(dirnames)<2 and filtername!='Halpha_div_r':
	print 'Not enough directories found'
   	print filtername
   	print ex_path
   	sys.exit()
	   	
	   	
   

#bin choice = y ie. want to mosaic binned pixels
if args.bin_choice == 'y':
	
	
      	fin_dir = fin_dir_path+'/'+filtername +'bin'+str(args.bin_level)+'_bg_fin'
 	fin_fpath = fin_dir + '/' + filtername +'bin'+str(args.bin_level)+ '_bg_mosaic.fits'
   	if not os.path.exists(fin_dir):
      	  	os.makedirs(fin_dir)
      	  	print "FinDir created"	
	
	all_bin_files = []	
	
	for dirpath in dirnames:
	
		if filtername == 'Halpha_div_r':
			bin_files = glob.glob( dirpath + '/bin'+args.bin_level + '/*/*.fit')
				

	       	else:
          		bin_files = glob.glob( dirpath + '/bin'+ str(args.bin_level) + '/*.fit')
	            		
	            	
        	for fpath in bin_files:
         		all_bin_files.append(fpath)	
	            		
	            		
	            		
	#end if no binned imgs are found
	if len(all_bin_files)==0:
		print 'No binned files found'
		print filtername
		print dirpath
		sys.exit()



	#copy all the files into the working directory
	for fpath in all_bin_files:
        	_, fname = fpath.rsplit('/', 1)
            	if not os.path.exists(workdir_path+'/'+fname):
               		shutil.copy(fpath, workdir_path)
      	print "Files copied to working directory"
      		

	#call the mosaicking bash script
	#$1=workdir_path $2 = colour choice $3=finaldir/fin_filename $4=vphas_num
	sp.call(["bgcheck_mosaic_clipped.sh", workdir_path, filtername+'bg_', fin_fpath, args.vphas_num])


	#copy files from the proj directory to the corr directory if they've been missed.
	#Sometimes CCDs don't need to be projected? Or somehow just don't get processed
	#the mosaics come out fine if we just move them ourselves
	projPath = os.getcwd() + '/' + args.vphas_num + '_' + filtername +'bg_' + 'projdir'
	corrPath = os.getcwd() + '/' + args.vphas_num + '_' + filtername +'bg_' + 'corrdir'
           
	projFiles = os.listdir(projPath)
     	corrFiles = os.listdir(corrPath)



	copiedFiles = []
	for pName in projFiles:
        	if pName not in corrFiles:
        	    copiedFiles.append(pName)
        	    pPath = projPath + '/' + pName
        	    cPath = corrPath + '/' + pName
        	    shutil.copy(pPath, cPath)
        
      
     	#print "Files copied from projdir to corrdir: "
      	#for line in copiedFiles:
      		#   print line
	
	#final mosaicking
      	sp.call(["bg_Add.sh", fin_fpath, filtername+'bg_', args.vphas_num])








#if bin_choice = 'n'  ie. use full resolution imgs
else:
	
	
	fin_dir = fin_dir_path+'/'+filtername +'_bg_fin'
     	fin_fpath = fin_dir + '/' + filtername + '_bg_mosaic.fits'
      	if not os.path.exists(fin_dir):
        	os.makedirs(fin_dir)
        	print "FinDir created"

     
	#break if the mosaic already exists      
	if os.path.exists(fin_fpath):
		print "Mosaic already exists"
		sys.exit()
			
			
			

         all_fpaths = []	
         for dirpath in dirnames:
         	
         	if filtername == 'Halpha_div_r':
         		file_paths = glob.glob( dirpath + '/no_bin/*/*.fit')
      
     
         	else:
         		#confidence corrected files :creates 'holes' in the image
         		#that can make things more difficult rather that easier         		

         		#file_paths = glob.glob( dirpath + '/confcorr/*.fits')	        
     			file_paths = glob.glob( dirpath + '/single/*_ccds/*.fit')
            			
            		for f in file_paths:
            			all_fpaths.append(f)
            			
            if len(all_fpaths)==0:
      		print 'No input CCD files found'
            	print dirpath
            	sys.exit()

         
         
	#copy all the files into the working directory
	for fpath in all_fpaths:
          	_, fname = fpath.rsplit('/', 1)
            	if not os.path.exists(workdir_path+'/'+fname):
               		shutil.copy(fpath, workdir_path)
      	print "Files copied to working directory"
      		
      		
      		
	#call the mosaicking bash script
      	#$1=workdir_path $2 = colour choice $3=finaldir/fin_filename $4 = vphas_num
      	sp.call(["bgcheck_mosaic_clipped.sh", workdir_path, filtername+'bg_', fin_fpath, args.vphas_num])
      	

	#copy files from the proj directory to the corr directory if they've been missed.
	#Sometimes CCDs don't need to be projected? Or somehow just don't get processed
	#the mosaics come out fine if we just move them ourselves
	projPath = os.getcwd() + '/' + args.vphas_num + '_' + filtername +'bg_' + 'projdir'
      	corrPath = os.getcwd() + '/' + args.vphas_num + '_' + filtername +'bg_' + 'corrdir'
         
      
      	#copy files that were not corrected
      	projFiles = os.listdir(projPath)
      	corrFiles = os.listdir(corrPath)

      	copiedFiles = []
      	for pName in projFiles:
        	if pName not in corrFiles:
           		copiedFiles.append(pName)
            		pPath = projPath + '/' + pName
            		cPath = corrPath + '/' + pName
            		shutil.copy(pPath, cPath)
        
      	#print "Files copied from projdir to corrdir: "
      	#for line in copiedFiles:
      	#   print line
      		
      		
	#final mosaic
      	sp.call(["bg_Add.sh", fin_fpath, filtername+'bg_', args.vphas_num])



finTime = datetime.now() - startTime
print "Finished: %s" %finTime
































