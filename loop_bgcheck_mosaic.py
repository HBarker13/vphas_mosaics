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



args = make_lists.get_args()



startTime = datetime.now()


print "Start time: %s" %startTime

ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'
#all_colours = ['u', 'g', 'r_r', 'r_b', 'i', 'NB', 'Halpha_div_R']
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)
all_colours = ['Halpha_div_R']


#Make finDir
fin_dir_path = os.getcwd() + '/vphas_' + args.vphas_num + '_fin'
if not os.path.exists(fin_dir_path):
	os.makedirs(fin_dir_path)



filternames = ['u', 'g', 'r', 'r2', 'i', 'NB', 'Halpha_div_R']

#loop through colours, taking one directory from each block
for filtername in filternames:
	print "Filter: ", filtername
	
	
	
	#Make workdir
	workdir_path = os.getcwd() + '/'+args.vphas_num + '_bg_' + filtername +'workdir'
	if not os.path.exists(workdir_path):
		os.makedirs(workdir_path)
		print "Working directory created"
	

	dirnames = glob.glob( ex_path + '/'+filtername+'*' )
	#break if there is only one directory: the background smoothing needs
	#at least 2 directories to work
	if len(dirnames)<2 and filtername!='Halpha_div_R':
		print 'Not enough directories found'
		continue
	   	
	   	
   

	#bin choice = y
	if args.bin_choice == 'y':
	
	
	      	fin_dir = fin_dir_path+'/'+filtername +'bin'+str(args.bin_level)+'_bg_fin'
      		fin_fpath = fin_dir + '/' + filtername +'bin'+str(args.bin_level)+ '_bg_mosaic.fits'
      		if not os.path.exists(fin_dir):
      		  	os.makedirs(fin_dir)
      		  	print "FinDir created"	
	
		all_bin_files = []	
	
		for dirpath in dirnames:
	
			if filtername == 'Halpha_div_R_':
				bin_files = glob.glob( dirpath + '/bin'+args.bin_level + '/*/*.fit')
				

	        	else:
	            		bin_files = glob.glob( dirpath + '/bin'+ str(args.bin_level) + '/*.fit')
	            		
	            	
	        	for fpath in bin_files:
	         		all_bin_files.append(fpath)	
	            		

		if len(all_bin_files)==0:
			print 'No binned files found'
			print filtername
			print dirpath
			continue



		#copy all the files into the working directory
		for fpath in all_bin_files:
            		_, fname = fpath.rsplit('/', 1)
            		if not os.path.exists(workdir_path+'/'+fname):
               			shutil.copy(fpath, workdir_path)
      		print "Files copied to working directory"
      		

	      	#$1=workdir_path $2 = colour choice $3=finaldir/fin_filename $4=vphas_num
	      	sp.call(["bgcheck_mosaic_clipped.sh", workdir_path, filtername+'bg_', fin_fpath, args.vphas_num])


	      	projPath = os.getcwd() + '/' + args.vphas_num + '_' + filtername +'bg_' + 'projdir'
	      	corrPath = os.getcwd() + '/' + args.vphas_num + '_' + filtername +'bg_' + 'corrdir'
           
		projFiles = os.listdir(projPath)
      		corrFiles = os.listdir(corrPath)


		#copy files from the proj directory to the corr directory if they've been missed
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








	#if bin_choice = 'n'
	else:
	
	
		fin_dir = fin_dir_path+'/'+filtername +'_bg_fin'
      		fin_fpath = fin_dir + '/' + filtername + '_bg_mosaic.fits'
      		if not os.path.exists(fin_dir):
        		os.makedirs(fin_dir)
        		print "FinDir created"

     
		#break if the mosaic already exists      
		if os.path.exists(fin_fpath):
			print "Mosaic already exists"
			continue
			
			
			

         	all_fpaths = []	
         	for dirpath in dirnames:
         	
         		if filtername == 'Halpha_div_R':
         			file_paths = glob.glob( dirpath + '/no_bin/*/*.fit')
      
     
         		else:
            			file_paths = glob.glob( dirpath + '/confcorr/*.fits')
            			
            		for f in file_paths:
            			all_fpaths.append(f)
            			
            	if len(all_fpaths)==0:
            		print 'No confcorr files found'
            		print dirpath
            		continue

         
         
		#copy all the files into the working directory
		for fpath in all_fpaths:
            		_, fname = fpath.rsplit('/', 1)
            		if not os.path.exists(workdir_path+'/'+fname):
               			shutil.copy(fpath, workdir_path)
      		print "Files copied to working directory"
      		
      		
      		

      		#$1=workdir_path $2 = colour choice $3=finaldir/fin_filename $4 = vphas_num
      		sp.call(["bgcheck_mosaic_clipped.sh", workdir_path, filtername+'bg_', fin_fpath, args.vphas_num])


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
































