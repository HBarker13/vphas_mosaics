#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python
#Make (unbackground checked) mosaics from confidence corrected ccds of the same colour

import subprocess as sp
import os
import glob
import shutil
from datetime import datetime
import sys

import make_lists


args = make_lists.get_args()


startTime = datetime.now()
print "Mosaicing start time: %s" %startTime

ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)

filternames = ['u', 'g', 'r', 'r2', 'i', 'NB', 'Halpha_div_R']
filtername = None
while filtername not in filternames:
	filtername = raw_input('Choose a filter: (u, g, r, r2, i, NB, Halpha_div_R) ' )




#Make finDir
fin_dir_path = os.getcwd() + '/vphas_' + args.vphas_num + '_fin'
if not os.path.exists(fin_dir_path):
	os.makedirs(fin_dir_path)



#Make workdir
workdir_path = os.getcwd() + '/'+args.vphas_num + '_' + filtername +'_workdir'
if not os.path.exists(workdir_path):
	os.makedirs(workdir_path)
	print "Working directory created"


   
dirnames = glob.glob( ex_path + '/'+filtername+'*' )
if len(dirnames)==0:
 	print 'No directories could be found'
   	print filtername
   	print ex_path
   	raw_input('Paused')
   
 

#bin choice = y
if args.bin_choice == 'y':
	

      	fin_dir = fin_dir_path+'/'+filtername +'bin'+str(args.bin_level)+'_fin'
  	fin_fpath = fin_dir + '/' + filtername +'bin'+str(args.bin_level)+ 'mosaic.fits'
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
		sys.exit()



	#copy all the files into the working directory
	for fpath in all_bin_files:
       		_, fname = fpath.rsplit('/', 1)
       		if not os.path.exists(workdir_path+'/'+fname):
       			shutil.copy(fpath, workdir_path)
	print "Files copied to working directory"
      		
	#call the mosaicking script
	sp.call(["mosaic.sh", workdir_path, filtername, finFile_path, args.vphas_num])  





#bin choice = 'n'
else:
   	
      	fin_dir = fin_dir_path+'/'+ filtername +'_fin'
	fin_fpath = fin_dir + '/' + filtername +'mosaic.fits'
	if not os.path.exists(fin_dir):
      		os.makedirs(fin_dir)
      		print "FinDir created"


	if os.path.exists( fin_fpath ):
      		print "Mosaic already exists"
      		sys.exit()
         		
       	all_fpaths = []	
       	print dirnames
       	for dirpath in dirnames:
       		print dirpath
        	
        	
       		if filtername == 'Halpha_div_R':
       			print 'HERE'
       			file_paths = glob.glob( dirpath + '/no_bin/*/*.fit')
      
     
       		else:
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

     		

	#$1=workdir_path $2=colour choice $3=finPath $4=vphas_num
      	sp.call(["mosaic.sh", workdir_path, filtername, fin_fpath, args.vphas_num])
     




finTime = datetime.now() - startTime
print "Mosaicking finished: %s" %finTime





#remove the working directory
if os.path.exists(workdir_path):
	os.remove(workdir_path)


























print '---------------Complete---------------'
print datetime.now() - startTime





