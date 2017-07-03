#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#script to search the for ccds of the same colour and bin the pixels in each ccd

import os
import glob
from astropy.io import fits
import numpy as np
import math
import subprocess as sp
import sys

import make_lists

args = make_lists.get_args()
ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'

filternames = ['u', 'g', 'r', 'r2', 'i', 'NB', 'Halpha_div_R']

for filtername in filternames:

	dirnames = glob.glob( ex_path + '/'+filtername+'*' )
	
	#loop over the block a,b (and c) directories
	for dirpath in dirnames:
	
		print dirpath
		
		#create a directory to store the binned pixels in
		bin_dir = dirpath + '/bin'+str(args.bin_level)
		if not os.path.exists(bin_dir):
			os.makedirs(bin_dir)
			
	
			
		if filtername!='Halpha_div_R':
		
			#loop over the confidence corrected ccds
			ccds = glob.glob( dirpath + '/confcorr/*.fits')
			if len(ccds)==0:
				print 'No ccds found'
				print dirpath
				raw_input('Paused')
			
			for ccd in ccds:
				
				_, ccdname = ccd.rsplit('/', 1)
				ccdname = ccdname[:-5] #remove '.fits
				
				newname = bin_dir + '/' + ccdname + '_bin' + str(args.bin_level)+'.fits'
				if os.path.exists(newname):
					print 'Binned file already exists'
					continue
				
				
				#call the bash scripts, that uses Montage, to bin the pixels
				#1=in.fits, 2=out.fits 3=bin factor
         			sp.call(["mShrink", ccd, newname, args.bin_level])
         			
         			
         			
         			
         			
         			
         	#if filtername == Halpha_div_R		
         	else:
         		#find the a and b block directories
			block_dirs = glob.glob( dirpath + '/no_bin/*')

			
			for block in block_dirs:
			
				#make a subdirectory to keep the binned files in
				_,block_name = 	block.rsplit('/', 1)
				block_name += '_bin'+str(args.bin_level)
				
				binned_block_dir = bin_dir + '/' + block_name
				if not os.path.exists(binned_block_dir):
					os.makedirs(binned_block_dir)
					
				
				#get all the unbinned ccds
				ccds = glob.glob( block + '/*.fit')
				if len(ccds)==0:
					print 'No ccds found'
					print block
					raw_input('Paused')

					
				for ccd in ccds:
				
					_, ccdname = ccd.rsplit('/', 1)
					ccdname = ccdname[:-4] #remove '.fit.
					
					newname = binned_block_dir + '/' + ccdname + '_bin' + str(args.bin_level)+'.fits'
					#if os.path.exists(newname):
				#	print 'Binned file already exists'
				#	continue
				
				
				
					#call the bash scripts, that uses Montage, to bin the pixels
					#1=in.fits, 2=out.fits 3=factor
              				sp.call(["mShrink", ccd, newname, args.bin_level])
					
		
		
		
		
			
print 'Binning completed'		










