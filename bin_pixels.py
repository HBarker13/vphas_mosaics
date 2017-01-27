#!/usr/local/anaconda/bin/python

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

current_path  = os.getcwd()
ex_path = current_path + '/vphas_' + args.vphas_num + '_ex'

#need to look for unbinned Halpha_div_R and Halpha_sub_R
available_colours = make_lists.list_colours_full(ex_path, 'n', '')

for colour in available_colours:
   colour_choice = colour
   colour_choice += '_'

   colour_list = make_lists.chosen_colour_pathlist(ex_path, colour_choice)

   #for each dir containing image ccds of selected colour, get paths to all confidence corrected ccds
   for dirpath in colour_list:
      if colour_choice == 'Halpha_div_R_' or colour_choice == 'Halpha_sub_R_':
         sub_dirpaths = glob.glob(dirpath+'/no_bin/*')
         for sub_path in sub_dirpaths:
            _, sub_dir_name = sub_path.rsplit('/', 1)
            bin_dirPath = dirpath +'/bin'+args.bin_level+'/'+sub_dir_name
            if not os.path.exists(bin_dirPath):
               os.makedirs(bin_dirPath)
            
            for fpath in glob.glob(sub_path+'/*.fit'):
               print fpath
               _, ccdName = fpath.rsplit('/',1)
               finPath = bin_dirPath + '/' + ccdName[:-4] + '_bin' + str(args.bin_level) + '.fit'
               if os.path.exists(finPath):
                  print "Binned file already exists"
                  continue
                       
               #1=in.fits, 2=out.fits 3=factor
               sp.call(["mShrink", fpath, finPath, args.bin_level])
            
      else: #if colour not Halpha_div_R
         bin_dirPath = dirpath + '/bin' + str(args.bin_level)  
         if not os.path.exists(bin_dirPath):
            os.makedirs(bin_dirPath)
            print "%s created" %bin_dirPath

         for fpath in glob.glob(dirpath +'/confcorr/*.fit'):
            print fpath
            _, ccdName = fpath.rsplit('/',1)
            finPath = bin_dirPath + '/' + ccdName[:-4] + '_bin' + str(args.bin_level) + '.fit'
            if os.path.exists(finPath):
               print "Binned file already exists"
               continue
        
            #1=in.fits, 2=out.fits 3=factor
            sp.call(["mShrink", fpath, finPath, args.bin_level])

print "Binning completed"



