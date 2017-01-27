#!/usr/local/anaconda/bin/python

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
current_path  = os.getcwd()
ex_path = current_path + '/vphas_' + args.vphas_num + '_ex'
all_colours = ['u', 'g', 'r_r', 'r_b', 'i', 'NB', 'Halpha_div_R']
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)

#loop through colours, taking one directory from each block
for colour_choice in all_colours:
   print "Colour: %s" %colour_choice
   if colour_choice=='u': colour_list = [a_block[0], b_block[0]]
   if colour_choice=='g': colour_list = [a_block[1], b_block[1], c_block[1]]
   if colour_choice=='r_r': colour_list = [a_block[2], b_block[2]]
   if colour_choice=='r_b': colour_list = [a_block[3], b_block[3]]
   if colour_choice=='i': colour_list = [a_block[4], b_block[4]] 
   if colour_choice=='NB': colour_list = [a_block[5], b_block[5], c_block[0]]
   if colour_choice=='Halpha_div_R': colour_list = [ex_path + '/Halpha_div_R']
   colour_choice = colour_choice+'_'
   
   #colour_list = make_lists.chosen_colour_pathlist(ex_path, colour_choice)
 
   if args.bin_choice == 'y':
      bin_subdirs = []
      for dirpath in colour_list:
         if colour_choice == 'Halpha_div_R_' or colour_choice == 'Halpha_sub_R_':
            parent_dir = dirpath + '/bin'+args.bin_level
            child_dirs = glob.glob(parent_dir+'/*')
            for bin_dir in child_dirs:
               bin_subdirs.append(bin_dir)
         else:
            bin_dir = dirpath + '/bin'+args.bin_level
            if os.path.exists(bin_dir):
               bin_subdirs.append(bin_dir)
  
      all_file_paths = []
      for binpath in bin_subdirs:
         file_paths = glob.glob(binpath + '/*.fit')
         all_file_paths.append(file_paths)

      if len(all_file_paths)<2:       #break if less than one directory of ccds (not enough for background smoothing)
         print "Not enough ccds for background matching"
         incomplete_colours.append(colour_choice[:-1])
         continue

      #Make workdirs
      workdir_path = current_path + '/'+ args.vphas_num + '_' + colour_choice +'bg_workdir'
      if not os.path.exists(workdir_path):
         os.makedirs(workdir_path)
         print "Working directory created"

      #Make finDir(s)
      root_fin_path = current_path + '/vphas_' + args.vphas_num + '_fin'
      if not os.path.exists(root_fin_path):
         os.makedirs(root_fin_path)

      finDir_path = root_fin_path+'/'+colour_choice+'bin'+args.bin_level+'bg_finDir'
      finFile_path = finDir_path+'/'+ colour_choice+'bin'+args.bin_level+'bg_mosaic.fits'
      if not os.path.exists(finDir_path):
         os.makedirs(finDir_path)
         print "FinDir created"
         
      for binpath in bin_subdirs:
         file_paths = glob.glob(binpath + '/*.fit')
         for fpath in file_paths:
            _, fname = fpath.rsplit('/', 1)
            if not os.path.exists(workdir_path+'/'+fname):
               shutil.copy(fpath, workdir_path)
      print "Files copied to working directory"
         
      #$1=workdir_path $2 = colour choice $3=finaldir/fin_filename $4=vphas_num
      sp.call(["bgcheck_mosaic_clipped.sh", workdir_path, colour_choice+'bg_', finFile_path, args.vphas_num])

      projPath = current_path + '/' + args.vphas_num + '_' + colour_choice+'bg_' + 'projdir'
      corrPath = current_path + '/' + args.vphas_num + '_' + colour_choice+'bg_' + 'corrdir'
            
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

      sp.call(["bg_Add.sh", finFile_path, colour_choice+'bg_', args.vphas_num])


   if args.bin_choice == 'n':
      #break if <33 ccds
      #all_file_paths = []
      #for rootdir in colour_list:
      #   if colour_choice == 'Halpha_div_R_'or colour_choice == 'Halpha_sub_R_':
      #      parent_dirs = glob.glob(rootdir + '/no_bin/*')
      #      for parent in parent_dirs:
      #         file_paths = glob.glob(parent+ '/*.fit')
       #        all_file_paths.append(file_paths)

        # else:
         #   file_paths = glob.glob(rootdir + '/confcorr/*.fit')
          #  all_file_paths.append(file_paths)

      #if len(all_file_paths)<2:
      #   print "Not enough ccds for background matching"
      #   incomplete_colours.append(colour_choice[:-1])
      #   continue

      #Make workdir
      workdir_path = current_path + '/'+ args.vphas_num + '_' + colour_choice +'bg_workdir'
      if not os.path.exists(workdir_path):
         os.makedirs(workdir_path)
         print "Working directory created"

      #Make finDir(s)
      root_fin_path = current_path + '/vphas_' + args.vphas_num + '_fin'
      if not os.path.exists(root_fin_path):
         os.makedirs(root_fin_path)
               
      finDir_path = root_fin_path + '/' + colour_choice  + 'bg_finDir'
      finFile_path = finDir_path + '/' + colour_choice + 'bg_mosaic.fits'
      
      
      if os.path.exists(finFile_path):
          continue
      
      if not os.path.exists(finDir_path):
         os.makedirs(finDir_path)
         print "FinDir created"

      for rootdir in colour_list:
         if colour_choice == 'Halpha_div_R_' or colour_choice == 'Halpha_sub_R_':
            file_paths = glob.glob(rootdir + '/no_bin/*/*.fit')
         else:
            file_paths = glob.glob(rootdir + '/confcorr/*.fit')
         for fpath in file_paths:
            _, fname = fpath.rsplit('/', 1)
            if not os.path.exists(fpath+'/'+fname):
               shutil.copy(fpath, workdir_path)
      print "Files copied to working directory"

      #$1=workdir_path $2 = colour choice $3=finaldir/fin_filename $4 = vphas_num
      sp.call(["bgcheck_mosaic_clipped.sh", workdir_path, colour_choice+'bg_', finFile_path, args.vphas_num])

      projPath = current_path + '/' + args.vphas_num + '_' + colour_choice+'bg_' + 'projdir'
      corrPath = current_path + '/' + args.vphas_num + '_' + colour_choice+'bg_' + 'corrdir'
         
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

      sp.call(["bg_Add.sh", finFile_path, colour_choice+'bg_', args.vphas_num])

#print "Incomplete colours: "
#print incomplete_colours

finTime = datetime.now() - startTime
print "Finished: %s" %finTime
