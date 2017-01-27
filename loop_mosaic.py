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

current_path  = os.getcwd()
ex_path = current_path + '/vphas_' + args.vphas_num + '_ex'
all_colours = ['u', 'g', 'r_r', 'r_b', 'i', 'NB', 'Halpha_div_r']
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
   if colour_choice=='Halpha_div_r': colour_list = ex_path + '/Halpha_div_R'
   colour_choice = colour_choice+'_'
   
   #colour_list = make_lists.chosen_colour_pathlist(ex_path, colour_choice)
 
   #bin choice
   if args.bin_choice == 'y':
      bin_subdirs = [] #list to all dirs contaning ccds of right colour and bin level
      for dirpath in colour_list:
         if colour_choice == 'Halpha_div_R_' or colour_choice == 'Halpha_sub_R_':
            parent_dir = dirpath + '/bin'+args.bin_level
            child_dirs = glob.glob(parent_dir+'/*')
            for bin_dir in child_dirs:
               bin_subdirs.append(bin_dir)

         else:
            bin_dir = dirpath + '/bin'+ str(args.bin_level)
            if os.path.exists(bin_dir):
               bin_subdirs.append(bin_dir)

    #Make workdir
      workdir_path = current_path + '/'+args.vphas_num + '_' + colour_choice +'workdir'
      if not os.path.exists(workdir_path):
         os.makedirs(workdir_path)
         print "Working directory created"

    #Make finDir(s)
      root_fin_path = current_path + '/vphas_' + args.vphas_num + '_fin'
      if not os.path.exists(root_fin_path):
         os.makedirs(root_fin_path)

      finDir_path = root_fin_path+'/'+colour_choice+'bin'+str(args.bin_level)+'_finDir'
      finFile_path = finDir_path + '/' + colour_choice+'bin'+str(args.bin_level)+ 'mosaic.fits'
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
      sp.call(["mosaic.sh", workdir_path, colour_choice, finFile_path, args.vphas_num])  


#bin choice = 'n'
   else:
    #Make workdir
      workdir_path = current_path + '/'+ args.vphas_num + '_' + colour_choice +'workdir'
      if not os.path.exists(workdir_path):
         os.makedirs(workdir_path)
         print "Working directory created"

    #Make finDir(s)
      root_fin_path = current_path + '/vphas_' + args.vphas_num + '_fin'
      if not os.path.exists(root_fin_path):
         os.makedirs(root_fin_path)

      finDir_path = root_fin_path +'/'+colour_choice  + 'finDir'
      finFile_path = finDir_path + '/' + colour_choice + 'mosaic.fits'

      if os.path.exists(finFile_path):
         print "Mosaic already exists"
         continue
      
      if not os.path.exists(finDir_path):
         os.makedirs(finDir_path)
         print "FinDir created"
      

      for rootdir in colour_list:
         if colour_choice == 'Halpha_div_R_' or colour_choice == 'Halpha_sub_R_':
            file_paths = glob.glob(rootdir +'/no_bin/*/*.fit')
             
         else:
            file_paths = glob.glob(rootdir + '/confcorr/*.fit')

         for fpath in file_paths:
            _, fname = fpath.rsplit('/', 1)
            if not os.path.exists(fpath+'/'+fname):
               shutil.copy(fpath, workdir_path)
      print "Files copied to working directory"

#$1=workdir_path $2=colour choice $3=finPath $4=vphas_num
      sp.call(["mosaic.sh", workdir_path, colour_choice, finFile_path, args.vphas_num])
     

finTime = datetime.now() - startTime
print "Mosaicing finished: %s" %finTime

