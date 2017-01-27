#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python
#Make (unbackground checked) mosaics from confidence corrected ccds of the same colour

import subprocess as sp
import os
import glob
import shutil
from datetime import datetime

import make_lists


startTime = datetime.now()
current_path  = os.getcwd()
args = make_lists.get_args()
ex_path = current_path + '/vphas_' + args.vphas_num + '_ex'

available_colours = make_lists.list_colours_full(ex_path, args.bin_choice, args.bin_level)
   
colour_choice = ''
while colour_choice not in available_colours:
   print "Available confidence corrected colours: ",
   for colour in available_colours:
      print colour,
   print
   colour_choice = raw_input("Enter the colour of mosaic you want to create: ")
colour_choice += '_'


colour_list = make_lists.chosen_colour_pathlist(ex_path, colour_choice)

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

   #Make workdir
   workdir_path = current_path+'/'+colour_choice+'bin'+args.bin_level+'_workdir'
   if not os.path.exists(workdir_path):
      os.makedirs(workdir_path)
      print "Working directory created"
        
   #Make finDir(s)
   root_fin_path = current_path +'/vphas_' + args.vphas_num + '_fin'
   if not os.path.exists(root_fin_path):
      os.makedirs(root_fin_path)

   finDir_path = root_fin_path+'/'+colour_choice+'bin'+args.bin_level+'_finDir'
   finFile_path = finDir_path + '/' + colour_choice+'bin'+args.bin_level+ 'mosaic.fits'
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

   #$1=workdir_path $2=colour choice $3=finaldir/fin_filename $4=vphas_num 
   sp.call(["mosaic.sh", workdir_path, colour_choice[:-1], finFile_path, args.vphas_num])   

   #delete workdir
   shutil.rmtree(workdir_path)

elif args.bin_choice == 'n':
   #Make workdir
   workdir_path = current_path + '/' + colour_choice + 'workdir'
   if not os.path.exists(workdir_path):
      os.makedirs(workdir_path)
      print "Working directory created"

   #Make finDir(s)
   root_fin_path = current_path +'/vphas_'+ args.vphas_num + '_fin'
   if not os.path.exists(root_fin_path):
      os.makedirs(root_fin_path)

   finDir_path = root_fin_path  + '/' + colour_choice + 'finDir'
   finFile_path = finDir_path + '/' + colour_choice + 'mosaic.fits'
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


   sp.call(["mosaic.sh", workdir_path, colour_choice[:-1], finFile_path, args.vphas_num])
  

   #delete workdir
   shutil.rmtree(workdir_path)

print datetime.now() - startTime
