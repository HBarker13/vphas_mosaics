#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#Complete set of scripts needed to download data and process and make un-background flattened mosaics

import subprocess as sp
import os
import getopt
import sys
from datetime import datetime


#request number is the number supplied by the CASU website
reqNo = raw_input('Enter request number: ')

#the pointing number is used to label the directories
vphas_num = raw_input('Enter the vphas pointing number: ')

#use Montage to even out the background level of the images
bg_choice = ''
while bg_choice != 'y' and bg_choice != 'n':
    bg_choice = raw_input("Do you want to blackground flatten the mosaics (y/n)? ")

bin_choice = ''
while bin_choice != 'y' and bin_choice != 'n':
    bin_choice = raw_input('Do you want to bin ccd pixels? (y/n) ') 
if bin_choice == 'y':
    bin_level = raw_input("Choose binning level: ")



startTime = datetime.now()


#A scipt I use to download requested files from CASU. This isn't uploaded as it includes my password
print "Downloading..."
sp.call(["download_vphas.sh", reqNo, vphas_num])
print
print "Sorting files..."
os.system("sort_vphas.py -v %s" %vphas_num)
print 'Please check the sorted files'
raw_input('Press any key to continue')
print
print "Decompressing files..."
os.system("decompress_vphas.py -v %s" %vphas_num)
print
print "Extracting ccds..."
os.system("extract_ccds.py -v %s" %vphas_num)	
print 
print "Confidence overlaying..."
os.system("loop_conf_overlay.py -v %s" %vphas_num)
print
print "Making Halpha divide by r mosaics..."
os.system("div_mSubimg.py -v %s" %vphas_num)
print


if bin_choice == 'y':
    print "Binning pixels..."
    os.system("bin_pixels.py -v %s -b %s -s %s" %(vphas_num, bin_choice, bin_level))
    print
    print "Mosaicking binned pixels..."
    os.system("loop_mosaic.py -v %s -b %s -s %s" %(vphas_num, bin_choice, bin_level))
 
     
print "Mosaicking full resolution..."
if bg_choice == 'n':
    os.system("loop_mosaic.py -v %s -b %s" %(vphas_num, 'n'))
elif bg_choice == 'y':
    os.system("loop_bgcheck_mosaic.py -v %s -b %s" %(vphas_num, 'n'))

finTime = datetime.now() - startTime
print "Elapsed time: %s" %finTime

