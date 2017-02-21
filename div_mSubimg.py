#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#delete pixels in red and nb ccds that don't overlap and save new cdds. Uses mSubimage from Montage
#red and nb directories containing ccds of the same pointing must first be paired. Pairing is attemped using ra/dec in crval comments, then using crpix headers, then manually.  The list of paired directories is then saved.

from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import glob
import shutil
import sys
import subprocess as sp
import math

import make_lists

args = make_lists.get_vphas_num()
current_path = os.getcwd()
ex_path = current_path + '/vphas_' + args.vphas_num + '_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)

#list pairs of nb and (red block) r dirs
pair_list = [[a_block[2], a_block[5]], [b_block[2], b_block[5]]]
      

for pair in pair_list:
        print pair
        red_ccds = glob.glob(pair[0]+'/confcorr/*.fit')
        nb_ccds = glob.glob(pair[1]+'/confcorr/*.fit')
        
        _, redName = pair[0].rsplit('/', 1)
        _, nbName = pair[1].rsplit('/', 1)
                
        red_finDir = ex_path+'/trimmed_red/R' + redName
        nb_finDir = ex_path +'/trimmed_nb/NB' + nbName
                
        if not os.path.exists(red_finDir):
                os.makedirs(red_finDir)
                print "Created final directory: %s" %red_finDir

        if not os.path.exists(nb_finDir):
                os.makedirs(nb_finDir)
                print "Created final directory: %s" %nb_finDir

        div_finDir = ex_path+'/Halpha_div_R/no_bin/'+nbName+redName
        if not os.path.exists(div_finDir):
                os.makedirs(div_finDir)
                print "Created final directory: %s" %div_finDir

        #sub_finDir = ex_path+'/Halpha_sub_R/no_bin/'+nbName+redName
        #if not os.path.exists(sub_finDir):
        #        os.makedirs(sub_finDir)
        #        print "Created final directory: %s" %sub_finDir


        for i in range(1,33):
                for fpath in red_ccds:
                        if '_ccd'+str(i)+'_' in fpath:
                                redpath = fpath
                for fpath in nb_ccds:
                        if '_ccd'+str(i)+'_' in fpath:
                                nbpath = fpath                 
                
                newRed = red_finDir+'/'+redName+'_ccd'+str(i)+'.fit'
                #break if ccd already exists
                #if os.path.exists(newRed):
                #        print "Trimmed red ccd already exists"
                #        continue

                newNB = nb_finDir+'/'+nbName+'_ccd'+str(i)+'.fit'
                #break if ccd already exists
                #if os.path.exists(newNB):
                #        print "Trimmed nb ccd already exists"
                #        continue

                divName = div_finDir+'/'+nbName+redName+'_ccd'+str(i)+'.fit'
                #break if ccd already exists
                if os.path.exists(divName):
                        print "Halpha/R ccd already exists"
                        continue
                        
                #subName = sub_finDir+'/NB'+nbName+'_R'+redName+'_ccd'+str(i)+'.fit'
                #break if ccd already exists
                #if os.path.exists(divName):
                #        print "Halpha/R ccd already exists"
                #        continue

        
                nb = fits.open(nbpath)
                red = fits.open(redpath)
                red_header = red[0].header
                nb_header = nb[0].header
                red_wcs = wcs.WCS(red_header)
                nb_wcs = wcs.WCS(nb_header)
                red_img = red[0].data
                nb_img = nb[0].data
                                
                
                #check for nb wcs shifts relative to red in raw images and shift nb ccd into wcs of the red
                cx = red_img.shape[0]/2
                cy = red_img.shape[1]/2
                cra, cdec = red_wcs.wcs_pix2world(cx, cy, 1) #red as reference
                refx, refy = nb_wcs.wcs_world2pix(cra, cdec, 1) #pixels in nb with same wcs as cra and cdec in red

                #shiftx = change of centre col, shifty = change of centre row
                shiftx = refx-cx
                if shiftx <1 :
                        shiftx = int(math.floor(shiftx))
                else:
                        shiftx = int(math.ceil(shiftx))
                        
                shifty = refy-cy
                if shifty <1:
                        shifty = int(math.floor(shifty))
                else:
                        shifty = int(math.ceil(shifty))
                        
                if shiftx<0:
                        red_x_start=str(-shiftx) #1 counting?
                        red_x_size = str(red_img.shape[1]+shiftx)
                        nb_x_start = '1'
                        nb_x_size =str(nb_img.shape[1]+shiftx)
                elif shiftx>=0:
                        red_x_start='1'
                        red_x_size = str(red_img.shape[1]-shiftx)
                        nb_x_start = str(shiftx)
                        nb_x_size =str(nb_img.shape[1]-shiftx)
                        

                if shifty>=0:
                        red_y_start = '1'
                        red_y_size= str(red_img.shape[0]-shifty)
                        nb_y_start = str(shifty)
                        nb_y_size = str(nb_img.shape[0]-shifty)  
                elif shifty<0:
                        red_y_start = str(-shifty)
                        red_y_size= str(red_img.shape[0]+shifty)
                        nb_y_start = '1'
                        nb_y_size = str(nb_img.shape[0]+shifty)                       
                        
                print newRed
                print newNB
                print red_img.shape
                print shiftx, shifty
                        
                print red_x_start, red_x_size
                print nb_x_start, nb_x_size
                print red_y_start, red_y_size
                print nb_y_start, nb_y_size
                
                sp.call(["mSubimage", "-p", redpath, newRed, red_x_start, red_y_start, red_x_size, red_y_size])
                sp.call(["mSubimage", "-p", nbpath, newNB, nb_x_start, nb_y_start, nb_x_size, nb_y_size])
                                

                #Reload trimmed mosaics and divide
                trimmed_red = fits.open(newRed)
                trimmed_nb = fits.open(newNB)
                trimmed_red_img = trimmed_red[0].data
                trimmed_nb_img = trimmed_nb[0].data
                trimmed_red_header = red[0].header

                #apply scaling factor to Halpha images
                if i in range(0,9): #A
                        scaling_factor = 2.923
                elif i in range(9,17): #D
                        scaling_factor = 2.877
                elif i in range(17, 25):#B
                        scaling_factor = 2.71
                elif i in range(25, 33): #C
                        scaling_factor =2.755
                scale_array = np.full(trimmed_nb_img.shape, scaling_factor)
                trimmed_nb_img = np.multiply(trimmed_nb_img, scale_array)
                
                #divide images
                divided_img = np.zeros((trimmed_red_img.shape[0], trimmed_red_img.shape[1]))
                divided_img = np.true_divide(trimmed_nb_img, trimmed_red_img)
                divided = fits.PrimaryHDU(divided_img, header=trimmed_red_header)
                divided.writeto(divName, clobber=True)

		"""
                #subracted images
                subtracted_img = np.zeros((trimmed_red_img.shape[0], trimmed_red_img.shape[1]))
                subtracted_img = np.subtract(trimmed_nb_img, trimmed_red_img)
                subtracted = fits.PrimaryHDU(subtracted_img, header=trimmed_red_header)
                subtracted.writeto(subName, clobber=True)
                """
 
                #close files
                red.close()
                nb.close()
                trimmed_red.close()
                trimmed_nb.close()



        
