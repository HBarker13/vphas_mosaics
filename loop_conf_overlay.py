#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python


from astropy.io import fits
import numpy as np
import sys, os
import shutil
import glob
import subprocess as sp

import make_lists

args = make_lists.get_vphas_num()
current_path  = os.getcwd()
ex_path = current_path + '/vphas_' + args.vphas_num + '_ex'
paths_tree =  make_lists.list_dirs(ex_path)
available_colours = make_lists.list_colours(ex_path)

print "Available colours: "
print available_colours
for colour_choice in available_colours:
    print "Colour choice: %s " %colour_choice
    colour_choice += '_'
    
    branches_to_process = make_lists.colour_branches_wanted(paths_tree, colour_choice)

#Make lists of paths to ccds for single and conf
    for branch in branches_to_process:
        print "Processing %s " %branch[0]  
        for leaf in branch:
            if 'single' in leaf and 'ccds' in leaf:
                imdirPath = leaf
                im_ccdPath_list = glob.glob(imdirPath + '/*')
            if 'conf' in leaf and 'ccds' in leaf:
                    confdirPath = leaf
                    conf_ccdPath_list = glob.glob(confdirPath + '/*')

        finDir = branch[0] + '/confcorr'
        if os.path.exists(finDir):
            print "Confidence corrected directory already exists"
            continue
        if not os.path.exists(finDir):
            os.makedirs(finDir)
            print "Created final directory: %s" %finDir

        text_fpath = branch[0] + '/ccd_confs.txt'
        textfile = open(text_fpath, 'w')
        confsum = 0.0
        textfile.write("File                         Confidence\n")

        #match images to conffiles
        for i in range(1,33):
            for impath in im_ccdPath_list:
                if 'ccd'+str(i)+'.fit' in impath:
                    singleFile = impath
                    single = fits.open(singleFile)
                    confName = single[1].header['CIR_CPM']
                    for confPath in conf_ccdPath_list:
                        _,fname = confPath.rsplit('/', 1)
                        if i<10:
                            if fname == confName[:-7]+'_ccd'+str(i)+'.fit':
                                confFile = confPath
                        elif i>=10:
                            if fname == confName[:-8]+'_ccd'+str(i)+'.fit':
                                confFile = confPath
                    conf = fits.open(confFile)
                    conf_array = conf[1].data
                    
                    #set confidence level and remove pixels below
                    stdev = np.std(conf_array)
                    print "Standard deviation: %f" %stdev
                    median = 100 #conf maps are normalised to median 100
                    #confLevel = 95  #limit used in Drew2014 vphas intro paper
                    confLevel = int(median-(3*stdev)) #use 3 * stdev. Rounds down
                    print "Confidence cut-off: %i " %confLevel
                    im_array = single[1].data
                    newData = np.where(conf_array>confLevel, im_array, float('NaN'))
                    #remove negative image pixels
                    newData = np.where(newData>=0, newData, float('NaN'))
                    
                    #save confidence corrected image
                    _, finName = branch[0].rsplit('/', 1)
                    finName = finName + '_ccd' + str(i) + '_confcorr.fit'
                    finPath = finDir + '/' + finName
                
                    newHeader = single[1].header
                    newHeader.append(('COMMENT', 'Low confidence pixels (median-3standard deviations) set to nan'))
                    confcorr = fits.PrimaryHDU(newData, newHeader)
                    confcorr.writeto(finPath, clobber=True)
                    
                    _, singleName = singleFile.rsplit('/', 1)
                    text_line = singleName + '         ' + str(confLevel) +'\n'
                    textfile.write(text_line)
                    confsum += confLevel

                    single.close()
                    conf.close()
        mean_conf = confsum/32.0
        textfile.write('\n')
        textline = 'Mean confidence: ' + str(mean_conf)
        textfile.write(textline)
        textfile.close()

