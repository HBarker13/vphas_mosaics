#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python
#functions to make lists of files/dirs in vphas directories


import glob
import os
import argparse
from astropy.io import fits
from astropy import wcs
import sys
import numpy as np
import math



#------------------------------------------------------------------------------------------------------#

#get args.vphas_num from command line
#use as: args = make_lists.get_vphas_num()
def get_vphas_num():
    parser = argparse.ArgumentParser(description="Collect vphas mosaicing and processing inputs")
    parser.add_argument('-v','--vphas_num', help="vphas pointing number", required=True)
    return parser.parse_args()
    
#------------------------------------------------------------------------------------------------------#    

#get all args from command line
def get_args():
    parser = argparse.ArgumentParser(description="Collect vphas mosaicing and processing inputs")
    parser.add_argument('-v','--vphas_num', help="vphas pointing number", required=True)
    parser.add_argument('-b','--bin_choice', help="bin pixels? y/n", required=True)
    parser.add_argument('-s','--bin_level', help="bin shrink factor (optional)", required=False)
    return parser.parse_args()

#------------------------------------------------------------------------------------------------------#


#add a column to a recarray
def append_table(table, name, arr, d_type):
    arr = np.asarray(arr)
    dtype = d_type
    newdtype = np.dtype(table.dtype.descr + [(name, d_type)])
    newtable = np.empty(table.shape, dtype=newdtype)
    for field in table.dtype.fields:
        newtable[field] = table[field]
    newtable[name] = arr
    return newtable



#------------------------------------------------------------------------------------------------------#


#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):	
	dtype_list = ['>f4' for item in title_list]
	str_list = ['vphas_num', 'png', 'pn_name','name', 'v_name', 'PN', 'data_release', 'sourceID', 'primaryID', 'warning_u', 'detectionID_u', 'warning_g', 'detectionID_g', 'warning_r2', 'detectionID_r2', 'warning_ha', 'detectionID_ha', 'warning_r', 'detectionID_r', 'warning_i', 'detectionID_i', 'field', 'spectype', 'Companion_SpecType_(i)', 'Lower_SpecType_(i)', 'Upper_SpecType_(i)', 'Companion SpecType (J)', 'Lower SpecType (J)', 'Upper SpecType (J)', 'Abundance', 'filtername', 'Filter_1', 'Filter_2', 'Filter_3', 'Filter_4', 'Filter_5', 'Filter_6', 'Fname_1', 'Fname_2', 'Fname_3', 'Fname_4', 'Fname_5', 'Fname_6', 'pn', 'block', 'Lum_class']
	for ind, val in enumerate(title_list):
		if val in str_list:
			dtype_list[ind]='|S20'
			
	if str_override==True:
		dtype_list = ['|S20' for item in title_list]
			
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]

	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]
		data_array.append(col)

	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array	


#-------------------------------------------------------------------------------------------------------#
#creates lists of the a,b and c pointings in one vphas directory. Returns in order u,g,r_r,r_b,i,NB
#NOTE: other scrpts rely on the paths being returned in this order
def bandmerge_list(ex_path):

	a_u = None
	a_g = None
	a_r = None
	a_r2 = None
	a_i = None
	a_NB = None
	
	b_u = None
	b_g = None
	b_r = None
	b_r2 = None
	b_i = None
	b_NB = None
	
	c_g = None
	c_NB = None



	#list of all the u, g, r, i, and NB sorted directories
	dirs = [dirpath for dirpath in glob.glob(ex_path+'/*') if os.path.isdir(dirpath) and 'trimmed' not in dirpath and 'Halpha' not in dirpath]
	
	for dirpath in dirs:
		_, dirname = dirpath.rsplit('/', 1)
		
		#dirname has the form: filtername_date_block
		filtername, _, _, block = dirname.split('_')
		
		
		if block=='a': 
			if filtername == 'u':
				a_u = dirpath
			if filtername == 'g':
				a_g = dirpath 
			if filtername == 'r':
				a_r = dirpath 
			if filtername == 'r2':
				a_r2 = dirpath 
			if filtername == 'i':
				a_i = dirpath 
			if filtername == 'NB':
				a_NB = dirpath 

		elif block=='b':
			if filtername == 'u':
				b_u = dirpath
			if filtername == 'g':
				b_g = dirpath 
			if filtername == 'r':
				b_r = dirpath 
			if filtername == 'r2':
				b_r2 = dirpath 
			if filtername == 'i':
				b_i = dirpath 
			if filtername == 'NB':
				b_NB = dirpath 


		elif block=='c':
			if filtername =='g':
				c_g = dirpath
			if filtername =='NB':
				c_NB = dirpath
		
			
	a_block = [a_u, a_g, a_r, a_r2, a_i, a_NB] 
	b_block = [b_u, b_g, b_r, b_r2, b_i, b_NB]
	c_block = [c_g, c_NB]		
			
			
	for val in a_block:
		if val==None:
			print 'a block is incomplete'
			for line in a_block:
				print line
			raw_input('Press any key to continue')
			
	for val in b_block:
		if val==None:
			print 'b block is incomplete'
			for line in b_block:
				print line
			raw_input('Press any key to continue')

	for val in c_block:
		if val==None:
			print 'c block is incomplete'
			for line in c_block:
				print line
			raw_input('Press any key to continue')

	
			
	return a_block, b_block, c_block


        
