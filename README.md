# vphas_mosaics

These scripts work assuming a request for a ~complete VPHAS+ pointing has been requested from CASU (http://casu.ast.cam.ac.uk/vstsp/).
It assumes there are two u, r (red block), r (blue block), i exposures, and three g and NB (narrowband Halpha) exposures, but should work without.


These were written early in my PhD. They aren't "good" codes, but they work assuming everything is installed and set up in the PATH.




Montage must be installed: http://montage.ipac.caltech.edu/
decompress_vphas.py also uses imcopy (iraf)



Step 0.  Download a pointing from CASU. I like to save it in a directory named: vphas_####  where #### is the pointing number.
CASU: http://casu.ast.cam.ac.uk/vstsp/


If you're feeling brave, run full_mosaic.py. Below is a description of all the scripts that it calls, which can also be run sepearately.
full_mosaic.py will ask for the pointing number (as it will look for the vphas_*pointing_number* directory), if you want to flatten the background of the mosaics (a Montage feature), and if you want to bin the pixels. 




sort_vphas.py
Sorts the single (image), catalogue and calibration files in ./vphas_### into ./vphas_####_ex
Each exposure in every filter gets its own subdirectory. They havea naming format: filtername_######_blocknumber
NB. Vphas pointings are split into blocks A, B, C (g and Halpha only) to account for gaps in the ccds.
At this point, you must check the names of these subdirectories. There should be:
u_####_a
u_####_b

g_####_a
g_####_b
g_####_c

r_####_a
r_####_b

r2_###_a
r2_###_b

i_####_a
i_####_b

NB_####_a
NB_####_b
NB_####_c


Common problems: 
The ESO grade hasn't been properly input into the image header, adn the script has rejected it. (VPHAS+ has finished now, but pointings with a grade D were planned to be repeated).
The block assignment in the g and NB filters is wrong/incomplete. I would recomend opening ccd 1 from each block in rgb (in ds9 ) to see how they lie on the sky. Blocks and B should not overlap, and C should be between them.

If the block assigment is wrong, just manually change the directory name (eg. Rename g_1234_b to g_1234_c)




decompress_vphas.py
VPHAS+ data from CASU is rice compressed. I find it easier to decompress everything so I have easy access to the files. This script loops over everything and calls imcopy (from iraf).
Deccompressed files have '_decom' added to their name.




extract_ccds.py
The '_decom' files contain the data from all 32 CCDS in headers. ie. filename[1] to filename[32].
This script creates a seperate file for all 32 CCDs (in each filter). They are put in a subdirectory in the same place as the file containing the data. This made thing easier as I was familiarising myself with the data, but if you were to rewrite these scripts, this step is likely not necessary.



DO NOT USE THIS ONE-----------------------------------------------------------------------
loop_conf_overlay.py
Loops over all the image files and corresponding calibration files and removes any low confidence pixels. 
DO NOT DO THIS - it makes it harder to look for new objects.
I haven't ammended this script or the ones that follow yet, so this is still being run with the removal of pixels commented out.
------------------------------------------------------------------------------------------


div_mSubimg.py
This is a bad name from the script.
It uses a Montage commmand (mSubimage) to create Halpha / r mosaics. It uses the red block r exposure.
This can likely be improved upon. At the minute, each Halpha and r CCD it trimmed (a few rows of pixels taken off each edge) because they don't line up perfectly. If you open a Halpha and r CCD in rgb in ds9 and zoom in on an edge you'll see what I mean.


(optional)
bin_pixles.py
Uses Montage command mShrink to bin pixels. The idea was that this might make it easier to see faint, nebulous objects. 


loop_(bgcheck)_mosaic.py
Loops over all the filters and creates mosaics using Montage. Calls mosaic.sh.  Files to be included in the mosaic are moved into a working directory. See the Montage website for details.

The background flattening version of the script is a littel more complicated. Again, see the website for full details. Two .sh files are needed because if not all the CCDs are corrected, they won't be moved into the next working directory and so in the end there's holes in the mosaic.

There is a similar version of this script (mosaic.py) where you choose a filter to process.
Moaics are put in a directory vphas_####_fin / filtername_fin / filename.fits





In the end, you'll have directories that look something like:


top_dir -  vphas_0001 - single
		      - catalogues
		      - calib
		      
	- vphas_0001_ex - g_a  - single - single_ccds
			       - calib - calib_ccds
			       - catalogues
			       - confcorr   (confidence corrected files: currently commented out)
			- g_b
			- g_c
			etc ...
	- vphas_0001_fin - g_bg_fin  (background flattened)
			 - g_fin 
			




























