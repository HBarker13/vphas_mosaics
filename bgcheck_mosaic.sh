#! /bin/bash
#Script to make mosaics of vphas CCDS (must be decompressed and extracted but not necessarily confidence mapped) including a background check
#$1=workdir_path $2=colour $3=fin_filepath $4=vphas_num

mImgtbl $1 $4_$2raw_images.tbl 
echo -e "Raw image table made"
mMakeHdr $4_$2raw_images.tbl $4_$2template.hdr
echo -e "Header file made"
#Make projdir
mkdir -p $4_$2projdir
mProjExec -p $1 $4_$2raw_images.tbl $4_$2template.hdr $4_$2projdir $4_$2stats.tbl
echo -e "Images projected"
mImgtbl $4_$2projdir $4_$2proj_images.tbl
echo -e "New image table made"
mOverlaps $4_$2proj_images.tbl $4_$2diffs.tbl
echo -e "Overlaps determined"
#Make diffdir
mkdir -p $4_$2diffdir
mDiffExec -p $4_$2projdir $4_$2diffs.tbl $4_$2template.hdr $4_$2diffdir
echo -e "Difference images created"
mFitExec $4_$2diffs.tbl $4_$2fits.tbl $4_$2diffdir
echo -e "Plane-fitting coeffs calculated"
mBgModel $4_$2proj_images.tbl $4_$2fits.tbl $4_$2corrections.tbl
echo -e "Corrections table created"
#Make corrdir
mkdir -p $4_$2corrdir
mBgExec -p $4_$2projdir $4_$2proj_images.tbl $4_$2corrections.tbl $4_$2corrdir
echo -e "Background matched, reprojected images made"
mAdd -p corrdir $4_$2proj_images.tbl $4_$2template.hdr $3
echo -e "Mosaic made"


rm -rfv $4_$2projdir
rm -rfv $4_$2corrdir
rm -rfv $4_$2diffdir
rm $4_$2raw_images.tbl
rm $4_$2template.hdr
rm $4_$2stats.tbl
rm $4_$2diffs.tbl
rm $4_$2fits.tbl
rm $4_$2proj_images.tbl
rm $4_$2corrections.tbl
echo -e "Mosaiking directories deleted"


