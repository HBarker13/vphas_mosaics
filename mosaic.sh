#! /bin/bash
#Script to make mosaics of vphas CCDS (must be decompressed and extracted but not necessarily confidence mapped) including a background check
#called in python script as mosaic.sh $1=workdir_path $2=colour $3=finpath $4=vphas_num


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
mAdd -p $4_$2projdir $4_$2proj_images.tbl $4_$2template.hdr $3
echo -e "Completed"

rm -rfv $4_$2workdir 
rm -rfv $4_$2projdir
rm $4_$2raw_images.tbl
rm $4_$2proj_images.tbl
rm $4_$2template.hdr
rm $4_$2stats.tbl
echo -e "Working files deleted"


