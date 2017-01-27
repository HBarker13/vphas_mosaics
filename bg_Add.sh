#! /bin/bash

#$1=finpath $2=colour $3=vphas_num

mAdd -p $3_$2corrdir $3_$2proj_images.tbl $3_$2template.hdr $1
echo -e "Mosaic made"

rm -rfv $3_$2projdir
rm -rfv $3_$2corrdir
rm -rfv $3_$2diffdir
rm -rfv $3_$2workdir
rm $3_$2raw_images.tbl
rm $3_$2template.hdr
rm $3_$2stats.tbl
rm $3_$2diffs.tbl
rm $3_$2fits.tbl
rm $3_$2proj_images.tbl
rm $3_$2corrections.tbl
echo -e "Mosaiking directories deleted"

