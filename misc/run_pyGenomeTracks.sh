#!/bin/bash
# scp wanghm@210.27.94.249:~/whm/mango/pb_cov/gz23_coverage.bw ./
# get HGT gene type
#cut -f 1,13 /Users/alexwang/0data/0mango/accessory/HGT/gz23_putative_HGT.blastp.tbl > hgt.list
# get color mapped bed file
#awk 'BEGIN{FS=OFS="\t"}NR==FNR{cc[$1]=$2}NR>FNR && $3~/gene/{split($9,id,"\"");if(id[2] in cc){if(cc[id[2]] == "Bacteria"){col="0,0,255"}else{col="255,0,0"}}else{col="190,190,190"};print $1,$4,$5,id[2],0,$7,".",".",col}' hgt.list ../../accessory/gz23_accessory.gtf > gz23_accessory.bed
# initiate track.ini
make_tracks_file --trackFiles gz23_coverage.bw gz23_accessory.bed -o tracks.ini

# plot
pyGenomeTracks --tracks tracks.ini --region gz23_12:1-601234 --width 120 --outFileName gz23_12_image.pdf
pyGenomeTracks --tracks tracks.ini --region gz23_13:1-594284 --width 120 --outFileName gz23_13_image.pdf
