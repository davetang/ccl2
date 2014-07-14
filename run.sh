#!/bin/bash

#download files
wget -O CNhs14214.bam http://fantom.gsc.riken.jp/5/datafiles/phase1.3/basic/human.timecourse.hCAGE/hIPS%252c%2520biol_rep1.CNhs14214.14380-156B6.hg19.nobarcode.bam
wget -O CNhs14215.bam http://fantom.gsc.riken.jp/5/datafiles/phase1.3/basic/human.timecourse.hCAGE/hIPS%252c%2520biol_rep2.CNhs14215.14381-156B7.hg19.nobarcode.bam
wget -O CNhs14216.bam http://fantom.gsc.riken.jp/5/datafiles/phase1.3/basic/human.timecourse.hCAGE/hIPS%252c%2520biol_rep3.CNhs14216.14382-156B8.hg19.nobarcode.bam
wget -O CNhs14217.bam http://fantom.gsc.riken.jp/5/datafiles/phase1.3/basic/human.timecourse.hCAGE/hIPS%2520%252bCCl2%252c%2520biol_rep1.CNhs14217.14383-156B9.hg19.nobarcode.bam
wget -O CNhs14218.bam http://fantom.gsc.riken.jp/5/datafiles/phase1.3/basic/human.timecourse.hCAGE/hIPS%2520%252bCCl2%252c%2520biol_rep2.CNhs14218.14384-156C1.hg19.nobarcode.bam
wget -O CNhs14219.bam http://fantom.gsc.riken.jp/5/datafiles/phase1.3/basic/human.timecourse.hCAGE/hIPS%2520%252bCCl2%252c%2520biol_rep3.CNhs14219.14385-156C2.hg19.nobarcode.bam

#filter out reads not passing quality control
#and filter out reads with a mapping quality of 9 or less
for file in `ls *.bam`;
   do echo $file
   base=`basename $file .bam`
   samtools view -F 512 -q 10 -b $file > ${base}_F512.bam
done
