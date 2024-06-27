#!/bin/bash
cd /data/huangzr/Spatial/Kras/version4/repeats/ST_KRAS/scTE
#clean data
for i in `ls file to bam file`
do
  a=${i:46}
  b=${a%%/*}
	samtools view ${i} -h | awk '/^@/ || /CB:/' | samtools view -h -b > ./${b}.clean.bam
done
#build index
(for i in `ls file to cleaned bam file`
do
	samtools index -b ${i} &
done; wait)
#mapping TE
for i in `ls file to cleaned bam file`
do
  a=${i##*/}
  scTE -i $i -o ${a//.clean.bam/} -x /data/huangzr/Spatial/Kras/version5/Example_data/hg38.exclusive.idx --hdf5 True -CB CB -UMI UB
done
