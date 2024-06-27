##### 20240109 regeneration data
@ /data/zhangff/ST/data/regeneration
ls ./*sp_count.tsv | while read id
do
  file=`basename $id`
  sample=${file%%_sp*}
  stereoscope run --sc_cnt ./regeneration_sc_count.tsv --sc_labels ./regeneration_meta.tsv \
  -sce 25000 -o ./stereo_res/$sample -n 5000 --st_cnt $id \
  -ste 25000 -stb 500 -scb 500 -gp
done
