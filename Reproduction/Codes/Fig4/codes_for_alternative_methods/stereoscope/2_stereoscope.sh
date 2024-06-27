### 240329 new sc data
ls ./*sp_count.tsv | while read id
do
  file=`basename $id`
  sample=${file%%_sp*}
  stereoscope run --sc_cnt ./240329_sc_count.tsv --sc_labels ./240329_meta.tsv \
  -sce 25000 -o ./stereo_res/240329_new_scdata/$sample -n 5000 --st_cnt $id \
  -ste 25000 -stb 500 -scb 500 -gp
done
