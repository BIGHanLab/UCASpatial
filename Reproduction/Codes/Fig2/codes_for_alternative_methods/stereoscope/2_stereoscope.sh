##### simu data 20231113
cd /data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Data/Fig2/stereoscope_result/
# mkdir stereoscope_res/
ls ../simu_TC/*tsv | while read id
do
  file=`basename $id`
  sample=${file%%_sp*}
  stereoscope run --sc_cnt ./sc_count.tsv --sc_labels ./meta.tsv \
  -sce 25000 -o ./stereoscope_res/simu_TC/$sample -n 5000 --st_cnt $id \
  -ste 25000 -stb 500 -scb 500 -gp
done

ls ../simu_TM/*tsv | while read id
do
  file=`basename $id`
  sample=${file%%_sp*}
  stereoscope run --sc_cnt ./sc_count.tsv --sc_labels ./meta.tsv \
  -sce 25000 -o ./stereoscope_res/simu_TM/$sample -n 5000 --st_cnt $id \
  -ste 25000 -stb 500 -scb 500 -gp
done

ls ../simu_TS/*tsv | while read id
do
  file=`basename $id`
  sample=${file%%_sp*}
  stereoscope run --sc_cnt ./sc_count.tsv --sc_labels ./meta.tsv \
  -sce 25000 -o ./stereoscope_res/simu_TS/$sample -n 5000 --st_cnt $id \
  -ste 25000 -stb 500 -scb 500 -gp
done