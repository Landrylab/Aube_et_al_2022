REM Script to repeat all the simulations with independant mutational effects using a wide range of standard deviations for the distributions of mutational effects.
REM Second step of the scanning of standard deviations

cd /d F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Final_simulations_v2\Mut_SD_WGD

mkdir 2_5_percent
cd 2_5_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.025)" --bp_dist "(0, 0.025)" --dupli_type WGD

cd ..

mkdir 5_percent
cd 5_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.05)" --bp_dist "(0, 0.05)" --dupli_type WGD

cd ..

mkdir 7_5_percent
cd 7_5_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.075)" --bp_dist "(0, 0.075)" --dupli_type WGD

cd ..

mkdir 10_percent
cd 10_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD

cd ..