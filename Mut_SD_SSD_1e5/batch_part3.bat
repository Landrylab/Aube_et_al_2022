REM Script to repeat all the simulations with independant mutational effects using a wide range of standard deviations for the distributions of mutational effects.
REM This time with a reduced efficacy of selection
REM Third step of the scanning of standard deviations

cd /d F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Final_simulations_v2\Mut_SD_SSD_1e5

mkdir 12_5_percent
cd 12_5_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type SSD --pop 1.0e5

cd ..

mkdir 15_percent
cd 15_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.15)" --bp_dist "(0, 0.15)" --dupli_type SSD --pop 1.0e5

cd ..

mkdir 20_percent
cd 20_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)"  --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.2)" --bp_dist "(0, 0.2)" --dupli_type SSD --pop 1.0e5

cd ..

mkdir 25_percent
cd 25_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.25)" --bp_dist "(0, 0.25)" --dupli_type SSD --pop 1.0e5

cd ..

mkdir 35_percent
cd 35_percent

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bm_dist "(0, 0.35)" --bp_dist "(0, 0.35)" --dupli_type SSD --pop 1.0e5

cd ..