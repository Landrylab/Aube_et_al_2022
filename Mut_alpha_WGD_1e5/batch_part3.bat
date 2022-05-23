REM Command lines to run the scanning of alpha parameter (asymmetry of the skew normal distributions of mutational effects) with an SD of mutational effects of 10 % (0.1)
REM For this simulation, independent mutations affecting either transcription or translation rate are considered. Final run started on 2022-01-17.
REM Third part of the script

cd /d F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Final_simulations_v2\Mut_alpha_WGD_1e5

mkdir alpha_0075_neg
cd alpha_0075_neg

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --skew -0.075 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

cd ..

mkdir alpha_025_neg
cd alpha_025_neg

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --skew -0.25 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

cd ..

mkdir alpha_035_neg
cd alpha_035_neg

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --skew -0.35 --bm_dist "(0, 0.1)" --bp_dist "(0, 0.1)" --dupli_type WGD --pop 1.0e5

cd ..