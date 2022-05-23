REM Command lines to run the scanning of the correlation coefficient of mutational effects with a reference SD of 12.5 % (0.125)
REM First half of the script.

cd /d F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Final_simulations_v2\Mut_correlations_WGD_1e5

mkdir corr_0
cd corr_0

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bivariate True --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

cd ..

mkdir corr_01
cd corr_01

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bivariate True --correlation 0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

cd ..

mkdir corr_01_neg
cd corr_01_neg

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bivariate True --correlation -0.1  --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bivariate True --correlation -0.1 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

cd ..

mkdir corr_02
cd corr_02

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bivariate True --correlation 0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

cd ..

mkdir corr_02_neg
cd corr_02_neg

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33  --mut_ratio 5 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
 
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bivariate True --correlation -0.2 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

cd ..

mkdir corr_03
cd corr_03

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter1 --couples 2000 --seed 22 --mut_ratio 2 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter2 --couples 2000 --seed 33 --mut_ratio 2 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm2_iter3 --couples 2000 --seed 44 --mut_ratio 2 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter1 --couples 2000 --seed 22 --mut_ratio 5 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter2 --couples 2000 --seed 33 --mut_ratio 5 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm5_iter3 --couples 2000 --seed 44 --mut_ratio 5 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter1 --couples 2000 --seed 22 --mut_ratio 6 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter2 --couples 2000 --seed 33 --mut_ratio 6 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm6_iter3 --couples 2000 --seed 44 --mut_ratio 6 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter1 --couples 2000 --seed 22 --mut_ratio 7 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter2 --couples 2000 --seed 33 --mut_ratio 7 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm7_iter3 --couples 2000 --seed 44 --mut_ratio 7 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter1 --couples 2000 --seed 22 --mut_ratio 9 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter2 --couples 2000 --seed 33 --mut_ratio 9 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm9_iter3 --couples 2000 --seed 44 --mut_ratio 9 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter1 --couples 2000 --seed 22 --mut_ratio 10 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter2 --couples 2000 --seed 33 --mut_ratio 10 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm10_iter3 --couples 2000 --seed 44 --mut_ratio 10 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --bivariate True --correlation 0.3 --bm_dist "(0, 0.125)" --bp_dist "(0, 0.125)" --dupli_type WGD --pop 1.0e5

cd ..