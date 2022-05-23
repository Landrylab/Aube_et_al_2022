REM Script to assess whether the postulated post-duplication expression optimum affects the reproduction of the empirical expression divergence patterns
REM Part 1

cd /d F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Final_simulations_v2\Dupli_effect

mkdir effect_16
cd effect_16

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --dupli_type WGD --dupli_effect 1.6
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --dupli_type WGD --dupli_effect 1.6
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --dupli_type WGD --dupli_effect 1.6

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --dupli_type WGD --dupli_effect 1.6
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --dupli_type WGD --dupli_effect 1.6
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --dupli_type WGD --dupli_effect 1.6

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --dupli_type WGD --dupli_effect 1.6
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --dupli_type WGD --dupli_effect 1.6
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --dupli_type WGD --dupli_effect 1.6

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --dupli_type WGD --dupli_effect 1.6
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --dupli_type WGD --dupli_effect 1.6
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --dupli_type WGD --dupli_effect 1.6

cd ..

mkdir effect_186
cd effect_186

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter1 --couples 2000 --seed 22 --mut_ratio 1 --dupli_type WGD --dupli_effect 1.86
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter2 --couples 2000 --seed 33 --mut_ratio 1 --dupli_type WGD --dupli_effect 1.86
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm1_iter3 --couples 2000 --seed 44 --mut_ratio 1 --dupli_type WGD --dupli_effect 1.86

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter1 --couples 2000 --seed 22 --mut_ratio 3 --dupli_type WGD --dupli_effect 1.86
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter2 --couples 2000 --seed 33 --mut_ratio 3 --dupli_type WGD --dupli_effect 1.86
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm3_iter3 --couples 2000 --seed 44 --mut_ratio 3 --dupli_type WGD --dupli_effect 1.86

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter1 --couples 2000 --seed 22 --mut_ratio 4 --dupli_type WGD --dupli_effect 1.86
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter2 --couples 2000 --seed 33 --mut_ratio 4 --dupli_type WGD --dupli_effect 1.86
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm4_iter3 --couples 2000 --seed 44 --mut_ratio 4 --dupli_type WGD --dupli_effect 1.86

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter1 --couples 2000 --seed 22 --mut_ratio 8 --dupli_type WGD --dupli_effect 1.86
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter2 --couples 2000 --seed 33 --mut_ratio 8 --dupli_type WGD --dupli_effect 1.86
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm8_iter3 --couples 2000 --seed 44 --mut_ratio 8 --dupli_type WGD --dupli_effect 1.86

cd ..