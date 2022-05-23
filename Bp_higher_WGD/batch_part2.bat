REM Command lines for the simulations to make sure that a high bias towards translational mutations cannot rescue the cost-precision model.
REM Second part

cd /d F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Final_simulations_v2\Bp_higher_WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter1 --couples 2000 --seed 22 --mut_ratio 0.5 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter2 --couples 2000 --seed 33 --mut_ratio 0.5 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm05_iter3 --couples 2000 --seed 44 --mut_ratio 0.5 --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm033_iter1 --couples 2000 --seed 22 --mut_ratio 0.33 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm033_iter2 --couples 2000 --seed 33 --mut_ratio 0.33 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm033_iter3 --couples 2000 --seed 44 --mut_ratio 0.33 --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm025_iter1 --couples 2000 --seed 22 --mut_ratio 0.25 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm025_iter2 --couples 2000 --seed 33 --mut_ratio 0.25 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm025_iter3 --couples 2000 --seed 44 --mut_ratio 0.25 --dupli_type WGD