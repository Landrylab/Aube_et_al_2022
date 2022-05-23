REM Command lines for the simulations to make sure that a high bias towards translational mutations cannot rescue the cost-precision model.
REM Third part

cd /d F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Final_simulations_v2\Bp_higher_WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm020_iter1 --couples 2000 --seed 22 --mut_ratio 0.20 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm020_iter2 --couples 2000 --seed 33 --mut_ratio 0.20 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm020_iter3 --couples 2000 --seed 44 --mut_ratio 0.20 --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm016_iter1 --couples 2000 --seed 22 --mut_ratio 0.16 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm016_iter2 --couples 2000 --seed 33 --mut_ratio 0.16 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm016_iter3 --couples 2000 --seed 44 --mut_ratio 0.16 --dupli_type WGD

F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm01_iter1 --couples 2000 --seed 22 --mut_ratio 0.1 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm01_iter2 --couples 2000 --seed 33 --mut_ratio 0.1 --dupli_type WGD
F:\Simon_Aube\Hausser_2019_reanalysis\Simulations_ADS_Hausser\pythonProject\Genome_script.py --name Bm01_iter3 --couples 2000 --seed 44 --mut_ratio 0.1 --dupli_type WGD