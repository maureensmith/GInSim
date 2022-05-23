#!/bin/bash
SIMU_DIR="/path/to/scripts"

RESULT_DIR="/path/to/save/results"

#sinusoidal growth rate
p_repl=(1.05)
p_repl_2=(0.95)
##### shredding the reads simulating sewage data
#probability for each position to be a cut position
p_shred=0.005
#probability for each segment to be dismissed (degraded)
p_dismiss=0.3

p_mut=(0.0001) #original!

L=1000
N=50 #original

T_FINAL=120
N_SIM=5

num_intro=(0) #5 10 50 100 500)


for pm in ${p_mut[@]}
do
	echo $pm
	for pr in ${p_repl[@]}
	do
		echo $pr
		for pr2 in ${p_repl_2[@]}
		do
			echo $pr2
				for i in $(seq $N_SIM)
				do
					 #i=3
					echo i: $i
          #for sw in ${T_SWITCH_ORIG[@]}
          #do
              #echo $pr
			  	SUB_RES_DIR="$RESULT_DIR"_pmut_"$pm"
					SUB_RES_DIR_PROPERROR="$SUB_RES_DIR"
					CONFIG_DIR="$SUB_RES_DIR"/config
					mkdir -p  $SUB_RES_DIR
					mkdir -p  $SUB_RES_DIR_PROPERROR


					printf "\n\n***********************************************\n\n"
					echo "# Running simulation..."
					printf "\n\n***********************************************\n\n"

					# conda activate py39 -> required version!
					python $SIMU_DIR/generate_simulation_data.py \
						-file_prefix sim_$i \
						-L $L \
						-N_init $N \
						-p_mut $pm \
						-t_final $T_FINAL \
						-output $SUB_RES_DIR \
						-p_rep $pr \
						-p_rep2 $pr2 \
						-p_shred $p_shred \
						-p_dismiss $p_dismiss

					printf "\n\n...finished simulation $i \n\n"
				done #end i
		done # end pr2
	done # end pr
done # end pm
