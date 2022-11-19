#!/bin/bash
gfortran linear_algebra_cholesky.f90 -o linear_algebra_cholesky -fopenmp
echo "Enter n_max"
read n_max
echo "Enter threads"
read max_num_threads
echo "Enter amount of experiments"
read experiments
for (( experiment=1; experiment<=$experiments; experiment++ ))
do
	echo "Experiment #$experiment"
	for (( num_threads=1; num_threads<=$max_num_threads; num_threads++ ))
	do
		export OMP_NUM_THREADS=$num_threads
		echo "num_threads=$num_threads"
		filename=../data/$experiment.$num_threads.txt
		echo "# n; Execution time, sec" > $filename
		for n in `seq 100 50 $n_max`
		do
			echo "$n 0" | ./linear_algebra_cholesky >> $filename
		done
	done
done