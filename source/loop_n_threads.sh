#!/bin/bash
# compiling the program
gfortran linear_algebra_cholesky.f90 -o linear_algebra_cholesky -fopenmp
# Entering the maximum number of threads and the amount of experiments
echo "Enter maximum number of threads"
read max_num_threads
echo "Enter amount of experiments"
read experiments
# Running the program in the following loop
for (( experiment=1; experiment<=$experiments; experiment++ ))
do
   echo "Experiment #$experiment"
   for (( num_threads=1; num_threads<=$max_num_threads; num_threads++ ))
   do
      export OMP_NUM_THREADS=$num_threads
      echo "num_threads=$num_threads"
      filename=../data/1/$experiment.$num_threads.txt
      echo "# n; Execution time, sec" > $filename
      for n in `seq 100 1100 2500`
      do
         echo "n=$n"
         # keeping time data
         echo "$n" | ./linear_algebra_cholesky >> $filename
      done
   done
done
