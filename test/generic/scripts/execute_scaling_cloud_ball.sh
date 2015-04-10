#!/bin/bash -x

taskslist="1 2 4"

for tasks in $taskslist; do
  formattasks=$(printf "%02d" $tasks)
  mpirun -np $tasks ./scafacos_test p2nfft halver_ref.xml.gz -c tolerance_field,1e-6 > id00-p${formattasks}.out
done



