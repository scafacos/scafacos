# @ job_name = pp3mg_fcs_fortran
# @ comment = "BGL Job by Size"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL;
# @ wall_clock_limit = 12:00:00
# @ notification = error 
# @ notify_user = l.westphal@fz-juelich.de
# @ job_type = bluegene
# @ bg_size = 4096
# @ queue 

dir_exe=./

for p in 4096 2048 1024 512 256 128 64 32
do

        /usr/local/bin/mpirun -exe ${dir_exe}/test_pp3mg > `/bin/pwd`/pp3mg_p${p}.log \
        -mode VN -np $p -verbose 1 -args "../inp_data/pp3mg/conf_1000000.in"

done
