#PBS -N EXMIRAS-submit
#PBS -A NEOL0014
#PBS -q main
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=1GB

# loop over runIdeal.sh and change the temperature in intervals of 5

sleep 28800
for T0 in $(seq 285 5 310)
# for T0 in 290 310
do
    for HUMIDITY in $(seq 0.1 0.1 0.9)
    do
        echo "Running for T0 = $T0, HUMIDITY = $HUMIDITY"
        # Create a job script for each temperature
        cat << EOF > runEXMIRAS_$T0-$HUMIDITY.sh
#PBS -N EMRS-$T0-$HUMIDITY
#PBS -A NEOL0014
#PBS -q main
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=32:mem=64GB

cd /glade/u/home/nbarron/workshop/EXMIRAS/test
export T0=$T0
export HUMIDITY=$HUMIDITY
matlab -nosplash -nodisplay -r "run 'runIdeal.m'; exit"
EOF

    # Submit the job
    qsub runEXMIRAS_$T0-$HUMIDITY.sh

    # Optional: wait for a few seconds before submitting the next job
    sleep 1
    done
done