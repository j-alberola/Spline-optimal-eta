#!/usr/bin/env qpsh

#SBATCH -J bu.p
#SBATCH -p zen3 -c 96  -n 1 -N 1
#SBATCH --mem=300000


# Set name of the molecule and basis set
#INPUT='N2_minus'
#BASIS=basis_N2_aug-cc-pvtz-3s3p3d


# Source QP2
module load intel
module load lcpq/gnu/gcc-7.5.0
source ~/PHD/qp2/quantum_package.rc


#LAUNCH DESIRED CALCUALTION WITH SPECIFIED PARAMETERS
#
##MODIFIY FORTRAN SCRIPT TO WRITE OUTPUT IN A GIVEN FILE

qp set_file $1
#qp edit -s [1,2]
qp set cap do_cap true
qp set cap eta_cap 0.001
qp set cap eta_step_size 0.002
qp set cap n_steps_cap  2
qp set davidson_keywords n_det_max_full 1000
qp set davidson_keywords state_following false
qp set determinants read_wf True
#qp set determinants n_states 3

qp run diagonalize_h_cap > $1.scan.out
cp $1.scan.out $1.scan.script.1st.out
./script_data_corrected.sh $1 2 1
./Spline.x << EOF
$1.data1.good
Minima_find
EOF
threshold=0.00003
qp set cap eta_step_size 0.
qp set cap n_steps_cap  1
old_energy=$(awk 'NR == 1 {print $2}' Minima.dat)



while true; do


new_min=$(awk 'NR == 1 || $6 < min {min = $6; val = $1} END {print val}' Minima.dat)
echo $new_min > min_file
rm -f Minima.dat
qp set cap eta_cap $new_min
qp run diagonalize_h_cap > $1.scan.out
cp $1.scan.out $1.scan.script.$new_min.out
./script_data_corrected.sh $1 1 2
awk {'print $2, $3, $4, $5, $6'} $1.data2.good > energy_file
paste min_file energy_file >> $1.data1.good
sort -n -k1,1 $1.data1.good > $1.data1.good.sorted
./Spline.x << EOF
$1.data1.good.sorted
Minima_find
EOF
new_energy=$(awk 'NR == 1 {print $2}' $1.data2.good)

# Calculate energy difference
energy_diff=$(awk "BEGIN {print $new_energy - $old_energy}")
echo $energy_diff
if [ $(awk "BEGIN {if ($energy_diff < 0) print 1; else print 0}") -eq 1 ]; then
    energy_diff=$(awk "BEGIN {print -1*$energy_diff}")
fi
echo $energy_diff
echo $new_min
# Compare with a threshold
# REVISE THIS PART
if [ $(awk "BEGIN {if ($energy_diff < $threshold) print 1; else print 0}") -eq 1 ]; then
    echo "Convergence reached. Energy difference: $energy_diff"
    new_min=$(awk 'NR == 1 || $6 < min {min = $6; val = $1} END {print val}' Minima.dat)
    echo $new_min > min_file
    qp set cap eta_cap $new_min
    qp run diagonalize_h_cap > $1.scan.out
    cp $1.scan.out $1.scan.script.$new_min.out
    ./script_data_corrected.sh $1 1 2
    awk {'print $2, $3, $4, $5, $6, "FINAL_MINIMA"'} $1.data2.good > energy_file
    paste min_file energy_file >> $1.data1.good
    sort -n -k1,1 $1.data1.good > $1.data1.good.sorted
    break
fi

old_energy=$new_energy
done

rm -f $1.data2.good
rm -f energy_file
rm -f min_file
rm -f Minima.dat



