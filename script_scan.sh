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
./script_data_corrected.sh $1 2 1
./a.out << EOF > Minimum.dat
$1.data1.good
Minima_find
EOF
threshold=0.00003
qp set cap eta_step_size 0.
qp set cap n_steps_cap  1
old_energy=$(awk 'NR == 3 {print $2}' Minimum.dat)



while true; do

new_min=$(awk 'NR == 3 {print $1}' Minimum.dat)
qp set cap eta_cap $new_min
qp run diagonalize_h_cap > $1.scan.out
./script_data_corrected.sh $1 1 2
cat $1.data2.good >> $1.data1.good
sort -n -k1,1 $1.data1.good > $1.data1.good.sorted
./a.out << EOF > Minimum.dat
$1.data1.good.sorted
Minima_find
EOF
new_energy=$(awk 'NR == 3 {print $2}' Minimum.dat)

# Calculate energy difference
energy_diff=$(awk "BEGIN {print $new_energy - $old_energy}")
echo $energy_diff
if [ $(awk "BEGIN {if ($energy_diff < 0) print 1; else print 0}") -eq 1 ]; then
    energy_diff=$(awk "BEGIN {print -1*$energy_diff}")
fi
echo $energy_diff
# Compare with a threshold
if [ $(awk "BEGIN {if ($energy_diff < $threshold) print 1; else print 0}") -eq 1 ]; then
    echo "Convergence reached. Energy difference: $energy_diff"
    new_min=$(awk 'NR == 3 {print $1}' Minimum.dat)
    qp set cap eta_cap $new_min
    qp run diagonalize_h_cap > $1.scan.out
    ./script_data_corrected.sh $1 1 2
    cat $1.data2.good >> $1.data1.good
    sort -n -k1,1 $1.data1.good > $1.data1.good.sorted
    break
fi

old_energy=$new_energy
done





