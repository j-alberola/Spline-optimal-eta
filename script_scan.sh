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
0
EOF

threshold=0.00003
qp set cap eta_step_size 0.
qp set cap n_steps_cap  1



#  ADD A LOOP HERE
#  MAYBE REMOVE THE DIFFERENT 2 CASES, JUST 1 IS ENOUGH

iter = 1
while iter < 20 ; do


if [ -f "Minima_converged.dat" ]; then
awk -v thresh=0.01 '
BEGIN {
    while ((getline < "Minima_converged.dat") > 0) {
        file2_values[$1] = 1;  # Store all values from file2
    }
}
{
    for (val in file2_values) {
        diff = ($1 > val) ? $1 - val : val - $1;  # Absolute difference
        if (diff <= $1 * thresh) {
            next;  # Skip the line if within 1% interval
        }
    }
    print $0;  # Print the line if no match within 1%
}' "Minima.dat" > aux_file

cat aux_file > Minima.dat
rm -f aux_file

fi
#
#  If there are no remaining minima, exit from the loop
#
if [ -z "$(head -n 1 Minima.dat)" ]; then
    break
fi

fit_real_energy=$(awk 'NR == 1 {print $2}' Minima.dat)
fit_imag_energy=$(awk 'NR == 1 {print $3}' Minima.dat)
fit_real_deriv=$(awk 'NR == 1 {print $4}' Minima.dat)
fit_imag_deriv=$(awk 'NR == 1 {print $5}' Minima.dat)
new_min=$(awk 'NR == 1 {print $1}' Minima.dat)
qp set cap eta_cap $new_min
echo $new_min > min_file
qp run diagonalize_h_cap > $1.scan.out
#REMOVE WHEN NOT DEBUGGINMG
cp $1.scan.out $1.scan.script.$new_min.out
./script_data_corrected.sh $1 1 2
qp_real_energy=$(awk 'NR == 1 {print $2}' $1.data2.good)
qp_imag_energy=$(awk 'NR == 1 {print $3}' $1.data2.good)
qp_real_deriv=$(awk 'NR == 1 {print $4}' $1.data2.good)
qp_imag_deriv=$(awk 'NR == 1 {print $5}' $1.data2.good)
#Calculating differences between fit and qp energies and derivatives in absolute value
real_energy_diff=$(awk "BEGIN {print $qp_real_energy - $fit_real_energy}")
if [ $(awk "BEGIN {if ($real_energy_diff < 0) print 1; else print 0}") -eq 1 ]; then
    real_energy_diff=$(awk "BEGIN {print -1*$real_energy_diff}")
fi
imag_energy_diff=$(awk "BEGIN {print $qp_imag_energy - $fit_imag_energy}")
if [ $(awk "BEGIN {if ($imag_energy_diff < 0) print 1; else print 0}") -eq 1 ]; then
    imag_energy_diff=$(awk "BEGIN {print -1*$imag_energy_diff}")
fi
real_deriv_diff=$(awk "BEGIN {print $qp_real_deriv - $fit_real_deriv}")
if [ $(awk "BEGIN {if ($real_deriv_diff < 0) print 1; else print 0}") -eq 1 ]; then
    real_deriv_diff=$(awk "BEGIN {print -1*$real_deriv_diff}")
fi
imag_deriv_diff=$(awk "BEGIN {print $qp_imag_deriv - $fit_imag_deriv}")
if [ $(awk "BEGIN {if ($imag_deriv_diff < 0) print 1; else print 0}") -eq 1 ]; then
    imag_deriv_diff=$(awk "BEGIN {print -1*$imag_deriv_diff}")
fi


# If all the conditions are fulfilled, exit from the loop
if [ $(awk "BEGIN {if ($real_energy_diff < $threshold && $imag_energy_diff < $threshold && $real_deriv_diff < $threshold && $imag_deriv_diff < $threshold) print 1; else print 0}") -eq 1 ]; then
        echo "Convergence reached"
        echo $new_min
        awk {'print $2, $3, $4, $5, $6, FINAL MINIMA'} $1.data2.good > energy_file
        paste min_file energy_file >> $1.data1.good
        paste min_file energy_file >> Minima_converged.dat
        sort -n -k1,1 $1.data1.good > $1.data1.good.sorted
	./Spline.x << EOF
        $1.data1.good
        Minima_find
        0
EOF

else
        awk {'print $2, $3, $4, $5, $6'} $1.data2.good > energy_file
        paste min_file energy_file >> $1.data1.good
        sort -n -k1,1 $1.data1.good > $1.data1.good.sorted
        ./Spline.x << EOF
        $1.data1.good
        Minima_find
        0
EOF
fi

iter = $iter + 1
done

rm -f $1.data2.good
rm -f energy_file
rm -f min_file
rm -f Minima.dat

