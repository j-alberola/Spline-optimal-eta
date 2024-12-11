#!/bin/bash


##
## ADD AFTER AN ADDITIONAL SCRIPT TO FIND THE COLUMNS CORRESPONDING TO THE LOWST
## IMAGINARY ENERGY, INSTEAD OF SELECTING THEM IN THE CODE
##

grep -A$(( $2+10 )) 'CAP energies' $1.scan.out | tail -$2 | awk {'print $1, $4, $5'} > $1.data


grep -A$(( $2+7 )) 'Energy corrections' $1.scan.out | tail -$2 | awk {'print $4, $5'} > $1.corrections.data


paste $1.data $1.corrections.data >$1.data$3.good

