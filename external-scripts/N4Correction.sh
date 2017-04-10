#!/bin/bash
# args 
#  1: Input - Absolute Path for Input Image
#  2: Input - Shrink factor (1/2/3/4...)
#  3: Output - Absolute Path for Corrected Image
#  4: Output - Absolute Path for Bias Field

input="$1"
shrink="$2"
corrected="$3"
biasField="$4"


#Calls N4Correction
N4BiasFieldCorrection -i $input -o [$corrected,$biasField] --shrink-factor $shrink
