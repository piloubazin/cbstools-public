#!/bin/bash
if [ -z "$4" ]; then 
	echo usage: $0 input_image level_parameter iter_parameter output_image
	exit
fi

# run the script assuming Lipsia is in the path
input="$1"
level="$2"
iter="$3"
output="$4"

# 1. convert to vista format
echo vvinidi -in $input -out $input.v -repn float
vvinidi -in $input -out $input.v -repn float

# 2. run the denoising
echo vdenoise -in $input.v -out $output.v -level $level
vdenoise -in $input.v -out $output.v -level $level

# 3. run the speckle removal
echo vsrad -in $output.v -out $output.v -iter $iter
vsrad -in $output.v -out $output.v -iter $iter

# 4. convert back to nifti format
echo vvinidi -in $output.v -out $output -repn float
vvinidi -in $output.v -out $output -repn float

# 5. clean-up intermediate files
rm $input.v
rm $output.v
