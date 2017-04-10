#!/bin/bash

# args 
#  1: Input - Absolute Path for Moving Image (Subject)
#  2: Input - Absolute Path for Fixed Image (Target)
#  3: Input - Contrast Weight
#  4: Input - Absolute Path for Coordinate Image (both Subject and Target)
#  5: Input - Number of Coarse Iterations
#  6: Input - Number of Medium Iterations
#  7: Input - Number of Fine Iterations
#  8: Input - Run Affine Step Before Registration
#  9: Output - Absolute Path for Registration Result (Deform Subject)
#  10: Output - Absolute Path for Forward Mapping (Subject to Target)
#  11: Output - Absolute Path for Inverse Mapping (Target to Subject)
#  12: Output - Absolute Path for working directory

moving="$1"
fixed="$2"
weight="$3"
coord="$4"
coarseItr="$5"
mediumItr="$6"
fineItr="$7"
affineFirst="$8"
outResult="$9"
outMap="${10}"
outInvMap="${11}"
destdir="${12}"

output="$destdir/syn.nii"
synWarp="$destdir/synWarp.nii"
synInverseWarp="$destdir/synInverseWarp.nii"
synAffine="$destdir/synAffine.txt"

# Perform affine registration first
if [ "$affineFirst" == "true" ]
then
		offAffineTag=""
else
    offAffineTag="--number-of-affine-iterations 0"
fi

metricTag="-m CC[$fixed,$moving,$weight,2]"

interpTag=""

# Calls ANTS 
ANTS 3 $metricTag -o $output -i $coarseItr"x"$mediumItr"x"$fineItr -r Gauss[3,3] -t SyN[.25] $offAffineTag 

# Transform the moving images
WarpImageMultiTransform 3 $moving $outResult -R $fixed $interpTag $synWarp $synAffine

# Transform the moving coordinates image
WarpTimeSeriesImageMultiTransform 4 $coord $outMap -R $fixed $synWarp $synAffine

# Transform the fixed coordinates image
WarpTimeSeriesImageMultiTransform 4 $coord $outInvMap -R $moving -i $synAffine $synInverseWarp

# delete temp files
rm -f $synAffine
rm -f $synWarp
rm -f $synInverseWarp

