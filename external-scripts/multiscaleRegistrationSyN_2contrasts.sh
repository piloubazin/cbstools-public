#!/bin/bash

# args 
#  1: Input - Absolute Path for Moving Image (Subject)
#  2: Input - Absolute Path for Fixed Image (Target)
#  3: Input - Contrast Weight
#  4: Input - Absolute Path for Moving Image 2 (Subject)
#  5: Input - Absolute Path for Fixed Image 2 (Target)
#  6: Input - Contrast Weight 2
#  7: Input - Absolute Path for Coordinate Image (both Subject and Target)
#  8: Input - Number of Coarse Iterations
#  9: Input - Number of Medium Iterations
#  10: Input - Number of Fine Iterations
#  11: Input - Run Affine Step Before Registration
#  12: Output - Absolute Path for Registration Result (Deform Subject)
#  13: Output - Absolute Path for Registration Result 2 (Deform Subject)
#  14: Output - Absolute Path for Forward Mapping (Subject to Target)
#  15: Output - Absolute Path for Inverse Mapping (Target to Subject)
#  16: Output - Absolute Path for working directory

moving="$1"
fixed="$2"
weight="$3"
moving2="$4"
fixed2="$5"
weight2="$6"
coord="$7"
coarseItr="$8"
mediumItr="$9"
fineItr="${10}"
affineFirst="${11}"
outResult="${12}"
outResult2="${13}"
outMap="${14}"
outInvMap="${15}"
destdir="${16}"

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
metricTag2="-m CC[$fixed2,$moving2,$weight2,2]"

interpTag=""

# Calls ANTS
ANTS 3 $metricTag $metricTag2 -o $output -i $coarseItr"x"$mediumItr"x"$fineItr -r Gauss[3,3] -t SyN[.25] $offAffineTag --use-all-metrics-for-convergence

# Transform the moving images
WarpImageMultiTransform 3 $moving $outResult -R $fixed $interpTag $synWarp $synAffine
WarpImageMultiTransform 3 $moving2 $outResult2 -R $fixed2 $interpTag $synWarp $synAffine

# Transform the moving coordinates image
WarpTimeSeriesImageMultiTransform 4 $coord $outMap -R $fixed $synWarp $synAffine

# Transform the fixed coordinates image
WarpTimeSeriesImageMultiTransform 4 $coord $outInvMap -R $moving -i $synAffine $synInverseWarp

# delete temp files
rm -f $synAffine
rm -f $synWarp
rm -f $synInverseWarp

