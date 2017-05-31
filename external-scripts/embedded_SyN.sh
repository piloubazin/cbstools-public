#!/bin/bash
# args 
#  2: Input - Absolute Path for Fixed Image (Target)
#  3: Input - Absolute Path for Moving Image (Subject)
#  4: Input - Absolute Path for Fixed Coordinate Image (Target)
#  5: Input - Absolute Path for Moving Coordinate Image (Subject)
#  6: Input - Number of Coarse Iterations
#  7: Input - Number of Medium Iterations
#  8: Input - Number of Fine Iterations
#  9: Input - Run Affine Step Before Registration
#  10: Input - Type of Cost Function 
#  11: Input - Type of Interpolation (moving image)
#  12: Output - Absolute Path for Registration Result (Deform Subject)
#  13: Output - Absolute Path for Forward Mapping (Subject to Target)
#  14: Output - Absolute Path for Inverse Mapping (Target to Subject)

dir="$1"
fixed="$2"
moving="$3"
coordf="$4"
coordm="$5"
coarseItr="$6"
mediumItr="$7"
fineItr="$8"
affineFirst="$9"
costFunction="${10}"
interpolation="${11}"
outResult="${12}"
outMap="${13}"
outInvMap="${14}"


if [ "$affineFirst" == "true" ]
then
        offAffineTag=""
else
        offAffineTag="--number-of-affine-iterations 0"
fi


if [ "$costFunction" == "Cross Correlation" ]
then
        metricTag="CC[$fixed,$moving,1,5]"
elif [ "$costFunction" == "Mutual Information" ]
then
        metricTag="MI[$fixed,$moving,1,64]"
fi

if [ "$interpolation" == "Nearest Neighbor" ]
then
        interpTag="--use-NN"
elif [ "$interpolation" == "Linear" ]
then
        interpTag=""
elif [ "$interpolation" == "BSpline" ]
then
        interpTag="--use-BSpline"
fi

# Make sure the used version of ANTS is correct (not clean; should be an afs version instead, but that needs to be tested)
export ANTSPATH=/home/piloubazin/Software/ants-2.1.0-redhat/ANTS/

export PATH=$ANTSPATH:$PATH

# Set working directory to output
cd $dir
#Calls ANT and then transform moving (subject) image using deformation field
echo ANTS 3 -m $metricTag -o syn -i $coarseItr"x"$mediumItr"x"$fineItr -r Gauss[3,1] -t SyN[.25] $offAffineTag
ANTS 3 -m $metricTag -o syn -i $coarseItr"x"$mediumItr"x"$fineItr -r Gauss[3,1] -t SyN[.25] $offAffineTag
# Transforms the moving image
echo WarpImageMultiTransform 3 $moving $outResult -R $fixed $interpTag synWarp.nii synAffine.txt
WarpImageMultiTransform 3 $moving $outResult -R $fixed $interpTag synWarp.nii synAffine.txt
# Transforms the moving coordinates image
echo WarpTimeSeriesImageMultiTransform 4 $coordm $outMap -R $fixed synWarp.nii synAffine.txt
WarpTimeSeriesImageMultiTransform 4 $coordm $outMap -R $fixed synWarp.nii synAffine.txt
# Transforms the fixed coordinates image
echo WarpTimeSeriesImageMultiTransform 4 $coordf $outInvMap -R $moving -i synAffine.txt synInverseWarp.nii
WarpTimeSeriesImageMultiTransform 4 $coordf $outInvMap -R $moving -i synAffine.txt synInverseWarp.nii
# Clean-up
rm -f synAffine.txt synWarp.* synInverseWarp.*

