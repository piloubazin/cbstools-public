#!/bin/bash
if [ -z "$9" ]; then 
	echo usage: $0 full_output_directory_path input_image isMPRAGE noSurface segmentation inner_surface_rh inner_surface_lh outer_surface_rh outer_surface_lh
	exit
fi
export SUBJECTS_DIR=$1
export BASENAME=$(basename $2)
export BASENAME=${BASENAME%.*}
export MPRAGE=$3
export NOSURFACE=$4
export SEGMENTATION=$5
export INNER_RH=$6
export INNER_LH=$7
export OUTER_RH=$8
export OUTER_LH=$9
export FREESURFER_HOME="/opt/freesurfer"
echo "Set up Freesurfer environment"
if [ ! -d $FREESURFER_HOME ];then
	echo "Could not locate correct FreeSurfer version"
	exit
fi
export FREESURFER_DIR="/opt/freesurfer/bin"
if [ ! -d $FREESURFER_DIR ];then
	echo "Could not locate correct FreeSurfer binaries"
	exit
fi
source $FREESURFER_HOME/SetUpFreeSurfer.sh
echo "subject dir:" $SUBJECTS_DIR
echo processing $BASENAME ....

if [ ! -d $SUBJECTS_DIR ];then
	echo $SUBJECTS_DIR not found	
	if [ -d ${SUBJECTS_DIR%/*} ];then
	
	echo $SUBJECTS_DIR created
		mkdir $SUBJECTS_DIR
	else
		echo can not creat output directory
		exit
	fi
fi

if [ ! -e $2 ]; then
	echo $2 not found
	exit
fi
$FREESURFER_DIR/mri_convert $2 $SUBJECTS_DIR/$BASENAME.mgz
cd $SUBJECTS_DIR
if [ "$MPRAGE" == "true" ]; then
	echo "MPRAGE"
	mprageTag="-mprage"	
else
	mprageRag=""
fi
if [ "$NOSURFACE" == "true" ]; then
	echo "Just Segmentation"
	$FREESURFER_DIR/recon-all -subject $BASENAME/ -i $BASENAME.mgz -autorecon1 -subcortseg -segstats $mprageTag
	cd $1/$BASENAME/mri
	$FREESURFER_DIR/mri_label2vol --seg aseg.mgz --temp rawavg.mgz --o $BASENAME_FS_seg.mgz --regheader aseg.mgz
	$FREESURFER_DIR/mri_convert $BASENAME_FS_seg.mgz $SEGMENTATION
else
	$FREESURFER_DIR/recon-all -subject $BASENAME/ -i $BASENAME.mgz -all $mprageTag
	$FREESURFER_DIR/mris_convert $1/$BASENAME/surf/rh.white $INNER_RH
	$FREESURFER_DIR/mris_convert $1/$BASENAME/surf/lh.white $INNER_LH
	$FREESURFER_DIR/mris_convert $1/$BASENAME/surf/rh.pial $OUTER_RH
	$FREESURFER_DIR/mris_convert $1/$BASENAME/surf/lh.pial $OUTER_LH
	cd $1/$BASENAME/mri
	$FREESURFER_DIR/mri_label2vol --seg aseg.mgz --temp rawavg.mgz --o $BASENAME_FS_seg.mgz --regheader aseg.mgz
	$FREESURFER_DIR/mri_convert $BASENAME_FS_seg.mgz $SEGMENTATION
fi
