# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 13:41:28 2016

@author: chris
"""

import numpy as np
import nibabel as nb
import os
import sys

sys.path.append('/home/chris/Documents/code/python/cbstools-python/cbstoolsjcc-3.1.0.1-py2.7-linux-x86_64.egg')
import cbstoolsjcc as cj

from defaults import * #ATLAS_DIR and TOPOLOGY_LUT_DIR

def normalise(img_d):
    return (img_d - np.min(img_d))/np.max(img_d)

def setup_JVM(JVM_initialheap = '4000M', JVM_maxheap = '4000M'):
    """
    initialise the JVM and set all of the base with reasonable defaults for memory
    :param JVM_initialheap:
    :param JVM_maxheap:
    :return:
    """
    try:
        res=cj.initVM(initialheap=JVM_initialheap,maxheap=JVM_maxheap)
        print(res)
        print("Java virtual machine successfully started.")
    except ValueError:
        print("A java virtual machine is already running.")


def ExtractBrainRegion():
    pass

def Mp2rageSkullStripping():
    pass

def IntensityBackgroundEstimator():
    pass

def SurfaceProbabilityToLevelset():
    pass

def MGDMBrainSegmentation(input_filename_type_list, output_dir = None, num_steps = 5, atlas_file=None, topology_lut_dir = None):
    """
    Perform MGDM segmentation
    :param input_filename_type_list: list of [[fname1,type1],[fname2,type2],...] - for a maximum of 4 inputs
    :param output_dir: full path to the output directory
    :param num_steps: number of steps for (default 5, set to 0 for testing)
    :param atlas_file: full path to the atlas file, default set in defaults.py
    :param topology_lut_dir: full path to the directory with the topology files, default set in defaults.py
    :return:
    """
    print("Thank you for choosing the MGDM segmentation from the cbstools for your brain segmentation needs")
    print("Sit back and relax, let the magic of algorithms happen...")
    print("")
    if output_dir is None:
        output_dir = os.path.dirname(input_filename_type_list[0][0])
    if atlas_file is None:
        atlas = os.path.join(ATLAS_DIR,'brain-atlas-3.0.3.txt')
    else:
        atlas = atlas_file

    if topology_lut_dir is None:
        topology_lut_dir = TOPOLOGY_LUT_DIR  # grabbing this from the default settings in defaults.py
    else:
        if not(topology_lut_dir[-1] == os.sep): #if we don't end in a path sep, we need to make sure that we add it
            topology_lut_dir += os.sep
    print("Atlas file: " + atlas)
    print("Topology LUT durectory: " + topology_lut_dir)
    print("")

    if not any(isinstance(el, list) for el in input_filename_type_list): #make into list of lists
        input_filename_type_list = [input_filename_type_list]


    #now we setup the mgdm specfic settings
    mgdm = cj.BrainMgdmMultiSegmentation2()
    mgdm.setAtlasFile(atlas)
    mgdm.setTopologyLUTdirectory(topology_lut_dir)

    mgdm.setOutputImages('segmentation');
    mgdm.setOrientations(mgdm.AXIAL, mgdm.R2L, mgdm.A2P, mgdm.I2S);
    mgdm.setAdjustIntensityPriors(False)  # default is True
    mgdm.setComputePosterior(False)
    mgdm.setDiffuseProbabilities(False)
    mgdm.setSteps(num_steps)
    mgdm.setTopology('wcs')  # {'wcs','no'} no=off for testing, wcs=default

    for idx,con in enumerate(input_filename_type_list):
        print("Input files and filetypes:")
        print("  " + str(idx+1) + " "),
        print(con)
        fname = con[0]
        type = con[1]
        img = nb.load(fname)
        d = img.get_data()
        if idx+1 == 1:
            # we use the first image to set the dimensions and resolutions
            res = img.header.get_zooms()
            res = [a1.item() for a1 in res]  # cast to regular python float type
            mgdm.setDimensions(d.shape[0], d.shape[1], d.shape[2])
            mgdm.setResolutions(res[0], res[1], res[2])

            # keep the shape and affine from the first image for saving
            d_shape = np.array(d.shape)
            d_aff = img.affine
            out_root_fname = os.path.basename(fname)[0:os.path.basename(fname).find('.')] #assumes no periods in filename, :-/

            mgdm.setContrastImage1(cj.JArray('float')((d.flatten('F')).astype(float)))
            mgdm.setContrastType1(type)
        elif idx+1 == 2:
            mgdm.setContrastImage2(cj.JArray('float')((d.flatten('F')).astype(float)))
            mgdm.setContrastType2(type)
        elif idx + 1 == 3:
            mgdm.setContrastImage3(cj.JArray('float')((d.flatten('F')).astype(float)))
            mgdm.setContrastType3(type)
        elif idx + 1 == 4:
            mgdm.setContrastImage4(cj.JArray('float')((d.flatten('F')).astype(float)))
            mgdm.setContrastType4(type)
    try:
        print("Executing MGDM on your inputs")
        print("Don't worry, the magic is happening!")
        mgdm.execute()
        print(os.path.join(output_dir, out_root_fname + '_seg_cjs.nii.gz'))


        # outputs
        # reshape fortran stype to convert back to the format the nibabel likes
        seg_im = np.reshape(np.array(mgdm.getSegmentedBrainImage(), dtype=np.uint32), d_shape,'F')
        lbl_im = np.reshape(np.array(mgdm.getPosteriorMaximumLabels4D(), dtype=np.uint32), d_shape, 'F')
        ids_im = np.reshape(np.array(mgdm.getSegmentedIdsImage(), dtype=np.uint32), d_shape, 'F')

        # save
        out_im = nb.Nifti1Image(seg_im, d_aff)
        nb.save(out_im, os.path.join(output_dir, out_root_fname + '_seg_cjs.nii.gz'))
        out_im = nb.Nifti1Image(lbl_im, d_aff)
        nb.save(out_im, os.path.join(output_dir, out_root_fname + '_lbl_cjs.nii.gz'))
        out_im = nb.Nifti1Image(ids_im, d_aff)
        nb.save(out_im, os.path.join(output_dir, out_root_fname + '_ids_cjs.nii.gz'))

        print("Data stored in: " + output_dir)
    except:
        print("--- MGDM failed. Go cry. ---")
    print("Execution completed")

    #return (img,out_im)