# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 13:41:28 2016

@author: Christopher J. Steele
steele{AT}cbs{dot}mpg{dot}de
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


def seg_erode(seg_d, iterations=1, background_idx=1,
                  structure=None, min_vox_count=5, seg_null_value=0,
                  VERBOSE=False):
    """
    Binary erosion of integer type segmentation data (np.array) with options

    :param seg_d:           np.array of segmentation, integers
    :param iterations:      number of erosion iterations
    :param background_idx:  value for background index, currently ignored (TODO: remove)
    :param structure:       binary structure for erosion from scipy.ndimage (ndimage.morphology.generate_binary_structure(3,1))
    :param min_vox_count:   minimun number of voxels to allow to be in a segmentation, if less, does not erode
    :param seg_null_value:  value to set as null for binary erosion step (i.e., a value NOT in your segmentation index)
    :param VERBOSE:         spit out loads of text to stdout, because you can.
    :return: seg_shrunk_d   eroded version of segmentation
    """

    import scipy.ndimage as ndi
    import numpy as np

    if structure is None:
        structure = ndi.morphology.generate_binary_structure(3, 1)
    if seg_null_value == 0:
        seg_shrunk_d = np.zeros_like(seg_d)
        temp_d = np.zeros_like(seg_d)
    else:
        seg_shrunk_d = np.ones_like(seg_d) * seg_null_value
        temp_d = np.ones_like(seg_d) * seg_null_value

    seg_idxs = np.unique(seg_d)

    if seg_null_value in seg_idxs:
        print("Shit, your null value is also an index. This will not work.")
        print("Set it to a suitably strange value that is not already an index. {0,999}")
        return None
    if VERBOSE:
        print("Indices:")
    for seg_idx in seg_idxs:
        print(seg_idx),
        if (background_idx is not None) and (background_idx == seg_idx):
            seg_shrunk_d[seg_d == seg_idx] = seg_idx  # just set the value to the bckgrnd value, and be done with it
            if VERBOSE:
                print("[bckg]"),
        else:
            temp_d[seg_d == seg_idx] = 1
            for idx in range(0, iterations):  # messy, does not exit the loop when already gone too far.
                temp_temp_d = ndi.binary_erosion(temp_d, iterations=1, structure=structure)
                if np.sum(temp_temp_d) >= min_vox_count:
                    temp_d = temp_temp_d
                    if VERBOSE:
                        print("[y]"),
                else:
                    if VERBOSE:
                        print("[no]"),
            seg_shrunk_d[temp_d == 1] = seg_idx
            temp_d[:, :, :] = seg_null_value
            if VERBOSE:
                print(seg_idx)
        if VERBOSE:
            print("")
    return seg_shrunk_d


def extract_metrics_from_seg(seg_d, metric_d, seg_idxs=None,norm_data=True,
                             background_idx=1, seg_null_value=0,
                             percentile_top_bot=[75, 25],
                             return_normed_metric_d=False):
    """
    Extract median and interquartile range from metric file given a co-registered segmentation

    :param seg_d:                   segmentation data (integers)
    :param metric_d:                metric data to extract seg-specific values from
    :param seg_idxs:                indices of segmentation, usually taken from LUT but can be generated based on seg_d
    :param norm_data:               perform data normalisation on metric_d prior to extracting values from metric
    :param background_idx:          index for background data, currently treated as just another index (TODO: remove)
    :param seg_null_value:          value to set as null for binary erosion step, not included in metric extraction
    :param percentile_top_bot:      top and bottom percentiles to extract from each seg region
    :param return_normed_metric_d:  return the normalised metric as an np matrix, must also set norm_data=True
    :return: seg_idxs, res          segmentation indices and results matrix of median, 75, 25 percentliles
             (metric_d)             optional metric_d scaled between 0 and 1
    """
    import numpy as np
    import scipy
    if seg_idxs is None:
        seg_idxs = np.unique(seg_d)
    if (seg_null_value is not None) and (seg_null_value in seg_idxs): #remove the null value from the idxs so we don't look
        np.delete(seg_idxs,np.where(seg_idxs==seg_null_value))
    res = np.zeros((len(seg_idxs), 3))

    if norm_data:  # rescale the data to 0
        if background_idx is not None:  # we need to exclude the background data from the norming
            metric_d[seg_d != background_idx] = (metric_d[seg_d != background_idx] - np.min(
                metric_d[seg_d != background_idx])) / (np.max(metric_d[seg_d != background_idx]) - np.min(
                metric_d[seg_d != background_idx]))
        else:
            metric_d = (metric_d - np.min(metric_d)) / (np.max(metric_d) - np.min(metric_d))

    for idx, seg_idx in enumerate(seg_idxs):
        # no longer required, since we also want data from the background
        # if (background_idx is not None) and ((seg_idx == background_idx) or (seg_idx == seg_null_value)):
        #     res[idx, :] = [0, 0, 0]
        # else:
            d_1d = np.ndarray.flatten(metric_d[seg_d == seg_idx])
            res[idx, :] = [np.median(d_1d),
                           np.percentile(d_1d, np.max(percentile_top_bot)),
                           np.percentile(d_1d, np.min(percentile_top_bot))]
    if return_normed_metric_d:
        return seg_idxs, res, metric_d
    else:
        return seg_idxs, res


def extract_lut_priors_from_atlas(atlas_file,contrast_name):
    """
    Given an MGDM segmentation priors atlas file, extract the lut and identify the start index (in the file) of the
    contrast of interest, and the number of rows of priors that it should have. Returns pandas dataframe of lut,
    contrast index, number of rows in prior definition, and pd.DataFrame of priors,

    :param atlas_file:      full path to atlas file for lut and metric index extraction
    :param contrast_name:   intensity prior contrast name as listed in the metric file
    :return: lut, con_idx, lut_rows, priors
    """
    import pandas as pd
    import os

    fp = open(os.path.join(ATLAS_DIR, atlas_file))
    for i, line in enumerate(fp):
        if "Structures:" in line:  # this is the beginning of the LUT
            lut_idx = i
            lut_rows = map(int, [line.split()[1]])[0]
            #    if "Topology Atlas:" in line: #the end of the LUT #OR you could use the value in the line?
            #        lut_idx[1] = i-2
        if "Intensity Prior:" in line:
            if contrast_name in line:
                con_idx = i
    fp.close()

    # dump into pandas dataframe
    lut = pd.read_csv(os.path.join(ATLAS_DIR, atlas_file), sep="\t+",
                      skiprows=lut_idx + 1, nrows=lut_rows, engine='python',
                      names=["Index", "Type"])

    # con_idx[1] = len(lut) #total number is the same length as the lut
    priors = pd.read_csv(os.path.join(ATLAS_DIR, atlas_file), sep="\t+",
                         skiprows=con_idx + 1, nrows=lut_rows, engine='python',
                         names=["Median", "Spread", "Weight"])
    return lut,con_idx,lut_rows,priors

def generate_group_intensity_priors(orig_seg_files,metric_files,metric_contrast_name,
                                    atlas_file,new_atlas_file_head,erosion_iterations=1,seg_iterations=1,
                                    output_dir=None):
    # generates group intensity priors for metric_files based on orig_seg files (i.e., orig_seg could be Mprage3T and metric_files could be DWIFA3T)
    # does not do the initial segmentation for you, that needs to be done first :-)
    # we assume that you already did due-diligence and have matched lists of inputs (orig_seg_files and metric_files)

    import os
    import nibabel as nb
    import numpy as np
    import pandas as pd

    [lut,con_idx,lut_rows,priors] = extract_lut_priors_from_atlas(atlas_file, metric_contrast_name)
    seg_idxs = lut.Index
    all_Ss_priors_list_median = np.array(seg_idxs)
    all_Ss_priors_list_spread = np.array(seg_idxs)
    seg_null_value = 0 #value to fill in when we are NOT using the voxels at all (not background and not other index)
    background_idx = 1
    min_quart_diff = 0.10 #minimun spread allowed in priors atlas
    new_atlas_files = []

    # make a list if we only input one dataset
    if len(orig_seg_files) == 1:
        orig_seg_files = [orig_seg_files]
    if len(metric_files) == 1:
        metric_files = [metric_files]

    if not(len(orig_seg_files) == len(metric_files)):
        print("You do not have the same number of segmentation and metric files. Bad!")
        print("Exiting")
        return

    for seg_iter in range(0,seg_iterations):
        seg_iter_text = str(seg_iter+1).zfill(3) #text for naming files etc?
        new_atlas_file = os.path.join(ATLAS_DIR,new_atlas_file_head+"_"+seg_iter_text)
        print("Running segmentation iteration: " + seg_iter_text)
        for idx, seg_file in enumerate(orig_seg_files):
            if seg_iter > 0: #TODO: this logic will need to be cleaned up when you have it iterating over segs
                #TODO: call segmentation code, using the priors that we just created, good job!
                pass
            else: #the first time, only extracting the metric, next time we actually do some segs on the file
                metric_file = metric_files[idx]
                img=nb.load(metric_file)
                d_metric = img.get_data()
                a_metric = img.affine #not currently using the affine and header, but could also output the successive steps
                a_header = img.header
                print(seg_file)
                d_seg = nb.load(seg_file).get_data()

                #erode our data
                if erosion_iterations>0:
                    d_seg_ero = seg_erode(d_seg,iterations=erosion_iterations,
                                          background_idx=background_idx,
                                          seg_null_value=seg_null_value)

                #extract summary metrics (median, 75 and 25 percentile) from metric file
                [seg_idxs, seg_stats] = extract_metrics_from_seg(d_seg_ero, d_metric, seg_idxs=seg_idxs,
                                                                 seg_null_value=seg_null_value,
                                                                 return_normed_metric_d=False)

                prior_medians = seg_stats[:, 0]
                prior_quart_diffs = np.squeeze(np.abs(np.diff(seg_stats[:, 1:3])))
                prior_quart_diffs[prior_quart_diffs < min_quart_diff] = min_quart_diff

            #now place this output into a growing array for use on the group level
            #print(np.shape(all_Ss_priors_list_median))
            #print(np.shape(prior_medians))
            all_Ss_priors_list_median = np.vstack((all_Ss_priors_list_median, prior_medians))
            all_Ss_priors_list_spread = np.vstack((all_Ss_priors_list_spread, prior_quart_diffs))

        #combine the output into a 3d image stack if we have more than one iteration to process
        if seg_iter == 0:
            iter_Ss_priors_median = all_Ss_priors_list_median
            iter_Ss_priors_spread = all_Ss_priors_list_spread
        elif seg_iter > 0:
            # stack to 3d if we have more than one iteration
            iter_Ss_priors_median = np.dstack((iter_Ss_priors_median, all_Ss_priors_list_median))
            iter_Ss_priors_spread = np.dstack((iter_Ss_priors_spread, all_Ss_priors_list_spread))
    #TODO: decide if you should do the looping over the segmentations here or in another function
    # priors_new = pd.DataFrame.copy(priors)
    # for idx in lut.Index:
    #     priors_new[lut["Index"] == idx] = [prior_medians[seg_idxs == idx], prior_quart_diffs[seg_idxs == idx],
    #                                        1]
    # priors_new.head()
    # priors_new_string = priors_new.to_csv(sep="\t", header=False, float_format="%.2f")
    # priors_new_string_lines = priors_new_string.split("\n")[
    #                           0:-1]  # convert to list of lines, cut the last empty '' line
    #
    # # write out the newly altered priors list to new_atlas_file
    # #TODO: give option to have it not in the ATLAS_DIR
    # fp = open(os.path.join(ATLAS_DIR, atlas_file))
    # fp_new = open(os.path.join(ATLAS_DIR, new_atlas_file), "w")
    # ii = 0
    # # only replace the lines that we changed
    # for i, line in enumerate(fp):
    #     if i > con_idx and i < con_idx + lut_rows:
    #         fp_new.write(priors_new_string_lines[ii] + "\n")
    #         ii += 1
    #     else:
    #         fp_new.write(line)
    # fp.close()
    # fp_new.close()
    # print("New atlas file written to: " + fp_new.name)

    new_atlas_files.append(new_atlas_file)
    return iter_Ss_priors_median, iter_Ss_priors_spread, new_atlas_files