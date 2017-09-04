import nibabel as nb
import numpy as np
import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('command')
parser.add_argument('-o', '--outdir', nargs=1, help='output directory (default: input directory)')
parser.add_argument('imgfiles', nargs='+', help='input image files')
print(parser.parse_args(sys.argv))

imgfiles = parser.parse_args(sys.argv).imgfiles
outdirs = parser.parse_args(sys.argv).outdir
	
for imgfile in imgfiles:
	print("file to convert: "+imgfile)
	if outdirs==None:
		outdir = os.path.dirname(imgfile)
		print("input/output directory: "+outdir)
	else:
		outdir=outdirs[0]
		print("output directory: "+outdir)

	imgdir = os.path.dirname(imgfile)

	print("opening file: "+imgfile)
	img = nb.load(imgfile)
	data = img.get_data()
	affine = img.affine #or img.get_affine(), which I think is being deprecated?
	header = img.header

	imgshape = np.shape(data)

	print("extracting magnitude and phase")
	dmag = np.split(data,2,axis=3)[0]
	dphs = np.split(data,2,axis=3)[1]

	imgname = os.path.basename(imgfile)
	basename = imgname.split(".")
	magname = outdir+basename[0]+"_mag."+'.'.join(basename[1:])
	phsname = outdir+basename[0]+"_phs."+'.'.join(basename[1:])

	print("saving to: "+magname+" and "+phsname)
	mag = nb.Nifti1Image(dmag, affine, header)
	mag.header['cal_min'] = np.min(dmag)
	mag.header['cal_max'] = np.max(dmag)
	mag.to_filename(magname)

	phs = nb.Nifti1Image(dphs, affine, header)
	phs.header['cal_min'] = np.min(dphs)
	phs.header['cal_max'] = np.max(dphs)
	phs.to_filename(phsname)
