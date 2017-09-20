package de.mpg.cbs.jist.utilities;

import javax.vecmath.Point3i;

import de.mpg.cbs.utilities.CropParameters;
import de.mpg.cbs.utilities.CubicVolumeCropper;
import de.mpg.cbs.libraries.*;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamPointInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;


/*
 * Image Cropper
 * @author Navid Shiee
 * @author Pierre-Louis Bazin
 *	(+minor modifications to allow the region to be defined from a mask)
 *
 */
public class JistUtilitiesCropToROI extends ProcessingAlgorithm{
	private ParamVolume input1Vol;
	private ParamVolume input2Vol;
	private ParamVolume input3Vol;
	private ParamVolume input4Vol;
	private ParamVolume maskVol;

	private ParamVolume output1Vol;
	private ParamVolume output2Vol;
	private ParamVolume output3Vol;
	private ParamVolume output4Vol;

	private ParamInteger borderSize;
	private ParamString roiName;
	private ParamOption		borderType;
	private static final String[] borderTypes = {"image", "zero", "min", "max"};
	
	private ParamPointInteger offset;
	private ParamPointInteger dim;
	private ParamPointInteger dimNew;
	

	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(input1Vol = new ParamVolume("Image 1", null,-1,-1,-1,-1));
		inputParams.add(input2Vol = new ParamVolume("Image 2 (opt)", null,-1,-1,-1,-1));
		inputParams.add(input3Vol = new ParamVolume("Image 3 (opt)", null,-1,-1,-1,-1));
		inputParams.add(input4Vol = new ParamVolume("Image 4 (opt)", null,-1,-1,-1,-1));
		input2Vol.setMandatory(false);
		input3Vol.setMandatory(false);
		input4Vol.setMandatory(false);
		input1Vol.setLoadAndSaveOnValidate(false);
		input2Vol.setLoadAndSaveOnValidate(false);
		input3Vol.setLoadAndSaveOnValidate(false);
		input4Vol.setLoadAndSaveOnValidate(false);
		
		inputParams.add(maskVol = new ParamVolume("ROI Mask"));
		inputParams.add(borderSize = new ParamInteger("BorderSize", 0, Integer.MAX_VALUE, 3));
		inputParams.add(borderType = new ParamOption("Border type", borderTypes));
		borderType.setValue("image");
		inputParams.add(roiName = new ParamString("ROI name"));

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Utilities");
		inputParams.setLabel("Crop To ROI");
		inputParams.setName("CropToROI");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Navid Shiee","","iacl.ece.jhu.edu"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Christine L Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Crops one or more images to include a given mask and a boundary around.");
		info.setLongDescription("Crops one or more images to include a given mask and a boundary around.\n" 
					+"The border type sets what value to use outside the mask.");
		info.setVersion("3.0.3");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}

	protected void createOutputParameters(ParamCollection outputParams) {

		outputParams.add(output1Vol = new ParamVolume("Cropped Image 1",null,-1,-1,-1,-1));
		outputParams.add(output2Vol = new ParamVolume("Cropped Image 2",null,-1,-1,-1,-1));
		outputParams.add(output3Vol = new ParamVolume("Cropped Image 3",null,-1,-1,-1,-1));
		outputParams.add(output4Vol = new ParamVolume("Cropped Image 4",null,-1,-1,-1,-1));
		output2Vol.setMandatory(false);
		output3Vol.setMandatory(false);
		output4Vol.setMandatory(false);
		output1Vol.setLoadAndSaveOnValidate(false);
		output2Vol.setLoadAndSaveOnValidate(false);
		output3Vol.setLoadAndSaveOnValidate(false);
		output4Vol.setLoadAndSaveOnValidate(false);
		
		offset=new ParamPointInteger("Offset");
		dim=new ParamPointInteger("Uncropped Dimensions");
		dimNew=new ParamPointInteger("Cropped Dimensions");
		
		outputParams.add(offset);
		outputParams.add(dim);
		outputParams.add(dimNew);
	}


	protected void execute(CalculationMonitor monitor) {
		CubicVolumeCropper cropper = new CubicVolumeCropper();
		monitor.observe(cropper);
		
		ImageDataFloat vol1=new ImageDataFloat(input1Vol.getImageData());
		input1Vol.dispose();
		
		ImageData mask=maskVol.getImageData();
		//maskVol.dispose();
		
		Point3i origDim = new Point3i(vol1.getRows(),vol1.getCols(),vol1.getSlices());
		CropParameters params = cropper.findBoundaries(mask, 0.0f, borderSize.getInt());
		ImageData cvol1= cropper.crop(vol1, params);
		
		
		offset.setValue(new Point3i(params.xmin,params.ymin,params.zmin));
		dimNew.setValue(new Point3i(params.getCroppedRows(),params.getCroppedCols(),params.getCroppedSlices()));
		dim.setValue(origDim);
		
		cvol1.setName(cvol1.getName()+roiName.getValue());
		
		// change the border if desired
		float val = 0.0f;
		if (borderType.getValue().equals("min")) {
			val = ImageStatistics.minimum(vol1.toArray3d(),	vol1.getRows(),vol1.getCols(),vol1.getSlices());
		} else if (borderType.getValue().equals("max")) {
			val = ImageStatistics.maximum(vol1.toArray3d(),	vol1.getRows(),vol1.getCols(),vol1.getSlices());
		}
		vol1 = null;
		if (!borderType.getValue().equals("image")) {
			for (int x=0; x<cvol1.getRows(); x++) for (int y=0; y<cvol1.getCols(); y++) for (int z = 0; z<cvol1.getSlices(); z++) for (int t = 0; t<cvol1.getComponents(); t++) {
				if (x<borderSize.getInt()|| x>=cvol1.getRows()-borderSize.getInt() 
					||	y<borderSize.getInt()|| y>=cvol1.getCols()-borderSize.getInt() 
					||	z<borderSize.getInt()|| z>=cvol1.getSlices()-borderSize.getInt() ) {
					cvol1.set(x, y, z, t, val);
				}
			}
		}
		
		
		output1Vol.setValue(cvol1);
		cvol1 = null;
		output1Vol.writeAndFreeNow(this);
		
		if (input2Vol.getImageData()!=null) {
			ImageDataFloat vol2=new ImageDataFloat(input2Vol.getImageData());
			input2Vol.dispose();
			ImageData cvol2 = cropper.crop(vol2, params);
			if (borderType.getValue().equals("min")) {
				val = ImageStatistics.minimum(vol2.toArray3d(),	vol2.getRows(),vol2.getCols(),vol2.getSlices());
			} else if (borderType.getValue().equals("max")) {
				val = ImageStatistics.maximum(vol2.toArray3d(),	vol2.getRows(),vol2.getCols(),vol2.getSlices());
			}
			vol2 = null;
			cvol2.setName(cvol2.getName()+roiName.getValue());
			if (!borderType.getValue().equals("image")) {
				for (int x=0; x<cvol2.getRows(); x++) for (int y=0; y<cvol2.getCols(); y++) for (int z = 0; z<cvol2.getSlices(); z++) for (int t = 0; t<cvol2.getComponents(); t++) {
					if (x<borderSize.getInt()|| x>=cvol2.getRows()-borderSize.getInt() 
						||	y<borderSize.getInt()|| y>=cvol2.getCols()-borderSize.getInt() 
						||	z<borderSize.getInt()|| z>=cvol2.getSlices()-borderSize.getInt() ) {
						cvol2.set(x, y, z, t, val);
					}
				}
			}
			output2Vol.setValue(cvol2);
			cvol2 = null;
			output2Vol.writeAndFreeNow(this);
		}
		
		if (input3Vol.getImageData()!=null) {
			ImageDataFloat vol3=new ImageDataFloat(input3Vol.getImageData());
			input3Vol.dispose();
			ImageData cvol3 = cropper.crop(vol3, params);
			if (borderType.getValue().equals("min")) {
				val = ImageStatistics.minimum(vol3.toArray3d(),	vol3.getRows(),vol3.getCols(),vol3.getSlices());
			} else if (borderType.getValue().equals("max")) {
				val = ImageStatistics.maximum(vol3.toArray3d(),	vol3.getRows(),vol3.getCols(),vol3.getSlices());
			}
			vol3 = null;
			cvol3.setName(cvol3.getName()+roiName.getValue());
			if (!borderType.getValue().equals("image")) {
				for (int x=0; x<cvol3.getRows(); x++) for (int y=0; y<cvol3.getCols(); y++) for (int z = 0; z<cvol3.getSlices(); z++) for (int t = 0; t<cvol3.getComponents(); t++) {
					if (x<borderSize.getInt()|| x>=cvol3.getRows()-borderSize.getInt() 
						||	y<borderSize.getInt()|| y>=cvol3.getCols()-borderSize.getInt() 
						||	z<borderSize.getInt()|| z>=cvol3.getSlices()-borderSize.getInt() ) {
						cvol3.set(x, y, z, t, val);
					}
				}
			}
			output3Vol.setValue(cvol3);
			cvol3 = null;
			output3Vol.writeAndFreeNow(this);
		}
		
		if (input4Vol.getImageData()!=null) {
			ImageDataFloat vol4=new ImageDataFloat(input4Vol.getImageData());
			input4Vol.dispose();
			ImageData cvol4 = cropper.crop(vol4, params);
			if (borderType.getValue().equals("min")) {
				val = ImageStatistics.minimum(vol4.toArray3d(),	vol4.getRows(),vol4.getCols(),vol4.getSlices());
			} else if (borderType.getValue().equals("max")) {
				val = ImageStatistics.maximum(vol4.toArray3d(),	vol4.getRows(),vol4.getCols(),vol4.getSlices());
			}
			vol4 = null;
			cvol4.setName(cvol4.getName()+roiName.getValue());
			if (!borderType.getValue().equals("image")) {
				for (int x=0; x<cvol4.getRows(); x++) for (int y=0; y<cvol4.getCols(); y++) for (int z = 0; z<cvol4.getSlices(); z++) for (int t = 0; t<cvol4.getComponents(); t++) {
					if (x<borderSize.getInt()|| x>=cvol4.getRows()-borderSize.getInt() 
						||	y<borderSize.getInt()|| y>=cvol4.getCols()-borderSize.getInt() 
						||	z<borderSize.getInt()|| z>=cvol4.getSlices()-borderSize.getInt() ) {
						cvol4.set(x, y, z, t, val);
					}
				}
			}
			output4Vol.setValue(cvol4);
			cvol4 = null;
			output4Vol.writeAndFreeNow(this);
		}
	}
}
