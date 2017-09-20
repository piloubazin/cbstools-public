package de.mpg.cbs.jist.cortex;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;
import org.apache.commons.math3.util.FastMath;

/*
 * @author Christine Lucas Tardif (ctardif@cbs.mpg.de)
 *
 */

/*
 * bug to fix: if the distance between the gwb and cgb is smaller than a voxel, there will be no output in this area.
 * 
 */
public class JistCortexIterativeSmoothing extends ProcessingAlgorithm{
	private ParamVolume dataImage;
	private ParamVolume gwbImage;
	private ParamVolume cgbImage;
	private ParamVolume maskImage;
	
	private ParamFloat fwhmParam;
	
	private ParamVolume smoothDataImage;
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(dataImage = new ParamVolume("Cortical Data"));
		inputParams.add(gwbImage = new ParamVolume("Grey-White or inner boundary"));
		inputParams.add(cgbImage = new ParamVolume("CSF-Grey or outer boundary"));
		inputParams.add(maskImage = new ParamVolume("Cortical Mask"));
		maskImage.setMandatory(false);
		
		inputParams.add(fwhmParam = new ParamFloat("FWHM (mm)", 0.0f, 50.0f, 5.0f));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing.devel");
		inputParams.setLabel("Iterative Cortical Smoothing");
		inputParams.setName("IterativeCorticalSmoothing");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Smooths cortical data using a heat kernel approach: an iterative algorithm that uses a weighted average of the data at each voxel with its nearest neighbours. May not work well in areas of thin cortex.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
		
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(smoothDataImage = new ParamVolume("Smoothed Data",VoxelType.FLOAT));
		
		outputParams.setName("SmoothedData");
		outputParams.setLabel("Smoothed Data");
	}

	
	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		
		// import the image data
		System.out.println("\n Import data");
		
		ImageHeader header = dataImage.getImageData().getHeader();
		String name = dataImage.getImageData().getName();
		
		ImageDataFloat dataImg = new ImageDataFloat(dataImage.getImageData());
		
		int nx = dataImg.getRows();
		int ny = dataImg.getCols();
		int nz = dataImg.getSlices();
		int nt = dataImg.getComponents();
		int nxyz = nx*ny*nz;
		float rx = dataImg.getHeader().getDimResolutions()[0];
		float ry = dataImg.getHeader().getDimResolutions()[1];
		float rz = dataImg.getHeader().getDimResolutions()[2];
		
		
		float[] data = null;
		if (nt>1) data = Interface.getFloatImage4D(dataImage);
		else data = Interface.getFloatImage3D(dataImage);
		
		dataImg.dispose();
		dataImage.dispose();
		
		float[] gwb = Interface.getFloatImage3D(gwbImage);
		float[] cgb = Interface.getFloatImage3D(cgbImage);
		
		gwbImage.dispose();
		cgbImage.dispose();

		byte[] mask = null;
		if (Interface.isValid(maskImage)) {
			mask = Interface.getUByteImage3D(maskImage);
		} else {
			mask = new byte[nxyz];
			for (int xyz=0;xyz<nxyz;xyz++) mask[xyz] = 1;
		}
		
		// remove image boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) {
				mask[xyz] = 0;
			} else {
				if (mask[xyz]>0 && gwb[xyz]>=(0.0f) && cgb[xyz]<=(0.0f)) {
					mask[xyz] = 1;
				} else {
					mask[xyz] = 0;
				}
			}
		}
		gwb = null;
		cgb = null;
			
		// 6-connected
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		
		BitSet white = new BitSet(nxyz);
		BitSet black = new BitSet(nxyz);
		
		//int count;
		float[] sdata = new float[nxyz];
		double sigma = 0.5f;
		float weight = (float) FastMath.exp((double) -(rx*rx)/(2*sigma*sigma));
		float sumweight;
		int iterations = (int) FastMath.round(FastMath.pow(fwhmParam.getValue().doubleValue()/(2*FastMath.sqrt(FastMath.log(4.0f))), 2)/sigma);
		System.out.println("\n Number of iterations with sigma = 1: "+iterations);
		System.out.println("\n Neighbour weight = "+weight);
		
		// set voxels in mask as either white or black checks.
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;	
			if (mask[xyz]>0) {
				for (int t=0;t<nt;t++) sdata[xyz+nxyz*t] = data[xyz+nxyz*t];
				if (xyz%2 == 0) {
					white.set(xyz);
				} else {
					black.set(xyz);
				}
			}
		}
		
		
		for (int itr=0; itr<iterations; itr++) {
			
			for (int xyz = white.nextSetBit(0); xyz >= 0; xyz = white.nextSetBit(xyz+1)) {
				//count = 1;
				sumweight = 1.0f;
				for (int k=0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (mask[xyzn]>0) {
						for (int t=0;t<nt;t++) sdata[xyz+nxyz*t] += weight * sdata[xyzn+nxyz*t];
						sumweight += weight;
						//count++;
					}
				}
				for (int t=0;t<nt;t++) sdata[xyz+nxyz*t] /= sumweight;
			}
			
			for (int xyz = black.nextSetBit(0); xyz >= 0; xyz = black.nextSetBit(xyz+1)) {
				//count = 1;
				sumweight = 1.0f;
				for (int k=0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (mask[xyzn]>0) {
						for (int t=0;t<nt;t++) sdata[xyz+nxyz*t] += weight * sdata[xyzn+nxyz*t];
						sumweight += weight;
						//count++;
					}
				}
				for (int t=0;t<nt;t++) sdata[xyz+nxyz*t] /= sumweight;
			}
		}
		
		
        System.out.println("\n done!");
        System.out.println("\n Outputting data");
		
		if (nt>1) Interface.setFloatImage4D(sdata, new int[]{nx,ny,nz}, nt, smoothDataImage, name+"_boundary", header);
		else Interface.setFloatImage3D(sdata, new int[]{nx,ny,nz}, smoothDataImage, name+"_boundary", header);
	}
}
