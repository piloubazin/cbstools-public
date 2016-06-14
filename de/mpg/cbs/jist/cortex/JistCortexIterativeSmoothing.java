package de.mpg.cbs.jist.cortex;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
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
	
	private ParamDouble fwhmParam;
	
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
		
		inputParams.add(fwhmParam = new ParamDouble("FWHM (mm)", 0.0, 50.0, 5.0));
		
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
		
		ImageDataFloat dataImg = new ImageDataFloat(dataImage.getImageData());
		
		int nx = dataImg.getRows();
		int ny = dataImg.getCols();
		int nz = dataImg.getSlices();
		int nlayers = dataImg.getComponents()-1;
		int nxyz = nx*ny*nz;
		float rx = dataImg.getHeader().getDimResolutions()[0];
		float ry = dataImg.getHeader().getDimResolutions()[1];
		float rz = dataImg.getHeader().getDimResolutions()[2];
		
		float[] data = new float[nxyz];
		float[][][] buffer = dataImg.toArray3d();
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			data[xyz] = buffer[x][y][z];
		}
		buffer = null;
		dataImg.dispose();
		dataImage.dispose();
		
		ImageDataFloat gwbImg = new ImageDataFloat(gwbImage.getImageData());
		ImageDataFloat cgbImg = new ImageDataFloat(cgbImage.getImageData());
		float[][][] gwb = gwbImg.toArray3d();
		float[][][] cgb = cgbImg.toArray3d();
		gwbImg.dispose();
		cgbImg.dispose();
		gwbImage.dispose();
		cgbImage.dispose();

		boolean[] mask = new boolean[nxyz];
		if(maskImage.getImageData() != null) {
			ImageDataUByte maskImg = new ImageDataUByte(maskImage.getImageData());
			byte[][][] bytebuffer = maskImg.toArray3d();
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				
				if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) {
					mask[xyz] = false;
				} else {
					if (gwb[x][y][z]>=(0.0f) && cgb[x][y][z]<=(0.0f)) {
						mask[xyz] = true;
					} else {
						mask[xyz] = false;
					}
				}
			}
			bytebuffer = null;
			maskImg.dispose();
			maskImage.dispose();
		} else {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				mask[xyz] = true;
			}
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) {
				mask[xyz] = false;
			} else {
				if (mask[xyz] && gwb[x][y][z]>=(0.0f) && cgb[x][y][z]<=(0.0f)) {
					mask[xyz] = true;
				} else {
					mask[xyz] = false;
				}
			}
		}
		gwb = null;
		cgb = null;
			
		// 6-connected
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		
		BitSet white = new BitSet(nx*ny*nz);
		BitSet black = new BitSet(nx*ny*nz);
		
		//int count;
		float[] sdata = new float[nxyz];
		double sigma = 0.5;
		float weight = (float) FastMath.exp((double) -(rx*rx)/(2*sigma*sigma));
		float sumweight;
		int iterations = (int) FastMath.round(FastMath.pow(fwhmParam.getValue().doubleValue()/(2*FastMath.sqrt(FastMath.log(4.0))), 2)/sigma);
		System.out.println("\n Number of iterations with sigma = 1: "+iterations);
		System.out.println("\n Neighbour weight = "+weight);
		
		// set voxels in mask as either white or black checks.
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;	
			if (mask[xyz]) {
				sdata[xyz] = data[xyz];
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
					if (mask[xyzn]) {
						sdata[xyz] += weight * sdata[xyzn];
						sumweight += weight;
						//count++;
					}
				}
				sdata[xyz] /= sumweight;
			}
			
			for (int xyz = black.nextSetBit(0); xyz >= 0; xyz = black.nextSetBit(xyz+1)) {
				//count = 1;
				sumweight = 1.0f;
				for (int k=0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (mask[xyzn]) {
						sdata[xyz] += weight * sdata[xyzn];
						sumweight += weight;
						//count++;
					}
				}
				sdata[xyz] /= sumweight;
			}
		}
		
		
        System.out.println("\n done!");
        System.out.println("\n Outputting data");
		
		
        float[][][] smoothData = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			smoothData[x][y][z] = sdata[xyz];
		}
		
		
		ImageDataFloat smoothDataImg = new ImageDataFloat(smoothData);		
		smoothDataImg.setHeader(dataImage.getImageData().getHeader());
		smoothDataImg.setName(dataImage.getImageData().getName()+"_boundary");
		smoothDataImage.setValue(smoothDataImg);
		
	}
}
