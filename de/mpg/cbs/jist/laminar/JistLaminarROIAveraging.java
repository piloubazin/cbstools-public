package de.mpg.cbs.jist.laminar;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
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

import java.io.File;
import java.io.IOException;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Juliane Dinse
 */
public class JistLaminarROIAveraging extends ProcessingAlgorithm {

	// jist containers
	// input
	private ParamVolume 	roiImage;
	private ParamVolume 	profSamples;
	private ParamString		roiName;
	private ParamVolume 	maskImage;
	
	// output
	private ParamFile		textFile;
	
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(profSamples = new ParamVolume("Intensity Profile Image",null,-1,-1,-1,-1));
		inputParams.add(roiImage = new ParamVolume("ROI Mask",null,-1,-1,-1,-1));
		inputParams.add(roiName = new ParamString("ROI Name"));
		inputParams.add(maskImage = new ParamVolume("Mask Image (opt, 3D or 4D)",null,-1,-1,-1,-1));
		maskImage.setMandatory(false);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis");
		inputParams.setLabel("Profile ROI Averaging");
		inputParams.setName("ProfileROIAveraging");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Juliane Dinse", "dinse@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute an average profile over a given ROI.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(textFile = new ParamFile("ROI Average"));
		
		outputParams.setName("ProfileROIAveraging");
		outputParams.setLabel("Profile ROI Averaging");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		
		ImageDataUByte	roiImageImg = new ImageDataUByte(roiImage.getImageData() );
		ImageDataFloat  profSamplesImg = new ImageDataFloat(profSamples.getImageData());

		
		int nx = profSamplesImg.getRows();
		int ny = profSamplesImg.getCols();
		int nz = profSamplesImg.getSlices();
		int nlayers = profSamplesImg.getComponents()-1;
		float rx = profSamplesImg.getHeader().getDimResolutions()[0];
		float ry = profSamplesImg.getHeader().getDimResolutions()[1];
		float rz = profSamplesImg.getHeader().getDimResolutions()[2];
		
		byte[][][] roiMask = roiImageImg.toArray3d();
		float[][][][] intensity = profSamplesImg.toArray4d();

		boolean[][][] ctxmask = new boolean[nx][ny][nz];
		boolean[][][][] layermask = null;
		if (maskImage.getImageData()==null) {
			// create a mask for regions of value 0 on center image (kinda arbitrary..)
			for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {	
				ctxmask[x][y][z] = (intensity[x][y][z][nlayers/2]!=0.0f);
			}
		} else {
			ImageDataUByte	maskImg = new ImageDataUByte(maskImage.getImageData());
			if (maskImg.getComponents()==1) {
				// 3D mask : replace ctxmask
				byte[][][] bytebuffer = maskImg.toArray3d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					ctxmask[x][y][z] = (bytebuffer[x][y][z]!=0);
				}
				bytebuffer = null;
				maskImg = null;
			} else {
				// 4D mask : create layer mask and replace ctxmask
				layermask = new boolean[nx][ny][nz][nlayers+1];
				byte[][][][] bytebuffer = maskImg.toArray4d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					ctxmask[x][y][z] = false;
					for (int l=0;l<=nlayers;l++) {
						layermask[x][y][z][l] = (bytebuffer[x][y][z][l]!=0);
						if (layermask[x][y][z][l]) ctxmask[x][y][z] = true;
					}
				}
				bytebuffer = null;
				maskImg = null;
			}
		}
		
		float[] average = new float[nlayers+1];
		int[] nrVoxelperLayer = new int[nlayers+1];
		float[] stddeviation = new float[nlayers+1];

		// iterating through all layers
		for(int t = 0; t <= nlayers; t++){
			double totalSumInLayer = 0.0f;
			double nrOfVoxelsInLayer = 0.0f;

			// get each single 3D volume
			for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
				if (ctxmask[x][y][z] && roiMask[x][y][z]!=0) {
					if (layermask==null || layermask[x][y][z][t]) {
						totalSumInLayer += intensity[x][y][z][t];
						nrOfVoxelsInLayer++;
					}
				}
			}

			average[t] = (float)(totalSumInLayer/nrOfVoxelsInLayer);
			nrVoxelperLayer[t] = (int)nrOfVoxelsInLayer;

			double stdvaluetotal = 0.0f;
			for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
				if (ctxmask[x][y][z] && roiMask[x][y][z]!=0) {
					if (layermask==null || layermask[x][y][z][t]) {
						float nvalue = intensity[x][y][z][t]-average[t];
						stdvaluetotal += nvalue*nvalue;
					}
				}
			}
			stddeviation[t] = (float)FastMath.sqrt(stdvaluetotal/nrVoxelperLayer[t]);
		}
		intensity = null;
		roiMask = null;

		String imgname = profSamples.getImageData().getName();

		File resultFile = null;

		String filename = imgname+"_roiavg_"+roiName.getValue()+".txt";
		File destdir = null;
		try{
			destdir = new File(this.getOutputDirectory().getCanonicalFile()+File.separator+this.getAlgorithmName());
			if(!destdir.isDirectory()){
				(new File(destdir.getCanonicalPath())).mkdir();
			}
			resultFile = new File(destdir+File.separator+filename);


		}catch(IOException e){
			e.printStackTrace();
		}


		FileWriter fw;
		try {
			fw = new FileWriter(resultFile);
		} catch (IOException error) {
			return;
		}
		PrintWriter pw = new PrintWriter( fw );
		pw.write("mean intensity per layer : standard deviation : nr of voxels per layer");
		pw.printf("\n");
		// open the file to write
		try {

			for(int i = 0; i < average.length; i++){

				float val = average[i];
				
				float stddev = stddeviation[i];
				int nrVoxLay = nrVoxelperLayer[i];

				pw.printf(" %.2f0f \t  %.8f \t ", val, stddev );
				pw.print(nrVoxLay + "\n");

			}
		}
		catch (Exception e) {
			System.out.println(e.getMessage());
		}	


		pw.close();
		try {
			fw.close();
		} catch (IOException error) {
			return;
		}

		fw = null;
		pw = null;
		stddeviation = null;
		nrVoxelperLayer = null;	
		average = null;

		textFile.setValue(resultFile);
		
	}

}
