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


/*
 * @author Pierre-Louis Bazin
 */
public class JistCortexSmoothData extends ProcessingAlgorithm {

	// jist containers
	
	// input 
	private ParamVolume 	gwblevelsetImage;
	private ParamVolume 	cgblevelsetImage;
	private ParamVolume 	dataImage;
	private ParamVolume 	maskImage;
	
	private ParamFloat		fwhmParam;
	
	// ouput
	private ParamVolume smoothdataImage;
		
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(gwblevelsetImage = new ParamVolume("GWB Levelset Image"));
		inputParams.add(cgblevelsetImage = new ParamVolume("CGB Levelset Image"));
		inputParams.add(dataImage = new ParamVolume("Data Image",VoxelType.FLOAT,-1,-1,-1,-1));
		inputParams.add(maskImage = new ParamVolume("ROI Image (opt)"));
		gwblevelsetImage.setLoadAndSaveOnValidate(false);		
		cgblevelsetImage.setLoadAndSaveOnValidate(false);		
		dataImage.setLoadAndSaveOnValidate(false);	
		maskImage.setMandatory(false);
		
		inputParams.add(fwhmParam = new ParamFloat("Gaussian FWHM (mm)", 0.0f, 50.0f, 5.0f));
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing");
		inputParams.setLabel("Smooth Cortical Data");
		inputParams.setName("SmoothCorticalData");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Smooths data using a Gaussian weighting as a function of the minimum marching distance in the cortex.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(smoothdataImage = new ParamVolume("Smoothed Data Image",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName("Smoothed Cortical Data");
		outputParams.setLabel("Smoothed Cortical Data");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the volume data into 1D arrays
		ImageDataFloat	gwbImg = new ImageDataFloat(gwblevelsetImage.getImageData());
		ImageDataFloat	cgbImg = new ImageDataFloat(cgblevelsetImage.getImageData());
		
		int nx = gwbImg.getRows();
		int ny = gwbImg.getCols();
		int nz = gwbImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = gwbImg.getHeader().getDimResolutions()[0];
		float ry = gwbImg.getHeader().getDimResolutions()[1];
		float rz = gwbImg.getHeader().getDimResolutions()[2];
		int xyz;
		
		float[][][] gwb = gwbImg.toArray3d();
		float[][][] cgb = cgbImg.toArray3d();
		
		// create mask
		byte[][][] roi = null;
		if (maskImage.getImageData()!=null) {
			ImageDataUByte	mskImg = new ImageDataUByte(maskImage.getImageData());
			roi = mskImg.toArray3d();
		}
		boolean[] msk = new boolean[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			xyz = x+nx*y+nx*ny*z;
			if (gwb[x][y][z]>=0 && cgb[x][y][z]<=0) {
				msk[xyz] = true;
				// inside the ROI?
				if (roi!=null && roi[x][y][z]==0) msk[xyz] = false;
			} else {
				msk[xyz] = false;
			}
		}
		gwb=null;
		cgb=null;
		gwbImg.dispose();
		gwblevelsetImage.dispose();
		cgbImg.dispose();
		cgblevelsetImage.dispose();
		
		System.out.println("mask created");
		
		ImageDataFloat	dataImg = new ImageDataFloat(dataImage.getImageData());
		int nd = dataImg.getComponents();
		String dataName = dataImage.getName();
		ImageHeader dataHeader = dataImg.getHeader();
		float[][] data = new float[nd][nxyz];
		if (nd>1) {
			float[][][][] buffer4d = dataImg.toArray4d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				xyz = x+nx*y+nx*ny*z;
				for (int d=0;d<nd;d++) {
					data[d][xyz] = buffer4d[x][y][z][d];
				}
			}
			buffer4d = null;
		} else {
			float[][][] buffer3d = dataImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				xyz = x+nx*y+nx*ny*z;
				data[0][xyz] = buffer3d[x][y][z];
			}
			buffer3d = null;
		}			
		dataImg.dispose();
		dataImage.dispose();
		
		System.out.println("data loaded");
		
		// output variables
		ImageDataFloat sdataImg;
		float[][][][] sdata = new float[nx][ny][nz][nd];
			
		System.out.println("output data created");
		
		// main algorithm
		
		CorticalRegionSmoothing4D regionSmoothing = new CorticalRegionSmoothing4D(data, fwhmParam.getValue().floatValue(), nx, ny, nz, nd, rx, ry, rz, msk);
		
		System.out.println("smoothing loop");
		
		// for debug
		int pts=0;
		for (int x=0;x<nx-1;x++) for (int y=0;y<ny-1;y++) for (int z=0;z<nz-1;z++) {
			xyz = x+nx*y+nx*ny*z;
			if (msk[xyz]) pts++;
		}
		int npt = 0;
		int iter = 0;
		long looptime = System.currentTimeMillis();
		for (int x=0;x<nx-1;x++) for (int y=0;y<ny-1;y++) for (int z=0;z<nz-1;z++) {
			xyz = x+nx*y+nx*ny*z;
			
			if (msk[xyz]) {
				// no need to compute outside the cortex
				npt++;
				if (npt%(pts/100)==0) {
					iter++;
					long newtime = System.currentTimeMillis();
					BasicInfo.displayMessage("iter "+iter+", t="+(newtime-looptime)+", "+npt+" pts\n");
					System.out.println("iter "+iter+", t="+(newtime-looptime)+", "+npt+" pts\n");
					looptime = newtime;
				}
				
				// minimal marching distance from voxel
				//sdata[x][y][z] = regionSmoothing.smoothFromPoint( xyz);
				regionSmoothing.dilateFromPoint(sdata[x][y][z], xyz);
			
			} else {
				for (int d=0;d<nd;d++) {
					sdata[x][y][z][d] = 0.0f;
				}
			}
				
		}
		data = null;
		
		System.out.println("output..");
		
		if (nd==1) {
			// collapse to 3D	
			float[][][] buffer3d = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				buffer3d[x][y][z] = sdata[x][y][z][0];
			}
			sdata = null;
			sdataImg = new ImageDataFloat(buffer3d);
		} else {
			sdataImg = new ImageDataFloat(sdata);
			sdata = null;
		}
		sdataImg.setHeader(dataHeader);
		sdataImg.setName(dataName+ "_smooth");
		smoothdataImage.setValue(sdataImg);
		
		// clean up
		sdataImg = null;
		regionSmoothing.finalize();
		
	}


}
