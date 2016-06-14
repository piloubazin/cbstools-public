package de.mpg.cbs.jist.laminar;

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


/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarProfileSampling extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume layersImage;
	private ParamVolume intensityImage;
	private ParamVolume maskImage;
		
	private ParamVolume mappedImage;
	private ParamVolume mappedmaskImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(layersImage = new ParamVolume("Profile Surface Image",null,-1,-1,-1,-1));
		imageParams.add(intensityImage = new ParamVolume("Intensity Image",null,-1,-1,-1,-1));
		imageParams.add(maskImage = new ParamVolume("Cortex Mask (opt)",null,-1,-1,-1,-1));
		maskImage.setMandatory(false);
		
		inputParams.add(imageParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis");
		inputParams.setLabel("Profile Sampling");
		inputParams.setName("ProfileSampling");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Juliane Dinse", "dinse@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Sample some intensity image along a cortical profile across layer surfaces.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(mappedImage = new ParamVolume("Profile-mapped Intensity Image",null,-1,-1,-1,-1));
		
		outputParams.add(mappedmaskImage = new ParamVolume("Profile 4D Mask",null,-1,-1,-1,-1));
		
		outputParams.setName("layers images");
		outputParams.setLabel("layers images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays : TO DO
		
		ImageDataFloat	layersImg = new ImageDataFloat(layersImage.getImageData());
		ImageDataFloat	intensImg = new ImageDataFloat(intensityImage.getImageData());
		
		int nx = layersImg.getRows();
		int ny = layersImg.getCols();
		int nz = layersImg.getSlices();
		int nlayers = layersImg.getComponents()-1;
		int nxyz = nx*ny*nz;
		float rx = layersImg.getHeader().getDimResolutions()[0];
		float ry = layersImg.getHeader().getDimResolutions()[1];
		float rz = layersImg.getHeader().getDimResolutions()[2];
		
		float[][] layers = new float[nlayers+1][nxyz];
		float[][][][] buffer4 = layersImg.toArray4d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			layers[l][xyz] = buffer4[x][y][z][l];
		}
		buffer4 = null;
		layersImg = null;
		
		float[] intensity = new float[nxyz];
		float[][][] buffer3 = intensImg.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			intensity[xyz] = buffer3[x][y][z];
		}
		buffer3 = null;
		intensImg = null;
		
		// create a mask for all the regions outside of the area where layer 1 is > 0 and layer 2 is < 0
		boolean[] ctxmask = new boolean[nxyz];
		if (maskImage.getImageData()!=null) {
			ImageDataUByte	maskImg = new ImageDataUByte(maskImage.getImageData());
			byte[][][] bufferbyte = maskImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				ctxmask[xyz] = (layers[0][xyz]>=0.0 && layers[nlayers][xyz]<=0.0 && bufferbyte[x][y][z]>0);
			}
			bufferbyte = null;
			maskImg = null;
		} else {
			for (int xyz=0;xyz<nxyz;xyz++) {
				ctxmask[xyz] = (layers[0][xyz]>=0.0 && layers[nlayers][xyz]<=0.0);
			}
		}
				
		// main algorithm
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		float maskval = 1e13f;
		float[][][][] mapping = new float[nx][ny][nz][nlayers+1];
		byte[][][][] mappingmask = new byte[nx][ny][nz][nlayers+1];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				profile.computeTrajectory(layers, x, y, z);
				
				for (int l=0;l<=nlayers;l++) {
					// interpolate the contrast
					mapping[x][y][z][l] = ImageInterpolation.linearInterpolation(intensity, ctxmask, maskval, 
																					profile.getPt(l)[X], profile.getPt(l)[Y], profile.getPt(l)[Z], 
																					nx, ny, nz);
					if (mapping[x][y][z][l]==maskval) {
						mappingmask[x][y][z][l] = (byte)0;
						mapping[x][y][z][l] = 0.0f;
					} else {
						mappingmask[x][y][z][l] = (byte)1;
					}
				}
			}
		}
		
		// output
		String imgname = intensityImage.getImageData().getName();
		
		ImageDataFloat mapData = new ImageDataFloat(mapping);		
		mapData.setHeader(layersImage.getImageData().getHeader());
		mapData.setName(imgname+"_profiles");
		mappedImage.setValue(mapData);
		mapData = null;
		mapping = null;
		
		ImageDataUByte mapmaskData = new ImageDataUByte(mappingmask);		
		mapmaskData.setHeader(layersImage.getImageData().getHeader());
		mapmaskData.setName(imgname+"_4dmask");
		mappedmaskImage.setValue(mapmaskData);
		mapmaskData = null;
		mappingmask = null;
	}


}
