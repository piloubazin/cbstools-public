package de.mpg.cbs.jist.segmentation;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.algorithms.PrinceGroupAuthors;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSegmentationFitMgdmEdges extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume labelImage;
	private ParamVolume mgdmImage;
	
	private ParamVolume intensityImage;
	private ParamVolume maskImage;
	
	private	ParamFloat 	balloonParam;
	private	ParamFloat 	edgeParam;
	private ParamFloat 	curvParam;
	private ParamFloat 	spatialParam;
	
	private ParamInteger 	iterationParam;
	private ParamFloat 	changeParam;
	
	private ParamOption 	topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume newlabelImage;
	private ParamVolume newmgdmImage;
	private ParamVolume gradientImage;
	//private ParamVolume debugImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(labelImage = new ParamVolume("Segmentation Image"));
		imageParams.add(mgdmImage = new ParamVolume("Levelset Boundary Image"));
		imageParams.add(intensityImage = new ParamVolume("New Intensity Image"));
		imageParams.add(maskImage = new ParamVolume("New Intensity Mask (opt)"));
		maskImage.setMandatory(false);
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		mainParams.add(edgeParam = new ParamFloat("Image edge weight", -1E10f, 1E10f, 0.5f));
		mainParams.add(balloonParam = new ParamFloat("Original boundary weight", -1E10f, 1E10f, 0.1f));
		mainParams.add(curvParam = new ParamFloat("Curvature weight", -1E10f, 1E10f, 0.1f));
		mainParams.add(spatialParam = new ParamFloat("Spatial extent (mm)", -1E10f, 1E10f, 2.0f));
		
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		mainParams.add(changeParam = new ParamFloat("Min change", 0, 1, 0.001f));
			
		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		inputParams.add(mainParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation");
		inputParams.setLabel("Fit Segmentation to Edges");
		inputParams.setName("FitSegementationToEdges");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Refines a MGDM-based segmentation based on the edges of a new image"
							+"\n (e.g. a slab of high-resolution data).");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(newlabelImage = new ParamVolume("Refined Labels",VoxelType.BYTE));
		outputParams.add(newmgdmImage = new ParamVolume("Refined Boundaries",VoxelType.FLOAT));
		outputParams.add(gradientImage = new ParamVolume("Image Gradient",VoxelType.FLOAT));
		//outputParams.add(debugImage = new ParamVolume("Debug",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName("refined images");
		outputParams.setLabel("refined images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataByte labelImg = new ImageDataByte(labelImage.getImageData());
		int nx = labelImg.getRows();
		int ny = labelImg.getCols();
		int nz = labelImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = labelImg.getHeader().getDimResolutions()[0];
		float ry = labelImg.getHeader().getDimResolutions()[1];
		float rz = labelImg.getHeader().getDimResolutions()[2];
		byte[][][] bytebuffer = labelImg.toArray3d();
		byte[] label = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			label[xyz] = bytebuffer[x][y][z];
		}
		bytebuffer = null;
		
		ImageDataFloat mgdmImg = new ImageDataFloat(mgdmImage.getImageData());
		float[][][] buffer = mgdmImg.toArray3d();
		float[] mgdm = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			mgdm[xyz] = buffer[x][y][z];
		}
		mgdmImg = null;
		buffer = null;
		
		ImageDataFloat intensityImg = new ImageDataFloat(intensityImage.getImageData());
		buffer = intensityImg.toArray3d();
		float[] intensity = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			intensity[xyz] = buffer[x][y][z];
		}
		buffer = null;
				
		boolean[] mask = new boolean[nxyz];
		if (maskImage.getImageData()!=null) {
			ImageDataByte maskImg = new ImageDataByte(maskImage.getImageData());
			bytebuffer = maskImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				mask[xyz] = (bytebuffer[x][y][z]!=0);
			}
			maskImg = null;
			bytebuffer = null;
		} else {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				mask[xyz] = (intensity[xyz]!=0);
			}
		}
		
		// spatial st.dev
		float stdev = 0.33f*spatialParam.getValue().floatValue()/rx;
		BasicInfo.displayMessage("spatial stdev: "+stdev+" voxels \n");
				
		// MGDM evolution		
		MgdmSegmentationRefinement refine = new MgdmSegmentationRefinement(intensity, mask, nx, ny, nz, rx, ry, rz,
																			label, mgdm, 4,																
																			balloonParam.getValue().floatValue(), 
																			curvParam.getValue().floatValue(), 
																			edgeParam.getValue().floatValue(), 
																			stdev,
																			topologyParam.getValue());
		
		if (iterationParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			refine.evolveNarrowBand(iterationParam.getValue().intValue(),changeParam.getValue().floatValue());
		}
		
		byte[][][] seg = new byte[nx][ny][nz];
		float[] img = refine.exportMap(refine.labelSegmentation());
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[x][y][z] = (byte)img[xyz];
		}
		img = null;
		ImageDataByte segData = new ImageDataByte(seg);	
		seg = null;
		segData.setHeader(intensityImg.getHeader());
		segData.setName(intensityImg.getName()+"_reseg");
		newlabelImage.setValue(segData);
		segData = null;
		BasicInfo.displayMessage("segmentation");
		
		float[][][] lvl = new float[nx][ny][nz];
		float[] fcn = refine.getFunctions()[0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			lvl[x][y][z] = fcn[xyz];
		}
		ImageDataFloat lvlData = new ImageDataFloat(lvl);	
		lvl = null;
		lvlData.setHeader(intensityImg.getHeader());
		lvlData.setName(intensityImg.getName()+"_remgdm");
		newmgdmImage.setValue(lvlData);
		lvlData = null;
		BasicInfo.displayMessage(".. boundaries");
		
		float[][][] grad = refine.exportEdgeNorm();
		ImageDataFloat gradData = new ImageDataFloat(grad);	
		grad = null;
		gradData.setHeader(intensityImg.getHeader());
		gradData.setName(intensityImg.getName()+"_edges");
		gradientImage.setValue(gradData);
		gradData = null;
		BasicInfo.displayMessage(".. edges");
		
		/*
		float[][][][] forces = refine.exportForces();
		ImageDataFloat forceData = new ImageDataFloat(forces);	
		forces = null;
		forceData.setHeader(intensityImg.getHeader());
		forceData.setName(intensityImg.getName()+"_forces");
		debugImage.setValue(forceData);
		forceData = null;
		BasicInfo.displayMessage(".. forces");
		*/
	}

}
