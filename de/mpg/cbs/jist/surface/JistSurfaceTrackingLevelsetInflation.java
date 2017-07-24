package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
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
public class JistSurfaceTrackingLevelsetInflation extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume levelsetImage;
	private ParamVolume intens1Image;
	private ParamVolume intens2Image;
	private ParamVolume intens3Image;
	
	private ParamFloat 	smoothingParam;
	private ParamFloat 	precisionParam;
	private ParamInteger 	iterationParam;
	private ParamInteger 	smoothiterParam;
	//private ParamInteger 	projectionParam;
	private ParamOption 	topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamOption 	mappingParam;
	private static final String[] mappingTypes = {"nearest_neighbor", "linear_interpolated", "closest_to_mean", "oversampled"};
	
	private ParamVolume inflatedImage;
	private ParamVolume mappingImage;
	private ParamVolume invmappingImage;
	private ParamVolume debugImage;
	private ParamVolume diffImage;
	private ParamVolume idiffImage;
	private ParamVolume mapped1Image;
	private ParamVolume mapped2Image;
	private ParamVolume mapped3Image;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(levelsetImage = new ParamVolume("Surface level set Image"));
		inputParams.add(intens1Image = new ParamVolume("Intensity Image 1 (opt)"));
		inputParams.add(intens2Image = new ParamVolume("Intensity Image 2 (opt)"));
		inputParams.add(intens3Image = new ParamVolume("Intensity Image 3 (opt)"));
		intens1Image.setMandatory(false);		
		intens2Image.setMandatory(false);		
		intens3Image.setMandatory(false);		

		
		inputParams.add(smoothingParam = new ParamFloat("Smoothing scale (voxels)", 0.0f, 50.0f, 1.5f));
		inputParams.add(smoothiterParam = new ParamInteger("Smoothing steps", 0, 100000, 10));
		inputParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		inputParams.add(precisionParam = new ParamFloat("Precision", 0, 1E9f, 0.001f));
		//inputParams.add(projectionParam = new ParamInteger("Max. projection", 0, 100000, 20));
			
		inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("no");

		inputParams.add(mappingParam = new ParamOption("Inverse mapping", mappingTypes));
		mappingParam.setValue("nearest_neighbor");

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces.devel");
		inputParams.setLabel("Tracking Levelset Inflation");
		inputParams.setName("TrackingLevelsetInflation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Inflates the cortex with a levelset approach and track the result.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(inflatedImage = new ParamVolume("Inflated level set surface",VoxelType.FLOAT));
		outputParams.add(mappingImage = new ParamVolume("Mapping function",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(invmappingImage = new ParamVolume("Inverse Mapping function",VoxelType.FLOAT,-1,-1,-1,-1));
		//outputParams.add(debugImage = new ParamVolume("Debug data",VoxelType.FLOAT,-1,-1,-1,-1));
		//outputParams.add(diffImage = new ParamVolume("Mapped level set",VoxelType.FLOAT));
		//outputParams.add(idiffImage = new ParamVolume("Inverse mapped levelset",VoxelType.FLOAT));
		outputParams.add(mapped1Image = new ParamVolume("Mapped intensity 1",VoxelType.FLOAT));
		outputParams.add(mapped2Image = new ParamVolume("Mapped intensity 2",VoxelType.FLOAT));
		outputParams.add(mapped3Image = new ParamVolume("Mapped intensity 3",VoxelType.FLOAT));
		mapped1Image.setMandatory(false);		
		mapped2Image.setMandatory(false);		
		mapped3Image.setMandatory(false);		
			
		outputParams.setName("inflated images");
		outputParams.setLabel("inflated images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		
		ImageDataFloat	levelsetImg = new ImageDataFloat(levelsetImage.getImageData());
		
		int nx = levelsetImg.getRows();
		int ny = levelsetImg.getCols();
		int nz = levelsetImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = levelsetImg.getHeader().getDimResolutions()[0];
		float ry = levelsetImg.getHeader().getDimResolutions()[1];
		float rz = levelsetImg.getHeader().getDimResolutions()[2];
		
		float[] levelset = new float[nxyz];
		float[][][] buffer = levelsetImg.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			levelset[xyz] = buffer[x][y][z];
		}
		buffer = null;
		levelsetImg = null;
		
		int nimg=0;
		float[][] intensity = new float[3][];
		intensity[0] = null;
		intensity[1] = null;
		intensity[2] = null;
		if (intens1Image.getImageData()!=null) {
			ImageDataFloat intensImg = new ImageDataFloat(intens1Image.getImageData());
			buffer = intensImg.toArray3d();
			intensity[nimg] = new float[nxyz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				intensity[nimg][xyz] = buffer[x][y][z];
			}
			buffer = null;
			intensImg = null;
			nimg++;
		}
		if (intens2Image.getImageData()!=null) {
			ImageDataFloat intensImg = new ImageDataFloat(intens2Image.getImageData());
			buffer = intensImg.toArray3d();
			intensity[nimg] = new float[nxyz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				intensity[nimg][xyz] = buffer[x][y][z];
			}
			buffer = null;
			intensImg = null;
			nimg++;
		}
		if (intens3Image.getImageData()!=null) {
			ImageDataFloat intensImg = new ImageDataFloat(intens3Image.getImageData());
			buffer = intensImg.toArray3d();
			intensity[nimg] = new float[nxyz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				intensity[nimg][xyz] = buffer[x][y][z];
			}
			buffer = null;
			intensImg = null;
			nimg++;
		}
		
		BasicInfo.displayMessage("mask creation...\n");

		// create a mask
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = (levelset[xyz]>=-5.0f && levelset[xyz]<=5.0f);
		}
		// important: increase the data range (useful for better smoothing)
		int delta = 30;
		ObjectMorphology.fastDilateObject(bgmask, nx, ny, nz, delta);
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x==0 || x==nx-1 || y==0 || y==ny-1 || z==0 || z==nz-1) bgmask[xyz] = false;
		}
		
		// create a mapping mask
		boolean[] intensitymask = new boolean[nxyz];
		if (intens1Image.getImageData()!=null) {
			for (int xyz=0;xyz<nxyz;xyz++) {
				intensitymask[xyz] = (intensity[0][xyz]>0.0f);
			} 
		} else {
			intensitymask = null;
		}
		
		BasicInfo.displayMessage("re-build levelset...\n");

			// main algorithm
			CorticalInflationGdm gdm = new CorticalInflationGdm(levelset, nx, ny, nz, rx, ry, rz, 
																	bgmask, 0.4f, 0.4f, topologyParam.getValue(),
																	intensity, intensitymask, nimg); 
			
			/* seems to work for the inflation part (still a few pbs) */
			double basis = 1.0f;
			double scale = smoothingParam.getValue().floatValue();
			if (smoothiterParam.getValue().floatValue()>1) {
				basis = Math.pow(2.0f*smoothingParam.getValue().floatValue(), 1.0f/(smoothiterParam.getValue().floatValue()-1.0f));
				scale = 0.5f;
			}
			for (int t=0;t<smoothiterParam.getValue().intValue();t++) {
				BasicInfo.displayMessage("step "+(t+1)+"\n");
				
				BasicInfo.displayMessage("input smoothing ("+scale+")...\n");
				//gdm.smoothLevelset((float)scale);
				gdm.computeSmoothedLevelset((float)scale, false);
	
				BasicInfo.displayMessage("inflated surface estimation...\n");
				if (mappingParam.getValue().equals("linear_interpolated"))
					gdm.evolveNarrowBandMappingMean(iterationParam.getValue().intValue(), precisionParam.getValue().floatValue());
				else if (mappingParam.getValue().equals("nearest_neighbor"))
					gdm.evolveNarrowBandMappingNearest(iterationParam.getValue().intValue(), precisionParam.getValue().floatValue());
				else if (mappingParam.getValue().equals("closest_to_mean"))
					gdm.evolveNarrowBandMappingClosestMean(iterationParam.getValue().intValue(), precisionParam.getValue().floatValue());
				
				//BasicInfo.displayMessage("copy target...\n");
				//gdm.updateTarget();
				scale *= basis;
			}
			gdm.cleanupForwardMapping();
			gdm.cleanupBackwardMapping();
			
			// output
			String imgname = levelsetImage.getImageData().getName();
			
			ImageDataFloat inflatedData = new ImageDataFloat(gdm.exportLevelset());		
			inflatedData.setHeader(levelsetImage.getImageData().getHeader());
			inflatedData.setName(imgname+"_inf");
			inflatedImage.setValue(inflatedData);
			inflatedData = null;		
			
			ImageDataFloat mappingData = new ImageDataFloat(gdm.exportMappingFD());
			//ImageDataFloat mappingData = new ImageDataFloat(gdm.generateDifferenceMapping());
			mappingData.setHeader(levelsetImage.getImageData().getHeader());
			mappingData.setName(imgname+"_map");
			mappingImage.setValue(mappingData);
			mappingData = null;
			
			ImageDataFloat invmappingData = new ImageDataFloat(gdm.exportMappingBD());
			//ImageDataFloat invmappingData = new ImageDataFloat(gdm.generateInverseDifferenceMapping());
			invmappingData.setHeader(levelsetImage.getImageData().getHeader());
			invmappingData.setName(imgname+"_invmap");
			invmappingImage.setValue(invmappingData);
			invmappingData = null;
			
			/*
			ImageDataFloat debugData = new ImageDataFloat(gdm.generateMappingDistances());
			debugData.setHeader(levelsetImage.getImageData().getHeader());
			debugData.setName(imgname+"_debug");
			debugImage.setValue(debugData);
			debugData = null;
			
			ImageDataFloat diffData = new ImageDataFloat(gdm.generateMappedLevelset());
			diffData.setHeader(levelsetImage.getImageData().getHeader());
			diffData.setName(imgname+"_infmap");
			diffImage.setValue(diffData);
			diffData = null;
			
			ImageDataFloat idiffData = new ImageDataFloat(gdm.generateInverseMappedLevelset());
			idiffData.setHeader(levelsetImage.getImageData().getHeader());
			idiffData.setName(imgname+"_infinv");
			idiffImage.setValue(idiffData);
			idiffData = null;
			*/
			if (intens1Image.getImageData()!=null) {
				ImageDataFloat mappedData=null;
				mappedData= new ImageDataFloat(gdm.generateMappedIntensity(0));
				mappedData.setHeader(intens1Image.getImageData().getHeader());
				mappedData.setName(intens1Image.getImageData().getName()+"_inf");
				mapped1Image.setValue(mappedData);
				mappedData = null;		
			}
			
			if (intens2Image.getImageData()!=null) {
				ImageDataFloat mappedData=null;
				mappedData= new ImageDataFloat(gdm.generateMappedIntensity(1));
				mappedData.setHeader(intens2Image.getImageData().getHeader());
				mappedData.setName(intens2Image.getImageData().getName()+"_inf");
				mapped2Image.setValue(mappedData);
				mappedData = null;		
			}
			
			if (intens3Image.getImageData()!=null) {
				ImageDataFloat mappedData=null;
				mappedData= new ImageDataFloat(gdm.generateMappedIntensity(2));
				mappedData.setHeader(intens3Image.getImageData().getHeader());
				mappedData.setName(intens3Image.getImageData().getName()+"_inf");
				mapped3Image.setValue(mappedData);
				mappedData = null;		
			}
		}	
	//}


}
