package de.mpg.cbs.jist.cortex;

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
public class JistCortexMyelinatedThickness extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume intensImage;
	private ParamVolume wmImage;
	private ParamVolume gmImage;
	private ParamVolume subctxImage;
	private ParamVolume pwmImage;
	
	private ParamBoolean	adjustwmParam;
	private ParamFloat 	offsetParam;
	private ParamBoolean	adjustgmParam;
	
	private ParamBoolean	fcmwithwmParam;
	private ParamFloat 	fcmsmoothParam;
	private ParamInteger 	fcmiterParam;
	private ParamFloat 	fcmdiffParam;
	
	private	ParamFloat 	gdmballoonParam;
	private ParamFloat 	gdmcurvParam;
	private ParamInteger 	gdmiterParam;
	private ParamFloat 	gdmdiffParam;
	private ParamOption 	topologyParam;
	private ParamFloat 	maxdistParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume cortexImage;
	private ParamVolume boundaryImage;
	private ParamVolume ratioImage;
	private ParamVolume thickImage;
	
	private static final boolean debug = true;
	private ParamVolume rfcmImage;
	private ParamVolume innerImage;
	private ParamVolume outerImage;
	private ParamVolume memsImage;
	
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(intensImage = new ParamVolume("Intensity Image"));
		imageParams.add(wmImage = new ParamVolume("WM Boundary Surface (levelset)"));
		imageParams.add(gmImage = new ParamVolume("GM Boundary Surface (levelset)"));
		imageParams.add(subctxImage = new ParamVolume("Subcortex Mask (GM nuclei,ventricles)"));
		imageParams.add(pwmImage = new ParamVolume("WM Probability Image (opt)"));
		pwmImage.setMandatory(false);
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		mainParams.add(adjustwmParam = new ParamBoolean("Adjust WM boundary", false));
		mainParams.add(offsetParam = new ParamFloat("WM probability boundary", 0.0f, 1.0f, 0.33f));
		mainParams.add(adjustgmParam = new ParamBoolean("Adjust outer GM boundary", false));
		
		mainParams.add(fcmwithwmParam = new ParamBoolean("Include WM in clustering", true));
		mainParams.add(fcmsmoothParam = new ParamFloat("Fantasm regularization weight", -1E10f, 1E10f, 0.1f));
		mainParams.add(fcmdiffParam = new ParamFloat("Fantasm min difference", -1E10f, 1E10f, 0.01f));
		mainParams.add(fcmiterParam = new ParamInteger("Fantasm max iterations", 0, 100000, 50));
		
		mainParams.add(gdmballoonParam = new ParamFloat("GDM balloon weight", -1E10f, 1E10f, 0.4f));
		mainParams.add(gdmcurvParam = new ParamFloat("GDM curvature weight", -1E10f, 1E10f, 0.05f));
		mainParams.add(gdmdiffParam = new ParamFloat("GDM min difference", -1E10f, 1E10f, 0.0001f));
		mainParams.add(gdmiterParam = new ParamInteger("GDM max iterations", 0, 100000, 500));
			
		mainParams.add(topologyParam = new ParamOption("GDM Topology", topoTypes));
		topologyParam.setValue("wcs");

		mainParams.add(maxdistParam = new ParamFloat("Max distance outside cortex (mm)", -1E10f, 1E10f, 1.0f));
		
		inputParams.add(mainParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing.devel");
		inputParams.setLabel("Myelinated Thickness Estimation");
		inputParams.setName("MyelinatedThicknessEstimation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Nicholas Bock", "bockn@mcmaster.ca","http://www.science.mcmaster.ca/medphys/faculty/45-nicholas-bock.html"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/staff/bazin-11500"));
		info.setAffiliation("McMaster university & Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Estimates the relative thickness of highly myelinated cortex in T1 and T1-weighted images");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(cortexImage = new ParamVolume("Myelinated cortex segmentation",VoxelType.UBYTE));
		outputParams.add(boundaryImage = new ParamVolume("Myelinated cortex levelset surface",VoxelType.FLOAT));
		outputParams.add(ratioImage = new ParamVolume("Myelinated cortex ratio",VoxelType.FLOAT));
		outputParams.add(thickImage = new ParamVolume("Myelinated cortex thickness",VoxelType.FLOAT));
		if (debug) {
			outputParams.add(rfcmImage = new ParamVolume("Fantasm Classification",VoxelType.UBYTE));
			outputParams.add(memsImage = new ParamVolume("Fantasm Memberships",VoxelType.FLOAT,-1,-1,-1,-1));
			outputParams.add(innerImage = new ParamVolume("Adjusted inner levelset surface",VoxelType.FLOAT));
			outputParams.add(outerImage = new ParamVolume("Adjusted outer levelset surface",VoxelType.FLOAT));
		}
			
		outputParams.setName("myelin thickness images");
		outputParams.setLabel("myelin thickness images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataFloat wmImg = new ImageDataFloat(wmImage.getImageData());
		int nx = wmImg.getRows();
		int ny = wmImg.getCols();
		int nz = wmImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = wmImg.getHeader().getDimResolutions()[0];
		float ry = wmImg.getHeader().getDimResolutions()[1];
		float rz = wmImg.getHeader().getDimResolutions()[2];
		
		float[][][] wm = wmImg.toArray3d();
		wmImg = null;
		
		ImageDataFloat gmImg = new ImageDataFloat(gmImage.getImageData());
		float[][][] gm = gmImg.toArray3d();
		gmImg = null;
				
		ImageDataFloat intensImg = new ImageDataFloat(intensImage.getImageData());
		float[][][] intens = intensImg.toArray3d();
		intensImg = null;
		
		ImageDataUByte subctxImg = new ImageDataUByte(subctxImage.getImageData());
		byte[][][] subctx = subctxImg.toArray3d();
		subctxImg = null;
		
		float outside = 0.05f; // allows the segmentation to go slightly out of the original one
		
		float[] wmadjusted = null;
		if (adjustwmParam.getValue().booleanValue() && pwmImage.getImageData()!=null) {
			ImageDataFloat pwmImg = new ImageDataFloat(pwmImage.getImageData());
			float[][][] pwm = pwmImg.toArray3d();
			pwmImg = null;
			// 1. evolve WM boundary (optional?)
			BasicInfo.displayMessage("refining WM boundary\n");
			float[] forces = new float[nxyz];
			float[] levelset = new float[nxyz];
			float factor = 1.0f/offsetParam.getValue().floatValue();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				levelset[xyz] = wm[x][y][z];
				if (gm[x][y][z]>outside) {
					// outside the GM : pull back
					forces[xyz] = -1.0f;
				} else if (subctx[x][y][z]>0) {
					// inside the subcortex : push out
					forces[xyz] = +1.0f;
				} else {
					// balloon = bounded( 1/offset * WM - 1 , -1, 1)
					forces[xyz] = Numerics.bounded(factor*pwm[x][y][z]-1.0f, -1.0f, 1.0f);
				}
			}
			wm = null;
			
			Gdm3d gdm = new Gdm3d(levelset, 20.0f/rx, nx, ny, nz, rx, ry, rz, null,
									forces, 0.0f, gdmballoonParam.getValue().floatValue(), 
									gdmcurvParam.getValue().floatValue(), topologyParam.getValue());
			
			if (gdmiterParam.getValue().intValue()>0) {
				BasicInfo.displayMessage("GDM adjustment...\n");
				
				gdm.evolveNarrowBand(gdmiterParam.getValue().intValue(), gdmdiffParam.getValue().floatValue());
			}
			
			wmadjusted = gdm.exportLevelset();
		} else {
			wmadjusted = new float[nxyz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				wmadjusted[xyz] = wm[x][y][z];
			}
			wm = null;
		}
		
		// 2. run Fantasm on cortex data
		BasicInfo.displayMessage("Cortical depth clustering...\n");
		boolean withwm = fcmwithwmParam.getValue().booleanValue();
		if (withwm) BasicInfo.displayMessage("(clustering WM, GMm, GM)\n");
		else BasicInfo.displayMessage("(clustering GMm, GM only)\n");
		
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			// include the cortex, possibly add region inside the WM with proba > threshold
			if (wmadjusted[xyz]>=0 && gm[x][y][z]<=outside) mask[x][y][z] = true;
			else if (withwm && subctx[x][y][z]==0 && wmadjusted[xyz]<=0) mask[x][y][z] = true;
			else mask[x][y][z] = false;
		}
		BasicRFCM rfcm = null;
		if (withwm) rfcm = new BasicRFCM(intens, mask, nx,ny,nz,3,fcmsmoothParam.getValue().floatValue());
		else rfcm = new BasicRFCM(intens, mask, nx,ny,nz,2,fcmsmoothParam.getValue().floatValue());
		
		rfcm.initCentroidsRange();
		int nt = fcmiterParam.getValue().intValue();
		float maxdiff = fcmdiffParam.getValue().floatValue();
		for (int t=0;t<nt;t++) {
			float diff = rfcm.computeMemberships();
			rfcm.computeCentroids();
			
			BasicInfo.displayMessage("Iteration "+(t+1)+" ("+diff+")\n");
			if (diff<maxdiff) t=nt;
		}
		// re-compute memberships for the entire image (to avoid empty regions at the WM/GM interface)
		rfcm.computeMembershipsEverywhere();
		
		// 3. use classification result to generate a myelinated boundary
		float dist1 = 0.0f, dist2 = 0.0f, dist3 = 0.0f;
		float sum1 = 0.0f, sum2 = 0.0f, sum3 = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			dist1 += -rfcm.getMemberships()[x][y][z][0]*gm[x][y][z];
			dist2 += -rfcm.getMemberships()[x][y][z][1]*gm[x][y][z];
			if (withwm) dist3 += -rfcm.getMemberships()[x][y][z][2]*gm[x][y][z];
			sum1 += rfcm.getMemberships()[x][y][z][0];
			sum2 += rfcm.getMemberships()[x][y][z][1];
			if (withwm) sum3 += rfcm.getMemberships()[x][y][z][2];
		}
		if (withwm) BasicInfo.displayMessage("classes distance to pial surface: 1 ("+(dist1/sum1)+"), 2 ("+(dist2/sum2)+"), 3 ("+(dist3/sum3)+")\n");
		else BasicInfo.displayMessage("classes distance to pial surface: 1 ("+(dist1/sum1)+"), 2 ("+(dist2/sum2)+")\n");
			
		int inner = -1, central = -1, outer = -1;
		if (withwm) {
			inner = Numerics.argmax(dist1/sum1,dist2/sum2,dist3/sum3);
			central = Numerics.argsec(dist1/sum1,dist2/sum2,dist3/sum3);
			outer = Numerics.argmin(dist1/sum1,dist2/sum2,dist3/sum3);
		} else {
			central = Numerics.argmax(dist1/sum1,dist2/sum2);
			outer = Numerics.argmin(dist1/sum1,dist2/sum2);
		}
		
		// 4. adjust all boundaries (WM if needed, GMm, GM if needed)
		float[] levelset = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			levelset[xyz] = wmadjusted[xyz];
		}
		float[] forces = new float[nxyz];
		if (withwm) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (gm[x][y][z]>outside) {
					// outside the GM : pull back
					forces[xyz] = -1.0f;
				} else if (subctx[x][y][z]>0) {
					// inside the subcortex : push out
					forces[xyz] = +1.0f;
				} else {
					// balloon = bounded( 2 * mems - 1 , -1, 1)
					//forces[xyz] = Numerics.bounded(2.0f*rfcm.getMemberships()[x][y][z][inner]-1.0f, -1.0f, 1.0f);
					// or local max membership (nicer conceptually)
					forces[xyz] = Numerics.bounded(rfcm.getMemberships()[x][y][z][inner]
													-Numerics.max(rfcm.getMemberships()[x][y][z][central],
																	rfcm.getMemberships()[x][y][z][outer]), -1.0f, 1.0f);
					
				}
			}

			Gdm3d gdm = new Gdm3d(levelset, 20.0f/rx, nx, ny, nz, rx, ry, rz, null,
									forces, 0.0f, gdmballoonParam.getValue().floatValue(), 
									gdmcurvParam.getValue().floatValue(), topologyParam.getValue());
			
			if (gdmiterParam.getValue().intValue()>0) {
				BasicInfo.displayMessage("GDM final WM adjustment...\n");
				
				gdm.evolveNarrowBand(gdmiterParam.getValue().intValue(), gdmdiffParam.getValue().floatValue());
			}
			levelset = gdm.exportLevelset();
			wmadjusted = gdm.exportLevelset();
		}
		
		// new forces for myelinated boundary
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (gm[x][y][z]>outside) {
				// outside the GM : pull back
				forces[xyz] = -1.0f;
			} else if (subctx[x][y][z]>0) {
				// inside the subcortex : push out
				forces[xyz] = +1.0f;
			} else if (wmadjusted[xyz]<0) {
				// inside the previous segmentation: push out
				forces[xyz] = +1.0f;
			} else {
				//forces[xyz] = Numerics.bounded(2.0f*rfcm.getMemberships()[x][y][z][central]-1.0f, -1.0f, 1.0f);
				// or local max membership (nicer conceptually)
				forces[xyz] = Numerics.bounded(Numerics.max(rfcm.getMemberships()[x][y][z][inner],
															 rfcm.getMemberships()[x][y][z][central])
															-rfcm.getMemberships()[x][y][z][outer], -1.0f, 1.0f);
			}
		}
		
		Gdm3d gdm = new Gdm3d(levelset, 20.0f/rx, nx, ny, nz, rx, ry, rz, null,
								forces, 0.0f, gdmballoonParam.getValue().floatValue(), 
								gdmcurvParam.getValue().floatValue(), topologyParam.getValue());
			
		if (gdmiterParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("GDM myelinated GM adjustment...\n");
			
			gdm.evolveNarrowBand(gdmiterParam.getValue().intValue(), gdmdiffParam.getValue().floatValue());
		}
		float[] myelinated = gdm.exportLevelset();
		levelset = gdm.exportLevelset();
		
		if (adjustgmParam.getValue().booleanValue()) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (gm[x][y][z]>outside) {
					// outside the GM : pull back
					forces[xyz] = -1.0f;
				} else if (subctx[x][y][z]>0) {
					// inside the subcortex : push out
					forces[xyz] = +1.0f;
				} else if (myelinated[xyz]<0) {
					// inside the previous segmentation: push out
					forces[xyz] = +1.0f;
				} else {
					// balloon = bounded( 2 * myelinated - 1 , -1, 1)
					//forces[xyz] = Numerics.bounded(2.0f*rfcm.getMemberships()[x][y][z][outer]-1.0f, -1.0f, 1.0f);
					// or local max membership (nicer conceptually)
					forces[xyz] = Numerics.bounded(Numerics.max(rfcm.getMemberships()[x][y][z][inner],
																 rfcm.getMemberships()[x][y][z][central],
																 rfcm.getMemberships()[x][y][z][outer]), -1.0f, 1.0f);
				}
			}

			gdm = new Gdm3d(levelset, 20.0f/rx, nx, ny, nz, rx, ry, rz, null,
							 forces, 0.0f, gdmballoonParam.getValue().floatValue(), 
							 gdmcurvParam.getValue().floatValue(), topologyParam.getValue());
			
			if (gdmiterParam.getValue().intValue()>0) {
				BasicInfo.displayMessage("GDM final GM adjustment...\n");
				
				gdm.evolveNarrowBand(gdmiterParam.getValue().intValue(), gdmdiffParam.getValue().floatValue());
			}
			// copy over original GM surface
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				gm[x][y][z] = gdm.getLevelSet()[xyz];
			}
		}
		
		// 4. compute thickness and ratio
		float[][][] thickness = new float[nx][ny][nz];
		float[][][] ratio = new float[nx][ny][nz];
		float max = maxdistParam.getValue().floatValue()/rx;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			//if (mask[x][y][z]) {
			if (wmadjusted[xyz] > -max && gm[x][y][z] < max) {
				
				// thickness estimate by the minimum distance to corresp. boundaries
				thickness[x][y][z] = wmadjusted[xyz] - myelinated[xyz];
				ratio[x][y][z] = (wmadjusted[xyz] - myelinated[xyz])/(wmadjusted[xyz] - gm[x][y][z]);
				ratio[x][y][z] = Numerics.bounded(ratio[x][y][z],0,1);
				thickness[x][y][z] *= rx;
			}
		}
		
		// output
		byte[][][] classification = new byte[nx][ny][nz];
		float[][][] boundary = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (wmadjusted[xyz]<=0) classification[x][y][z] = 3;
			else if (myelinated[xyz]<=0) classification[x][y][z] = 2;
			else if (gm[x][y][z]<=0) classification[x][y][z] = 1;
			boundary[x][y][z] = myelinated[xyz];
		}
		
		ImageDataUByte cortexData = new ImageDataUByte(classification);		
		cortexData.setHeader(intensImage.getImageData().getHeader());
		cortexData.setName(intensImage.getImageData().getName()+"_myeseg");
		cortexImage.setValue(cortexData);
		cortexData = null;
		classification = null;
		
		ImageDataFloat boundaryData = new ImageDataFloat(boundary);		
		boundaryData.setHeader(intensImage.getImageData().getHeader());
		boundaryData.setName(intensImage.getImageData().getName()+"_myebound");
		boundaryImage.setValue(boundaryData);
		boundaryData = null;
		boundary = null;
		
		ImageDataFloat ratioData = new ImageDataFloat(ratio);		
		ratioData.setHeader(intensImage.getImageData().getHeader());
		ratioData.setName(intensImage.getImageData().getName()+"_myeratio");
		ratioImage.setValue(ratioData);
		ratioData = null;
		ratio = null;
		
		ImageDataFloat thickData = new ImageDataFloat(thickness);		
		thickData.setHeader(intensImage.getImageData().getHeader());
		thickData.setName(intensImage.getImageData().getName()+"_myethick");
		thickImage.setValue(thickData);
		thickData = null;
		thickness = null;
		
		if (debug) {
			ImageDataUByte rfcmData = new ImageDataUByte(rfcm.exportHardClassification());		
			rfcmData.setHeader(intensImage.getImageData().getHeader());
			rfcmData.setName(intensImage.getImageData().getName()+"_myerfcm");
			rfcmImage.setValue(rfcmData);
			rfcmData = null;
			
			ImageDataFloat memsData = new ImageDataFloat(rfcm.exportMemberships());		
			memsData.setHeader(intensImage.getImageData().getHeader());
			memsData.setName(intensImage.getImageData().getName()+"_myemems");
			memsImage.setValue(memsData);
			memsData = null;
			
			wm = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				wm[x][y][z] = wmadjusted[xyz];
			}
			ImageDataFloat innerData = new ImageDataFloat(wm);		
			innerData.setHeader(intensImage.getImageData().getHeader());
			innerData.setName(intensImage.getImageData().getName()+"_myeinner");
			innerImage.setValue(innerData);
			innerData = null;
			wm = null;
		
			ImageDataFloat outerData = new ImageDataFloat(gm);		
			outerData.setHeader(intensImage.getImageData().getHeader());
			outerData.setName(intensImage.getImageData().getName()+"_myeouter");
			outerImage.setValue(outerData);
			outerData = null;
			gm = null;
		}
	}

}
