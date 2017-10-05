package de.mpg.cbs.jist.registration;

import edu.jhu.ece.iacl.jist.io.ImageDataReaderWriter;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.*;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
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

import org.apache.commons.math3.util.FastMath;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/*
 * @authors Christine Lucas Tardif and Pierre-Louis Bazin
 */

/*
 * source curvature is calculated before the levelset is deformed 
 */

/*
 * TO DO:
 * change surface inflation method (weighted by angle with surface normal)
 * inflation scales determined by global curvature, not by number of iterations?
 * try to speedup
 */

public class JistRegistrationMultimodalSurface extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	sourceLevelsetImage;
	private ParamVolume 	targetLevelsetImage;
	private ParamFloat 	wcurvParam;
	
	private ParamVolume 	sourceIntensityImage;
	private ParamVolume 	targetIntensityImage;
	private ParamFloat 	wintensParam;
	
	private	 ParamFile		synScript1Param;
	private	 ParamFile		synScript2Param;
	private	 ParamFile		synScript3Param;
	private	 ParamFile		synScript4Param;
	
	//private ParamVolume deformedImage;
	private ParamVolume mappingImage;
	private ParamVolume invmappingImage;
	
	private ParamBoolean atlasParam;
	private ParamBoolean phase2Param;
	
	private ParamBoolean	debugParam;
	
	// global variables
	private int nx, ny, nz, nxyz, nimg;
	private float rx, ry, rz;
	
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(sourceLevelsetImage = new ParamVolume("Source Levelset Surface"));
		inputParams.add(targetLevelsetImage = new ParamVolume("Target Levelset Surface"));
		inputParams.add(sourceIntensityImage = new ParamVolume("Source Contrasts (optional)"));
		sourceIntensityImage.setMandatory(false);
		inputParams.add(targetIntensityImage = new ParamVolume("Target Contrasts (optional)"));
		targetIntensityImage.setMandatory(false);
		
		inputParams.add(synScript1Param = new ParamFile("external SyN script 1 contrasts"));
		inputParams.add(synScript2Param = new ParamFile("external SyN script 2 contrasts"));
		inputParams.add(synScript3Param = new ParamFile("external SyN script 3 contrasts"));
		inputParams.add(synScript4Param = new ParamFile("external SyN script 4 contrasts"));
		inputParams.add(wintensParam = new ParamFloat("Intensity weighting (optional)", 0.0f, 2.0f, 1.0f));
		inputParams.add(wcurvParam = new ParamFloat("Curvature weighting", 0.0f, 2.0f, 1.0f));
		inputParams.add(atlasParam = new ParamBoolean("Is the target an atlas?", false));
		inputParams.add(phase2Param = new ParamBoolean("Align original curvature?", true));
		inputParams.add(debugParam = new ParamBoolean("Keep intermediate files", true));
				
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration");
		inputParams.setLabel("Multimodal Surface Registration");
		inputParams.setName("MultimodalSurfaceRegistration");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Registers a source cortical surface to a target using geometry (including the levelset, (optional) curvedness and shape index) and (optional) additional contrasts (eg T1). The multi-scale approach works with partially inflated levelset surfaces.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		//outputParams.add(deformedImage = new ParamVolume("Mapped source levelset",VoxelType.FLOAT));
		outputParams.add(mappingImage = new ParamVolume("Mapping function",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(invmappingImage = new ParamVolume("Inverse Mapping function",VoxelType.FLOAT,-1,-1,-1,-1));
			
		outputParams.setName("surface registration");
		outputParams.setLabel("surface registration");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data
		//-----------------------
		BasicInfo.displayMessage("Importing data.\n");
		ImageDataFloat	slevelsetImg = new ImageDataFloat(sourceLevelsetImage.getImageData());
		ImageDataFloat	tlevelsetImg = new ImageDataFloat(targetLevelsetImage.getImageData());
		
		nx = slevelsetImg.getRows();
		ny = slevelsetImg.getCols();
		nz = slevelsetImg.getSlices();
		nxyz = nx*ny*nz;
		rx = slevelsetImg.getHeader().getDimResolutions()[0];
		ry = slevelsetImg.getHeader().getDimResolutions()[1];
		rz = slevelsetImg.getHeader().getDimResolutions()[2];
		
		float[] slevelset = new float[nxyz];
		float[][][] buffer = slevelsetImg.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			slevelset[xyz] = buffer[x][y][z];
		}
		buffer = null;
		slevelsetImg = null;
		
		float[] tlevelset = new float[nxyz];
		buffer = tlevelsetImg.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			tlevelset[xyz] = buffer[x][y][z];
		}
		buffer = null;
		tlevelsetImg = null;
		
		float[][] scurv = null;
		float[][] tcurv = null;
		
		nimg=0;
		float[][] sintensity = null;
		float[][] tintensity = null;
		if (sourceIntensityImage.getImageData()!=null && targetIntensityImage.getImageData()!=null) {
			ImageDataFloat sintensImg = new ImageDataFloat(sourceIntensityImage.getImageData());
			ImageDataFloat tintensImg = new ImageDataFloat(targetIntensityImage.getImageData());
			nimg = Numerics.min(sintensImg.getComponents(),tintensImg.getComponents());
			
			sintensity = new float[nimg][nxyz];
			tintensity = new float[nimg][nxyz];
			if (nimg>1) {
				// 4D data set
				float[][][][] buffer4d = sintensImg.toArray4d();
				for (int d=0;d<nimg;d++) for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					sintensity[d][xyz] = buffer4d[x][y][z][d];
				}
				buffer4d = tintensImg.toArray4d();
				for (int d=0;d<nimg;d++) for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					tintensity[d][xyz] = buffer4d[x][y][z][d];
				}
				buffer4d = null;
			} else {
				// 3D data set
				buffer = sintensImg.toArray3d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					sintensity[0][xyz] = buffer[x][y][z];
				}
				buffer = tintensImg.toArray3d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					tintensity[0][xyz] = buffer[x][y][z];
				}
				buffer = null;
			}
		} 
		
		// Initializing deformations
		//----------------------------
		BasicInfo.displayMessage("Initializing deformations\n");
		
		// deform = Id; invdeform = Id;
		float[][] mapping = new float[3][nxyz];
		float[][] invmapping = new float[3][nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			mapping[X][xyz] = x;		
			mapping[Y][xyz] = y;		
			mapping[Z][xyz] = z;		
			invmapping[X][xyz] = x;		
			invmapping[Y][xyz] = y;		
			invmapping[Z][xyz] = z;		
		}
		
		// init coordinates for SyN processing
		float[][][][] coord = new float[nx][ny][nz][3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			coord[x][y][z][X] = x;		
			coord[x][y][z][Y] = y;		
			coord[x][y][z][Z] = z;		
		}
		
		// reinit the source and target levelset (only needed once)
		//----------------------------------------------------------
		BasicInfo.displayMessage("Levelset reinitialization\n");
		
		boolean[] bgmask = new boolean[nxyz];
		boolean[] bgmask2 = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) bgmask[xyz] = true; 
		InflateGdm reinit = new InflateGdm(slevelset, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no",null);
		reinit.fastMarchingReinitialization(false);
		slevelset = reinit.getLevelSet();

		reinit = new InflateGdm(tlevelset, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no",null);
		reinit.fastMarchingReinitialization(false);
		tlevelset = reinit.getLevelSet();		
		
		// pointers to the original data (to be kept unchanged)
		float[] initslevelset = slevelset.clone();
		float[] inittlevelset = tlevelset.clone();
		
		float[][] initsintensity = null; 
		float[][] inittintensity = null;
		
		if (sintensity != null) {
			initsintensity = sintensity.clone();
			inittintensity = tintensity.clone();
		}
		
		// Multiscale parameters
		//-----------------------
		int nbscales = 20;
		
		int[] smoothiter = new int[nbscales];
		/*
		smoothiter[0] = 45;
		smoothiter[1] = 40; 
		smoothiter[2] = 37; 
		smoothiter[3] = 34; 
		smoothiter[4] = 31; 
		smoothiter[5] = 28; 
		smoothiter[6] = 26; 
		smoothiter[7] = 24; 
		smoothiter[8] = 22; 
		smoothiter[9] = 20;
		smoothiter[10] = 18;
		smoothiter[11] = 16;
		smoothiter[12] = 14;
		smoothiter[13] = 12;
		smoothiter[14] = 10;
		smoothiter[15] = 9; 
		smoothiter[16] = 8; 
		smoothiter[17] = 7;
		smoothiter[18] = 6;
		smoothiter[19] = 5;
		smoothiter[20] = 4;
		smoothiter[21] = 3;
		smoothiter[22] = 2;
		smoothiter[23] = 1;
		smoothiter[24] = 0;
		*/
		smoothiter[0] = 40; 
		smoothiter[1] = 36; 
		smoothiter[2] = 32; 
		smoothiter[3] = 28; 
		smoothiter[4] = 25; 
		smoothiter[5] = 22;
		smoothiter[6] = 19;
		smoothiter[7] = 17;
		smoothiter[8] = 15;
		smoothiter[9] = 13;
		smoothiter[10] = 11;
		smoothiter[11] = 9; 
		smoothiter[12] = 7;
		smoothiter[13] = 6;
		smoothiter[14] = 5;
		smoothiter[15] = 4;
		smoothiter[16] = 3;
		smoothiter[17] = 2;
		smoothiter[18] = 1;
		smoothiter[19] = 0;
		
		boolean affine = true;
		int coarseiter = 100; 
		int mediter = 100;
		int fineiter = 0;
		int cutoff = 14;
		
		float smoothing = 1.5f/rx;
		
		double globalcurv = 0.0f;
		float maxdist = 0.0f;
		float maxdist2 = 0.0f;
		
		boolean[] smask = new boolean[nxyz];
		float[][][] slvlProb=null, tlvlProb=null;
		float[][][][] scurvProb=null, tcurvProb=null;
		float[] smoothslevelset = null;
		float[][][][] sintensProb=null, tintensProb=null;
	
		float kscale = 2.0f;
		float[][] kernel = ImageFilters.separableGaussianKernel(kscale, kscale, kscale);
		int ks = (kernel[0].length-1)/2;
		
		
		for (int sc=0; sc<nbscales; sc++) {
											
			BasicInfo.displayMessage("Registration scale "+(sc+1)+" of "+nbscales+"\n");
			
			slevelset = initslevelset.clone();
			tlevelset = inittlevelset.clone();
			
			if (initsintensity != null) {
				sintensity = initsintensity.clone();
				tintensity = inittintensity.clone();
			}
			
			// 1. inflate source and target levelsets
			//----------------------------------------
			BasicInfo.displayMessage("Surface inflation ("+smoothiter[sc]+" iter)\n");
			
			if(sintensity != null) {
				if (sc < nbscales-1) {
					maxdist = levelsetInflationMapping(slevelset, sintensity, smoothing, smoothiter[sc], smoothiter[sc+1], 0.0f);
				} else {
					maxdist = levelsetInflationMapping(slevelset, sintensity, smoothing, smoothiter[sc], smoothiter[sc], 0.0f);
				}
			} else {
				if (sc < nbscales-1) {
					maxdist = levelsetInflation(slevelset, smoothing, smoothiter[sc], smoothiter[sc+1], 0.0f);
				} else {
					maxdist = levelsetInflation(slevelset, smoothing, smoothiter[sc], smoothiter[sc], 0.0f);
				}
			}
			
			BasicInfo.displayMessage("max. distance = "+maxdist+" voxels\n"); 

			for (int xyz=0;xyz<nxyz;xyz++) {
				if(FastMath.abs(slevelset[xyz])<2.0f*1.7321f) {
					smask[xyz] = true;
				} else {
					smask[xyz] = false;
				}
			}
			
			// smooth the levelset to calculate global curvature
			smoothslevelset = ImageFilters.separableConvolution(slevelset, nx, ny, nz, kernel, ks, ks, ks);
			globalcurv = ImageStatistics.mean(ImageGeometry.squaredMeanCurvature(smoothslevelset, nx, ny, nz, smask), smask, nx, ny, nz);
			smoothslevelset = null;
			
			BasicInfo.displayMessage("nominalglobalcurv = "+globalcurv+"\n");
			
			if(sintensity != null) {
				if (atlasParam.getValue().booleanValue()) {
					levelsetInflationMapping(tlevelset, tintensity, smoothing, smoothiter[sc], smoothiter[sc], globalcurv);
				} else {
					levelsetInflationMapping(tlevelset, tintensity, smoothing, smoothiter[sc], smoothiter[sc], 0.0f);
				}
			} else {
				if (atlasParam.getValue().booleanValue()) {
					levelsetInflation(tlevelset, smoothing, smoothiter[sc], smoothiter[sc], globalcurv);
				} else {
					levelsetInflation(tlevelset, smoothing, smoothiter[sc], smoothiter[sc], 0.0f);
				}
			}
			
			// 2. calculate curvature metrics
			//--------------------------------
			if (wcurvParam.getValue().floatValue()>0.0f) {
				
				if (!phase2Param.getValue().booleanValue() && smoothiter[sc]<cutoff) {

					scurv = null;
					tcurv = null;

				} else {

					scurv = new float[2][nxyz];
					tcurv = new float[2][nxyz];

					// smooth the levelsets
					float[] sslevelset = ImageFilters.separableConvolution(slevelset, nx, ny, nz, kernel, ks, ks, ks);
					float[] stlevelset = ImageFilters.separableConvolution(tlevelset, nx, ny, nz, kernel, ks, ks, ks);

					// calculate the curvature metrics
					for (int x=1; x<nx-1; x++) for (int y=1; y<ny-1; y++) for (int z = 1; z<nz-1; z++) {
						int xyz = x + nx*y + nx*ny*z;

						if(FastMath.abs(slevelset[xyz]) < 2.0f*1.7321f) {
							scurv[0][xyz] = ImageGeometry.curvednessValue(sslevelset, xyz, nx, ny, nz);
							scurv[1][xyz] = ImageGeometry.shapeIndexValue(sslevelset, xyz, nx, ny, nz);	
						} else {
							scurv[0][xyz] = 0.0f;
							scurv[1][xyz] = 0.0f;
						}

						if(FastMath.abs(tlevelset[xyz]) < 2.0f*1.7321f) {
							tcurv[0][xyz] = ImageGeometry.curvednessValue(stlevelset, xyz, nx, ny, nz);
							tcurv[1][xyz] = ImageGeometry.shapeIndexValue(stlevelset, xyz, nx, ny, nz);	
						} else {
							tcurv[0][xyz] = 0.0f;
							tcurv[1][xyz] = 0.0f;
						}

					}
					sslevelset = null;
					stlevelset = null;
				}
			}

			// 3. deform the source data
			//----------------------------
			if (sc>0) {
				BasicInfo.displayMessage("combine with current deformation\n");
				
				float[] dlevelset = new float[nxyz];
				float[][] dcurv = null;
				if (scurv != null) dcurv = new float[2][nxyz];
				float[][] dintensity = null;
				if (sintensity != null) dintensity = new float[nimg][nxyz];
				boolean[] dbgmask = new boolean[nxyz];
				
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (mapping[X][xyz]>=1 && mapping[Y][xyz]>=1 && mapping[Z][xyz]>=1) {
						dlevelset[xyz] = ImageInterpolation.linearClosestInterpolation(slevelset, mapping[X][xyz], mapping[Y][xyz], mapping[Z][xyz], nx, ny, nz);
						dbgmask[xyz] = ImageInterpolation.nearestNeighborInterpolation(bgmask, false, mapping[X][xyz], mapping[Y][xyz], mapping[Z][xyz], nx, ny, nz);
						
						if (scurv != null) {
							dcurv[0][xyz] = ImageInterpolation.linearClosestInterpolation(scurv[0], mapping[X][xyz], mapping[Y][xyz], mapping[Z][xyz], nx, ny, nz);
							dcurv[1][xyz] = ImageInterpolation.linearClosestInterpolation(scurv[1], mapping[X][xyz], mapping[Y][xyz], mapping[Z][xyz], nx, ny, nz);
						} 

						if (sintensity != null) {
							for (int n=0;n<nimg;n++)
								dintensity[n][xyz] = ImageInterpolation.linearClosestInterpolation(sintensity[n], mapping[X][xyz], mapping[Y][xyz], mapping[Z][xyz], nx, ny, nz);
						}
					}
				}
				
				slevelset = dlevelset;
				if (scurv != null) scurv = dcurv;
				if (sintensity != null) sintensity = dintensity;
								
				// reinit the source levelset
				dbgmask = Morphology.erodeObject(dbgmask, nx, ny, nz, 1, 1, 1);
				reinit = new InflateGdm(slevelset, nx, ny, nz, rx, ry, rz, dbgmask, 0.4f, 0.4f, "no",null);
				reinit.fastMarchingReinitialization(false);
				slevelset = reinit.getLevelSet();
			}
			
			// 4. Find the scale of distance between levelsets
			//--------------------------------------------------
			BasicInfo.displayMessage("compute distance parameter:\n");
			maxdist2 = 0.0f;
			float dist;
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
			
				if (x>=1 || x<=nx-2 || y>=1 || y<=ny-2 || z>=1 || z<=nz-2) {
					//if (Numerics.abs(slevelset[xyz]) < 1.7321f || Numerics.abs(tlevelset[xyz]) < 1.7321f) {
					if (FastMath.abs(tlevelset[xyz]) < 1.7321f) {
						dist = Numerics.abs(slevelset[xyz]-tlevelset[xyz]);
						if (dist > maxdist2) maxdist2 = dist;
					}
				}
			}
			
			BasicInfo.displayMessage("max. distance 2 = "+maxdist2+" voxels\n"); 
			maxdist = Numerics.max(maxdist, maxdist2);		
			BasicInfo.displayMessage("final max. distance = "+maxdist+" voxels\n"); 

			// 5. create leveset-based contrasts
			//------------------------------------
			BasicInfo.displayMessage("normalized registration contrasts\n");
			
			slvlProb = new float[nx][ny][nz];
			if(scurv != null) scurvProb = new float[2][nx][ny][nz];
			if(sintensity != null) 	sintensProb = new float[nimg][nx][ny][nz];
			generateRegistrationContrasts(slevelset, scurv, sintensity, slvlProb, scurvProb, sintensProb, maxdist);
						
			tlvlProb = new float[nx][ny][nz];
			if(tcurv != null) tcurvProb = new float[2][nx][ny][nz];
			if(tintensity != null) 	tintensProb = new float[nimg][nx][ny][nz];
			generateRegistrationContrasts(tlevelset, tcurv, tintensity, tlvlProb, tcurvProb, tintensProb, maxdist);	
			
			
			// 6. run SyN + update deformations
			//-----------------------------------
			runSyNscript(sc, slvlProb, scurvProb, sintensProb, tlvlProb, tcurvProb, tintensProb, coord,
					mapping, invmapping, coarseiter, mediter, fineiter, affine);	
			
			affine = false;
		}
		
		
		// Output
		//----------
		BasicInfo.displayMessage("Output...\n");
			
		float[][][][] map = new float[nx][ny][nz][3];
		float[][][][] invmap = new float[nx][ny][nz][3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int n=0;n<3;n++) {
			int xyz = x+nx*y+nx*ny*z;
			map[x][y][z][n] = mapping[n][xyz];
			invmap[x][y][z][n] = invmapping[n][xyz];
		}
		
		String imgname = sourceLevelsetImage.getImageData().getName();
		
		ImageDataFloat mappingData = new ImageDataFloat(map);
		mappingData.setHeader(targetLevelsetImage.getImageData().getHeader());
		mappingData.setName(imgname+"_map");
		mappingImage.setValue(mappingData);
		mappingData = null;
		map = null;
		
		ImageDataFloat invmappingData = new ImageDataFloat(invmap);
		invmappingData.setHeader(sourceLevelsetImage.getImageData().getHeader());
		invmappingData.setName(imgname+"_invmap");
		invmappingImage.setValue(invmappingData);
		invmappingData = null;
		invmap = null;
	}

	private final float levelsetInflation (float[] levelset, float smoothing, int smoothiter, int smoothiter2, double nominalglobalcurv) {
		
		int itr = smoothiter;
		double realglobalcurv = 0.0f;
		float maxdist = 0.0f;
		boolean[] surfmask = new boolean[nxyz];
		
		// create a mask of a narrowband
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = (FastMath.abs(levelset[xyz]) <= 5.0f);
		}
		
		// important: increase the data range (useful for better smoothing)
		int delta = 30;
		ObjectMorphology.fastDilateObject(bgmask, nx, ny, nz, delta);

		// main algorithm
		BasicInfo.displayMessage("re-build levelset...\n");
		InflateGdm gdm = new InflateGdm(levelset, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no",null);
		
		double basis = 1.0f;
		double scale = smoothing;
		if (smoothiter>1) {
			basis = Math.pow(2.0f*smoothing, 1.0f/(smoothiter-1.0f));
			scale = 0.5f;
		}
		
		float[] newlevelset = new float[nxyz];
		float[] smoothlevelset = new float[nxyz];
		float[] nextlevelset = new float[nxyz];
		
		float kscale = 2.0f;
		float[][] kernel = ImageFilters.separableGaussianKernel(kscale, kscale, kscale);
		int ks = (kernel[0].length-1)/2;
		
		for (int t=0;t<smoothiter;t++) {
			BasicInfo.displayMessage(".");
			
			gdm.smoothLevelset((float) scale);
			gdm.evolveNarrowBand(500, 0.001f);
			gdm.updateTarget();
			scale *= basis;
			
			newlevelset = gdm.getLevelSet();
			
			if (t == smoothiter2) {
				for (int xyz=0;xyz<nxyz;xyz++) {
					nextlevelset[xyz] = newlevelset[xyz];
				}
				BasicInfo.displayMessage("t = "+t+"\n");
			}
			
			// check if the global shape metric is the same
			if (nominalglobalcurv > 0.0f) {
				
				for (int xyz=0;xyz<nxyz;xyz++) {
					if(FastMath.abs(newlevelset[xyz]) < 1.7321f) {
						surfmask[xyz] = true;
					} else {
						surfmask[xyz] = false;
					}
				}
				
				// smooth the surface
				smoothlevelset = ImageFilters.separableConvolution(newlevelset, nx, ny, nz, kernel, ks, ks, ks);
								
				realglobalcurv = ImageStatistics.mean(ImageGeometry.squaredMeanCurvature(smoothlevelset, nx, ny, nz, surfmask), surfmask, nx, ny, nz);
				BasicInfo.displayMessage("realglobalcurv = "+realglobalcurv+"\n");
				
				if (realglobalcurv<nominalglobalcurv) {
					itr = t+1;
					t = smoothiter;
				} 
			} 
			
		}
		
		newlevelset = gdm.getLevelSet();
		
		// calculate max distance between current and next smoothing scale levelsets
		if (smoothiter2 < smoothiter) {
			BasicInfo.displayMessage("calculate max distance \n");
			float dist;
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x + nx*y + nx*ny*z;
			
				if (FastMath.abs(nextlevelset[xyz]) < 2.0f*1.7321f) {
					dist = FastMath.abs(newlevelset[xyz]-nextlevelset[xyz]);
					if (dist>maxdist) maxdist = dist;
				}
				
			}	
		}
		
		// replace original values
		for (int xyz=0;xyz<nxyz;xyz++) {
			levelset[xyz] = newlevelset[xyz];
		}
		
		return maxdist;
	}
	
	private final float levelsetInflationMapping(float[] levelset, float[][] intensity, float smoothing, int smoothiter, int smoothiter2, double nominalglobalcurv) {
		
		int itr = smoothiter;
		double realglobalcurv = 0.0f;
		float maxdist = 0.0f;
		boolean[] surfmask = new boolean[nxyz];
		
		// create a mask of a narrowband
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = (FastMath.abs(levelset[xyz])<=5.0f);
		}
		// important: increase the data range (useful for better smoothing)
		int delta = 30;
		ObjectMorphology.fastDilateObject(bgmask, nx, ny, nz, delta);
		
		// create a mappping mask
		boolean[] intensitymask = null;
		if (intensity != null) {
			intensitymask = new boolean[nxyz];
			for (int xyz=0;xyz<nxyz;xyz++) {
				intensitymask[xyz] = (intensity[0][xyz]>0.0f);
			}
		}
		
		BasicInfo.displayMessage("re-build levelset...\n");

		// main algorithm
		CorticalInflationGdm gdm = new CorticalInflationGdm(levelset, nx, ny, nz, rx, ry, rz, 
																bgmask, 0.4f, 0.4f, "no",
																intensity, intensitymask, nimg);
		
		double basis = 1.0f;
		double scale = smoothing;
		if (smoothiter>1) {
			basis = Math.pow(2.0f*smoothing, 1.0f/(smoothiter-1.0f));
			scale = 0.5f;
		}
		
		float[] newlevelset = new float[nxyz];
		float[] smoothlevelset = new float[nxyz];
		float[] nextlevelset = new float[nxyz];
		
		float kscale = 2.0f;
		float[][] kernel = ImageFilters.separableGaussianKernel(kscale, kscale, kscale);
		int ks = (kernel[0].length-1)/2;
		
		for (int t=0;t<smoothiter;t++) {
			BasicInfo.displayMessage(".");
			
			gdm.computeSmoothedLevelset((float)scale, false);
			//gdm.evolveNarrowBandMapping(500, 0.001f);
			gdm.evolveNarrowBandMappingMean(500, 0.001f);
			
			scale *= basis;
			
			newlevelset = gdm.getLevelSet();
			
			if (t == smoothiter2) {
				for (int xyz=0;xyz<nxyz;xyz++) {
					nextlevelset[xyz] = newlevelset[xyz];
				}
				BasicInfo.displayMessage("t = "+t+"\n");
			}
			
			// check if the global shape metric is the same
			if (nominalglobalcurv > 0.0f) {
				
				for (int xyz=0;xyz<nxyz;xyz++) {
					if(FastMath.abs(newlevelset[xyz]) < 1.7321f) {
						surfmask[xyz] = true;
					} else {
						surfmask[xyz] = false;
					}
				}
				
				// smooth the surface
				smoothlevelset = ImageFilters.separableConvolution(newlevelset, nx, ny, nz, kernel, ks, ks, ks);
								
				realglobalcurv = ImageStatistics.mean(ImageGeometry.squaredMeanCurvature(smoothlevelset, nx, ny, nz, surfmask), surfmask, nx, ny, nz);
				BasicInfo.displayMessage("realglobalcurv = "+realglobalcurv+"\n");
				
				if (realglobalcurv<nominalglobalcurv) {
					itr = t+1;
					t = smoothiter;
				} 
			} 
			
		}
		
		newlevelset = gdm.getLevelSet();
		
		gdm.cleanupForwardMapping();
		gdm.cleanupBackwardMapping();
		
		// calculate max distance between current and next smoothing scale levelsets
		if (smoothiter2 < smoothiter) {
			BasicInfo.displayMessage("calculate max distance \n");
			float dist;
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x + nx*y + nx*ny*z;
			
				if (FastMath.abs(nextlevelset[xyz]) < 2.0f*1.7321f) {
					dist = FastMath.abs(newlevelset[xyz]-nextlevelset[xyz]);
					if (dist>maxdist) maxdist = dist;
				}
			
			}	
		}
		
		// replace original values
		if (intensity != null) {
			float[][] newintensity = gdm.getMappedIntensity();
			for (int xyz=0;xyz<nxyz;xyz++) {
				levelset[xyz] = newlevelset[xyz];
				for (int n=0;n<nimg;n++)
					intensity[n][xyz] = newintensity[n][xyz];
			}
		} else {
			for (int xyz=0;xyz<nxyz;xyz++) {
				levelset[xyz] = newlevelset[xyz];
			}
		}
		
		return maxdist;
	}
	
	private final void generateRegistrationContrasts(float[] levelset, float[][] curv, float[][] intensity, 
														float[][][] lvlProb, float[][][][] curvProb, float[][][][] intensProb, 
														float distance) {
		
		if (distance == 0.0f) distance = 2.0f; 
				
		// create a mask
		boolean[] mask = new boolean[nxyz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) {
				mask[xyz] = false;
			} else {
				mask[xyz] = (FastMath.abs(levelset[xyz]) <= 1.0f*distance/rx);
			}
		}
		
		// calculate levelset probability
		System.out.println("\n Calculating levelset probability");
		float scale = 0.5f*distance;
		
		for (int x=0;x<nx;x++) {
			System.out.print(".");
			for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				lvlProb[x][y][z] = (float)(1.0f/(1.0f+FastMath.exp(levelset[xyz]*rx/scale)));
			}
		}
		
		float m, b;
		CorticalSurfaceDilation surfdil = null;
		
		// calculate curvature probability
		if (curv != null) {
			System.out.println("\n Calculating curvature probability");
						
			//dilate curv values from levelset
			surfdil = new CorticalSurfaceDilation(curv[0], null, levelset, mask, nx, ny, nz, rx, ry, rz, "26");
			surfdil.dilateData();
			curv[0] = surfdil.getDilatedData();
			
			surfdil = new CorticalSurfaceDilation(curv[1], null, levelset, mask, nx, ny, nz, rx, ry, rz, "26");
			surfdil.dilateData();
			curv[1] = surfdil.getDilatedData();
			
			float curvMin = ImageStatistics.robustMinimum(curv[0], 0.005f, 4, nx, ny, nz);
			float curvMax = ImageStatistics.robustMaximum(curv[0], 0.005f, 4, nx, ny, nz);
			
			System.out.println("\n robust min curvedness = "+curvMin);
			System.out.println("\n robust max curvedness = "+curvMax);
			
			m = 1.0f/(curvMax-curvMin);
			b = 0.5f-m*curvMax;
	
			for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
				int xyz = x + nx*y + nx*ny*z;
				if(mask[xyz]) {
					curvProb[0][x][y][z] = m*curv[0][xyz]+b; 
					curvProb[0][x][y][z] = Numerics.bounded(curvProb[0][x][y][z], -0.5f, 0.5f);
				} else {
					curvProb[0][x][y][z] = 0.0f;
				}
			}
			
			curvMin = ImageStatistics.robustMinimum(curv[1], 0.001f, 4, nx, ny, nz);
			curvMax = ImageStatistics.robustMaximum(curv[1], 0.001f, 4, nx, ny, nz);
			
			System.out.println("\n robust min shape index = "+curvMin);
			System.out.println("\n robust max shape index = "+curvMax);
			
			m = 1.0f/(curvMax-curvMin);
			b = 0.5f-m*curvMax;
	
			for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
				int xyz = x + nx*y + nx*ny*z;
				if(mask[xyz]) {
					curvProb[1][x][y][z] = m*curv[1][xyz]+b; 
					curvProb[1][x][y][z] = Numerics.bounded(curvProb[1][x][y][z], -0.5f, 0.5f);
				} else {
					curvProb[1][x][y][z] = 0.0f;
				}
			}
		}
		
		// calculating intensity probability image
		if (intensity != null) {
			System.out.println("\n Calculating T1 probability");
			
			for (int n=0; n<nimg; n++) {

				float[] intens = intensity[n];
				
				//dilate values from levelset
				surfdil = new CorticalSurfaceDilation(intens, null, levelset, mask, nx, ny, nz, rx, ry, rz, "26");
				surfdil.dilateData();
				intens = surfdil.getDilatedData();	
				
				float intensMin = ImageStatistics.robustMinimum(intens, 0.001f, 4, nx, ny, nz);
				float intensMax = ImageStatistics.robustMaximum(intens, 0.001f, 4, nx, ny, nz);
				
				System.out.println("\n robust min intensity = "+intensMin);
				System.out.println("\n robust max intensity = "+intensMax);	
		
				m = 1.0f/(intensMax-intensMin);
				b = 0.5f-m*intensMax;
				for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
					int xyz = x + nx*y + nx*ny*z;
					if(mask[xyz]) {
						intensProb[n][x][y][z] = m*intens[xyz]+b;
						intensProb[n][x][y][z] = Numerics.bounded(intensProb[n][x][y][z], -0.5f, 0.5f);
					} else {
						intensProb[n][x][y][z] = 0.0f;
					}
				}
			}
		}
		
		return;
	}
	
	private final void runSyNscript(int sc, float[][][] slvl, float[][][][] scurv, float[][][][] sintens, 
			 float[][][] tlvl, float[][][][] tcurv, float[][][][] tintens,
			 float[][][][] coord,
			 float[][] map, float[][] invmap, 
			 int coarseIt, int medIt, int fineIt, boolean affine) {

		try {
			ImageDataReaderWriter readerWriter = new ImageDataReaderWriter();
			File destdir = new File(this.getOutputDirectory().getCanonicalFile()+File.separator+this.getAlgorithmName());

			// create the directory if it doesn't exist
			if(!destdir.isDirectory()){
				(new File(destdir.getCanonicalPath())).mkdir();
			}
			// create basic headers
			ImageHeader sHeader = new ImageHeader();
			sHeader.setDimResolutions(sourceLevelsetImage.getImageData().getHeader().getDimResolutions());

			ImageHeader tHeader = new ImageHeader();
			tHeader.setDimResolutions(targetLevelsetImage.getImageData().getHeader().getDimResolutions());

			// write the temporary files to be processed
			String inputSource1Loc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc) + "_slvl.nii";
			File inputSource1 = new File(inputSource1Loc);
			ImageDataFloat slvlImg = new ImageDataFloat(slvl);
			slvlImg.setHeader(sHeader);
			readerWriter.write(slvlImg, inputSource1);

			String inputTarget1Loc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc) + "_tlvl.nii";
			File inputTarget1 = new File(inputTarget1Loc);
			ImageDataFloat tlvlImg = new ImageDataFloat(tlvl);
			tlvlImg.setHeader(tHeader);
			readerWriter.write(tlvlImg, inputTarget1);

			String inputCoordLoc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc)  + "_coord.nii";
			File inputCoord = new File(inputCoordLoc);
			ImageDataFloat coordImg = new ImageDataFloat(coord);
			coordImg.setHeader(sHeader);
			readerWriter.write(coordImg, inputCoord);

			String outputDeform1Loc = destdir.getAbsolutePath() + File.separator + "output_scale" + Integer.toString(sc) + "_dlvl.nii";
			File outputDeform1 = new File(outputDeform1Loc);

			File inputSource2 = null;
			File inputTarget2 = null;
			File outputDeform2 = null;
			File inputSource3 = null;
			File inputTarget3 = null;
			File outputDeform3 = null;

			if (scurv != null) {
				String inputSource2Loc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc) + "_scurvedness.nii";
				inputSource2 = new File(inputSource2Loc);
				ImageDataFloat scurvImg = new ImageDataFloat(scurv[0]);
				scurvImg.setHeader(sHeader);
				readerWriter.write(scurvImg, inputSource2);			

				String inputSource3Loc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc) + "_sshapeindex.nii";
				inputSource3 = new File(inputSource3Loc);
				scurvImg = new ImageDataFloat(scurv[1]);
				scurvImg.setHeader(sHeader);
				readerWriter.write(scurvImg, inputSource3);			

				String inputTarget2Loc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc) + "_tcurvedness.nii";
				inputTarget2 = new File(inputTarget2Loc);
				ImageDataFloat tcurvImg = new ImageDataFloat(tcurv[0]);
				tcurvImg.setHeader(tHeader);
				readerWriter.write(tcurvImg, inputTarget2);

				String inputTarget3Loc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc) + "_tshapeindex.nii";
				inputTarget3 = new File(inputTarget3Loc);
				tcurvImg = new ImageDataFloat(tcurv[1]);
				tcurvImg.setHeader(tHeader);
				readerWriter.write(tcurvImg, inputTarget3);

				String outputDeform2Loc = destdir.getAbsolutePath() + File.separator + "output_scale" + Integer.toString(sc) + "_dcurvedness.nii";
				outputDeform2 = new File(outputDeform2Loc);

				String outputDeform3Loc = destdir.getAbsolutePath() + File.separator + "output_scale" + Integer.toString(sc) + "_dshapeindex.nii";
				outputDeform3 = new File(outputDeform3Loc);
			}


			File[] inputSource4 = new File[nimg];
			File[] inputTarget4 = new File[nimg];
			File[] outputDeform4 = new File[nimg];

			for (int n=0; n<nimg; n++) {
				String inputSource4Loc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc) + "_sintens"+n+".nii";
				inputSource4[n] = new File(inputSource4Loc);
				ImageDataFloat sintensImg = new ImageDataFloat(sintens[n]);
				sintensImg.setHeader(sHeader);
				readerWriter.write(sintensImg, inputSource4[n]);

				String inputTarget4Loc = destdir.getAbsolutePath() + File.separator + "input_scale" + Integer.toString(sc) + "_tintens"+n+".nii";
				inputTarget4[n] = new File(inputTarget4Loc);
				ImageDataFloat tintensImg = new ImageDataFloat(tintens[n]);
				tintensImg.setHeader(tHeader);
				readerWriter.write(tintensImg, inputTarget4[n]);

				String outputDeform4Loc = destdir.getAbsolutePath() + File.separator + "output_scale" + Integer.toString(sc) + "_dintens"+n+".nii";
				outputDeform4[n] = new File(outputDeform4Loc);
			}


			String outputMappingLoc = destdir.getAbsolutePath() + File.separator + "output_scale" + Integer.toString(sc) + "_map.nii";
			File outputMapping = new File(outputMappingLoc);

			String outputInvMappingLoc = destdir.getAbsolutePath() + File.separator + "output_scale" + Integer.toString(sc) + "_invmap.nii";
			File outputInvMapping = new File(outputInvMappingLoc);


			//build command array
			String[] comArray = new String[40];

			int totalContrasts = 1+nimg;
			if (scurv != null) totalContrasts += 2;

			if (totalContrasts ==1) comArray[0] = synScript1Param.getValue().getAbsolutePath();
			else if (totalContrasts ==2) comArray[0] = synScript2Param.getValue().getAbsolutePath();
			else if (totalContrasts ==3) comArray[0] = synScript3Param.getValue().getAbsolutePath();
			else if (totalContrasts ==4) comArray[0] = synScript4Param.getValue().getAbsolutePath();
			int count = 1;

			comArray[1] = inputSource1.getAbsolutePath();
			comArray[2] = inputTarget1.getAbsolutePath();
			comArray[3] = Float.toString(1.0f);
			count += 3;

			if (scurv != null) {
				comArray[4] = inputSource2.getAbsolutePath();
				comArray[5] = inputTarget2.getAbsolutePath();
				comArray[6] = Float.toString(wcurvParam.getValue().floatValue()/2.0f);
				//comArray[6] = Float.toString(wcurvParam.getValue().floatValue()*(nbscales-(sc/2.0f))/nbscales);
				comArray[7] = inputSource3.getAbsolutePath();
				comArray[8] = inputTarget3.getAbsolutePath();
				comArray[9] = Float.toString(wcurvParam.getValue().floatValue()/2.0f);
				//comArray[9] = Float.toString(wcurvParam.getValue().floatValue()*(nbscales-(sc/2.0f))/nbscales);
				count += 6;
			}

			for (int n=0; n<nimg; n++) {
				comArray[count] = inputSource4[n].getAbsolutePath();
				count++;
				comArray[count] = inputTarget4[n].getAbsolutePath();
				count++;
				comArray[count] = Float.toString(wintensParam.getValue().floatValue()/nimg);
				count++;
			}

			comArray[count] = inputCoord.getAbsolutePath();
			comArray[count+1] = Integer.toString(coarseIt);
			comArray[count+2] = Integer.toString(medIt);
			comArray[count+3] = Integer.toString(fineIt);
			comArray[count+4] = Boolean.toString(affine);
			comArray[count+5] = outputDeform1.getAbsolutePath();
			count += 6;

			if (scurv != null) {
				comArray[count] = outputDeform2.getAbsolutePath();
				comArray[count+1] = outputDeform3.getAbsolutePath();
				count += 2;
			}

			for (int n=0; n<nimg; n++) {
				comArray[count] = outputDeform4[n].getAbsolutePath();
				count++;
			}

			comArray[count] = outputMapping.getAbsolutePath();
			comArray[count+1] = outputInvMapping.getAbsolutePath();
			comArray[count+2] = destdir.getAbsolutePath();
			count += 3;

			String[] comArray2 = new String[count];

			System.out.println("Current Command:");
			for(int i=0; i<comArray2.length; i++) {
				comArray2[i] = comArray[i];
				System.out.println(comArray2[i]);
			}
			System.out.println("\n");

			// run the script
			Runtime rt = Runtime.getRuntime();
			Process pr = rt.exec(comArray2);

			// collect the command line output
			BufferedReader inputOut = new BufferedReader(new InputStreamReader(pr.getInputStream()));
			BufferedReader inputError = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
			String lineOut=null;
			String lineError=null;
			while((lineOut=inputOut.readLine()) != null || (lineError=inputError.readLine()) != null) {
				System.out.println(lineOut);
				System.err.println(lineError);
			}
			int exitVal = pr.waitFor();

			if(exitVal==0){
				System.out.println("Executed successfully");
			}else{
				System.out.println("Executed with error(s)");
			}

			// clean up : load the result images and update mapping with new result
			ImageData mapImg = readerWriter.read(outputMapping);//convert output into ImageData;
			float[][][][] deformation = (new ImageDataFloat(mapImg)).toArray4d();
			float[][] composed = new float[3][nxyz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (deformation[x][y][z][X]>=1 && deformation[x][y][z][Y]>=1 && deformation[x][y][z][Z]>=1) {
					int xyz = x + nx*y + nx*ny*z;
					composed[X][xyz] = ImageInterpolation.linearClosestInterpolation(map[X], deformation[x][y][z][X], deformation[x][y][z][Y], deformation[x][y][z][Z], nx, ny, nz);
					composed[Y][xyz] = ImageInterpolation.linearClosestInterpolation(map[Y], deformation[x][y][z][X], deformation[x][y][z][Y], deformation[x][y][z][Z], nx, ny, nz);
					composed[Z][xyz] = ImageInterpolation.linearClosestInterpolation(map[Z], deformation[x][y][z][X], deformation[x][y][z][Y], deformation[x][y][z][Z], nx, ny, nz);
				}
			}
			for (int xyz=0;xyz<nxyz;xyz++) {
				map[X][xyz] = composed[X][xyz];
				map[Y][xyz] = composed[Y][xyz];
				map[Z][xyz] = composed[Z][xyz];
			}				
			composed = null;
			mapImg = null;
			deformation = null;

			ImageData invmapImg = readerWriter.read(outputInvMapping);//convert output into ImageData;
			deformation = (new ImageDataFloat(invmapImg)).toArray4d();
			composed = new float[3][nxyz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (deformation[x][y][z][X]>=1 && deformation[x][y][z][Y]>=1 && deformation[x][y][z][Z]>=1) {
					int xyz = x + nx*y + nx*ny*z;
					composed[X][xyz] = ImageInterpolation.linearClosestInterpolation(invmap[X], deformation[x][y][z][X], deformation[x][y][z][Y], deformation[x][y][z][Z], nx, ny, nz);
					composed[Y][xyz] = ImageInterpolation.linearClosestInterpolation(invmap[Y], deformation[x][y][z][X], deformation[x][y][z][Y], deformation[x][y][z][Z], nx, ny, nz);
					composed[Z][xyz] = ImageInterpolation.linearClosestInterpolation(invmap[Z], deformation[x][y][z][X], deformation[x][y][z][Y], deformation[x][y][z][Z], nx, ny, nz);
				}
			}
			for (int xyz=0;xyz<nxyz;xyz++) {
				invmap[X][xyz] = composed[X][xyz];
				invmap[Y][xyz] = composed[Y][xyz];
				invmap[Z][xyz] = composed[Z][xyz];
			}				
			composed = null;
			invmapImg = null;
			deformation = null;


			//removes all temporary files created for/by SyN
			if (!debugParam.getValue().booleanValue()) {
				ImageDataReaderWriter.deleteImageFile(inputSource1);
				ImageDataReaderWriter.deleteImageFile(inputTarget1);
				ImageDataReaderWriter.deleteImageFile(inputCoord);
				ImageDataReaderWriter.deleteImageFile(outputDeform1);

				if (scurv != null) {
					ImageDataReaderWriter.deleteImageFile(inputSource2);
					ImageDataReaderWriter.deleteImageFile(inputSource3);
					ImageDataReaderWriter.deleteImageFile(inputTarget2);
					ImageDataReaderWriter.deleteImageFile(inputTarget3);
					ImageDataReaderWriter.deleteImageFile(outputDeform2);
					ImageDataReaderWriter.deleteImageFile(outputDeform3);
				}	

				for (int n=0; n<nimg; n++) {
					ImageDataReaderWriter.deleteImageFile(inputSource4[n]);
					ImageDataReaderWriter.deleteImageFile(inputTarget4[n]);
					ImageDataReaderWriter.deleteImageFile(outputDeform4[n]);
				}

				ImageDataReaderWriter.deleteImageFile(outputMapping);
				ImageDataReaderWriter.deleteImageFile(outputInvMapping);
			}

		} catch(Exception e) {
			System.out.println("Executed with error(s)");
			e.printStackTrace();
		}

		return;
	}


	
}
