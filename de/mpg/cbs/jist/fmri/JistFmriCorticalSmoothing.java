package de.mpg.cbs.jist.fmri;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
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

import gov.nih.mipav.model.file.FileInfoBase;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.model.structures.ModelStorageBase;

/*
 * @author Pierre-Louis Bazin
 */
public class JistFmriCorticalSmoothing extends ProcessingAlgorithm {

	// jist containers
	
	// input 
	private ParamVolume 	gwblevelsetImage;
	private ParamVolume 	ctrlevelsetImage;
	private ParamVolume 	cgblevelsetImage;
	private ParamVolume 	dataImage;
	private ParamVolume 	mapImage;
	
	private ParamFloat		fwhmParam;
	private ParamOption		regionParam;
	private static final String[] regions = {"full_cortex","central_surface","exclude_pv","inner_band","central_band","outer_band"};
	private ParamFloat		rangeParam;
	private ParamOption		spaceParam;
	private static final String[] spaces = {"anatomical","functional"};
	//private static final String[] spaces = {"anatomical","subsampled_anatomical","functional","oversampled_functional"};
	private ParamOption		typeParam;
	private static final String[] types = {"volumetric","anatomy-guided"};
	private ParamOption		interpParam;
	private static final String[] interps = {"linear","nearest"};
	//private ParamOption		subfactorParam;
	//private ParamOption		overfactorParam;
	//private static final String[] factors = {"1","2","3","4"};
	
	private static final byte X=0;
	private static final byte Y=1;
	private static final byte Z=2;
	
	// ouput
	private ParamVolume smoothdataImage;
	private ParamVolume smoothmaskImage;
		
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(dataImage = new ParamVolume("fMRI Data Image",VoxelType.FLOAT,-1,-1,-1,-1));
		dataImage.setLoadAndSaveOnValidate(false);
		inputParams.add(gwblevelsetImage = new ParamVolume("GM-WM Surface Levelset"));
		gwblevelsetImage.setLoadAndSaveOnValidate(false);
		inputParams.add(cgblevelsetImage = new ParamVolume("CSF-GM Surface Levelset"));
		cgblevelsetImage.setLoadAndSaveOnValidate(false);
		inputParams.add(ctrlevelsetImage = new ParamVolume("Central Surface Levelset (opt)"));
		ctrlevelsetImage.setMandatory(false);
		ctrlevelsetImage.setLoadAndSaveOnValidate(false);
		inputParams.add(mapImage = new ParamVolume("Anatomy to Function Voxel Mapping (opt)",VoxelType.FLOAT,-1,-1,-1,-1));
		mapImage.setMandatory(false);
		mapImage.setLoadAndSaveOnValidate(false);
				
		inputParams.add(fwhmParam = new ParamFloat("Gaussian FWHM (mm)", 0.0f, 50.0f, 5.0f));
		inputParams.add(typeParam = new ParamOption("Smoothing mode", types));
		typeParam.setValue("anatomy-guided");
					
		inputParams.add(regionParam = new ParamOption("Computation region", regions));
		regionParam.setValue("full_cortex");
		inputParams.add(rangeParam = new ParamFloat("min boundary distance (mm)", 0.0f, 50.0f, 1.0f));
		
		inputParams.add(interpParam = new ParamOption("fMRI interpolation method", interps));
		interpParam.setValue("nearest");
		
		inputParams.add(spaceParam = new ParamOption("Output space", spaces));
		spaceParam.setValue("functional");
		
		//inputParams.add(subfactorParam = new ParamOption("anatomy sub-sampling factor", factors));
		//subfactorParam.setValue("1");
		//inputParams.add(overfactorParam = new ParamOption("functional over-sampling factor", factors));
		//overfactorParam.setValue("2");
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("fMRI");
		inputParams.setLabel("fMRI Cortical Smoothing");
		inputParams.setName("fMRICorticalSmoothing");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Smooths data using a Gaussian weighting as a function of the minimum marching distance in the cortex.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(smoothdataImage = new ParamVolume("Smoothed Data Image",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(smoothmaskImage = new ParamVolume("Smoothing Mask Image",VoxelType.BYTE));
		
		outputParams.setLabel("fMRI Cortical Smoothing");
		outputParams.setName("fMRICorticalSmoothing");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the volume data into 1D arrays : faster approach?
		// get the reader/writer from the ParamVolume and then pull the ModelImage directly from it?
		//ImageDataMipav	dataImg = (ImageDataMipav)dataImage.getImageData();
		
		ImageDataFloat	dataImg = new ImageDataFloat(dataImage.getImageData());
		int nfx = dataImg.getRows();
		int nfy = dataImg.getCols();
		int nfz = dataImg.getSlices();
		int nfd = dataImg.getComponents();
		int nfxyz = nfx*nfy*nfz;
		float rfx = dataImg.getHeader().getDimResolutions()[0];
		float rfy = dataImg.getHeader().getDimResolutions()[1];
		float rfz = dataImg.getHeader().getDimResolutions()[2];
		int xyzf;

		String dataName = dataImg.getName();
		ImageHeader dataHeader = dataImg.getHeader();
				/*
		// from the mask, create an indexed 1D array for the fMRI (less memory usage => faster)
		int[][][] index = new int[nfx][nfy][nfz];
		

		// fMRI data
		float[] fmri = new float[nfx*nfy*nfz*nfd];
		ModelImage img = dataImg.getModelImageDirect();
		img.exportData(0, nfx*nfy*nfz*nfd, fmri);
		img.disposeLocal();
		img = null;
		dataImg.dispose();
		*/
		
		float[][][][] data = dataImg.toArray4d();
		dataImg.dispose();
		dataImage.dispose();
		
		System.out.println("fmri data loaded");

		// get anatomical information
		int nax,nay,naz;
		float rax,ray,raz;
		ImageDataFloat gwbImg=null, cgbImg=null, ctrImg=null;
		
		gwbImg = new ImageDataFloat(gwblevelsetImage.getImageData());
		cgbImg = new ImageDataFloat(cgblevelsetImage.getImageData());
		if (regionParam.getValue().equals("central_band") 
			|| regionParam.getValue().equals("inner_band") 
			|| regionParam.getValue().equals("outer_band") 
			|| regionParam.getValue().equals("central_surface"))
			ctrImg = new ImageDataFloat(ctrlevelsetImage.getImageData());
		
		nax = gwbImg.getRows();
		nay = gwbImg.getCols();
		naz = gwbImg.getSlices();
		rax = gwbImg.getHeader().getDimResolutions()[0];
		ray = gwbImg.getHeader().getDimResolutions()[1];
		raz = gwbImg.getHeader().getDimResolutions()[2];

		int naxyz = nax*nay*naz;
		int xyza;
		
		// distance to boundary
		float mindist = rangeParam.getValue().floatValue()/Numerics.min(rax,ray,raz);
		
		// create mask
		boolean[][][] msk = new boolean[nax][nay][naz];
		if (regionParam.getValue().equals("central_surface")) {
			float[][][] gwb = gwbImg.toArray3d();
			float[][][] cgb = cgbImg.toArray3d();
			float[][][] ctr = ctrImg.toArray3d();
			
			for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
				if (gwb[x][y][z]>=0 && cgb[x][y][z]<=0 && ctr[x][y][z]>-mindist && ctr[x][y][z]<mindist) {
					msk[x][y][z] = true;
				} else {
					msk[x][y][z] = false;
				}
			}
			gwb=null;
			cgb=null;
			ctr=null;
			
		} else if (regionParam.getValue().equals("central_band")) {
			float[][][] gwb = gwbImg.toArray3d();
			float[][][] cgb = cgbImg.toArray3d();
			float[][][] ctr = ctrImg.toArray3d();
			
			for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
				if (gwb[x][y][z]>=0 && cgb[x][y][z]<=0 && ctr[x][y][z]>=-gwb[x][y][z] && ctr[x][y][z]<=-cgb[x][y][z]) {
					msk[x][y][z] = true;
				} else {
					msk[x][y][z] = false;
				}
			}
			gwb=null;
			cgb=null;
			ctr=null;
			
		} else if (regionParam.getValue().equals("inner_band")) {
			float[][][] gwb = gwbImg.toArray3d();
			float[][][] cgb = cgbImg.toArray3d();
			float[][][] ctr = ctrImg.toArray3d();
			
			for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
				if (gwb[x][y][z]>=0 && ctr[x][y][z]<=0) {
					msk[x][y][z] = true;
				} else {
					msk[x][y][z] = false;
				}
			}
			gwb=null;
			cgb=null;
			ctr=null;
			
		} else if (regionParam.getValue().equals("outer_band")) {
			float[][][] gwb = gwbImg.toArray3d();
			float[][][] cgb = cgbImg.toArray3d();
			float[][][] ctr = ctrImg.toArray3d();
			
			for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
				if (ctr[x][y][z]>=0 && cgb[x][y][z]<=0) {
					msk[x][y][z] = true;
				} else {
					msk[x][y][z] = false;
				}
			}
			gwb=null;
			cgb=null;
			ctr=null;
			
		} else if (regionParam.getValue().equals("exclude_pv")) {
			float[][][] gwb = gwbImg.toArray3d();
			float[][][] cgb = cgbImg.toArray3d();
			
			for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
				if (gwb[x][y][z]>=mindist && cgb[x][y][z]<=-mindist) {
					msk[x][y][z] = true;
				} else {
					msk[x][y][z] = false;
				}
			}
			gwb=null;
			cgb=null;
		} else {
			float[][][] gwb = gwbImg.toArray3d();
			float[][][] cgb = cgbImg.toArray3d();
			
			for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
				if (gwb[x][y][z]>=0 && cgb[x][y][z]<=0) {
					msk[x][y][z] = true;
				} else {
					msk[x][y][z] = false;
				}
			}
			gwb=null;
			cgb=null;
		}
		gwbImg.dispose();
		cgbImg.dispose();
		if (ctrImg!=null) ctrImg.dispose();
		gwblevelsetImage.dispose();
		cgblevelsetImage.dispose();
		if (ctrlevelsetImage.getImageData()!=null) ctrlevelsetImage.dispose();
		System.out.println("mask created");
		
		// mapping
		float[][][][] mapping = null;
		if (mapImage.getImageData()!=null) {
			ImageDataFloat	mapImg = new ImageDataFloat(mapImage.getImageData());
			mapping = mapImg.toArray4d();
			mapImg.dispose();
		}
				
		// output variables
		ImageDataFloat sdataImg;
		float[][][][] sdata = null;
		float[][][] sden = null;
		//int substep = Integer.decode(subfactorParam.getValue()).intValue();
		//int overstep = Integer.decode(overfactorParam.getValue()).intValue();
		int substep = 1;
		int overstep = 1;
		if (spaceParam.getValue().equals("anatomical")) {
			sdata = new float[nax][nay][naz][nfd];
			sden = null;
			substep = 1;
			overstep = 1;
		} 
		else if (spaceParam.getValue().equals("subsampled_anatomical")) {
			sdata = new float[nax/substep][nay/substep][naz/substep][nfd];
			sden = null;
			overstep = 1;
		}
		else if (spaceParam.getValue().equals("functional")) {
			sdata = new float[nfx][nfy][nfz][nfd];
			sden = new float[nfx][nfy][nfz];
			overstep = 1;
		} 
		else if (spaceParam.getValue().equals("oversampled_functional")) {
			sdata = new float[overstep*nfx][overstep*nfy][overstep*nfz][nfd];
			sden = new float[overstep*nfx][overstep*nfy][overstep*nfz];
		}
			
		System.out.println("output data created");
		
		// main algorithm
		
		CorticalFmriSmoothing regionSmoothing = new CorticalFmriSmoothing(data, mapping, msk, fwhmParam.getValue().floatValue(), interpParam.getValue(),
																			nax, nay, naz, rax, ray, raz, nfx, nfy, nfz, nfd, rfx, rfy, rfz);
		
		System.out.println("smoothing loop");

		// for debug
		int pts=0;
		for (int x=0;x<nax-1;x++) for (int y=0;y<nay-1;y++) for (int z=0;z<naz-1;z++) {
			//xyz = x+nx*y+nx*ny*z;
			if (msk[x][y][z]) pts++;
		}
		int npt = 0;
		int iter = 0;
		long looptime = System.currentTimeMillis();
		float[] smoothed = new float[nfd];
		int xp,yp,zp;
		double alpha,beta,gamma;
		float weight;
		for (int x=0;x<nax-1;x+=substep) for (int y=0;y<nay-1;y+=substep) for (int z=0;z<naz-1;z+=substep) {
			//xyz = x+nx*y+nx*ny*z;
			
			if (msk[x][y][z]) {
				// no need to compute outside the cortex
				npt++;
				if (npt%(pts/100)==0) {
					iter++;
					long newtime = System.currentTimeMillis();
					BasicInfo.displayMessage("iter "+iter+", t="+(newtime-looptime)+", "+npt+" pts\n");
					looptime = newtime;
				}
				
				// minimal marching distance from voxel
				//sdata[x][y][z] = regionSmoothing.smoothFromPoint( xyz);
				if (typeParam.getValue().equals("volumetric")) {
					regionSmoothing.dilateVolumetrically(smoothed, x, y, z);
				} else {
					regionSmoothing.dilateFromPoint(smoothed, x, y, z);
				}
				
				if (spaceParam.getValue().equals("anatomical") || 
					spaceParam.getValue().equals("subsampled_anatomical") || 
					mapImage.getImageData()==null) {
					for (int n=0;n<nfd;n++) sdata[x/substep][y/substep][z/substep][n] = smoothed[n];
				} else {
					xp = Numerics.bounded(Numerics.floor(mapping[x][y][z][X]*overstep),0,nfx-2);	
					yp = Numerics.bounded(Numerics.floor(mapping[x][y][z][Y]*overstep),0,nfy-2);	
					zp = Numerics.bounded(Numerics.floor(mapping[x][y][z][Z]*overstep),0,nfz-2);
					alpha = Numerics.bounded(mapping[x][y][z][X]*overstep-xp,0,1);
					beta  = Numerics.bounded(mapping[x][y][z][Y]*overstep-yp,0,1);
					gamma = Numerics.bounded(mapping[x][y][z][Z]*overstep-zp,0,1);
					
					weight = (float)( (1.0f-alpha)*(1.0f-beta)*(1.0f-gamma) );
					for (int n=0;n<nfd;n++) sdata[xp][yp][zp][n] += weight*smoothed[n];
					sden[xp][yp][zp] += weight;
					weight = (float)( alpha*(1.0f-beta)*(1.0f-gamma) );
					for (int n=0;n<nfd;n++) sdata[xp+1][yp][zp][n] += weight*smoothed[n];
					sden[xp+1][yp][zp] += weight;
					weight = (float)( (1.0f-alpha)*beta*(1.0f-gamma) );
					for (int n=0;n<nfd;n++) sdata[xp][yp+1][zp][n] += weight*smoothed[n];
					sden[xp][yp+1][zp] += weight;
					weight = (float)( (1.0f-alpha)*(1.0f-beta)*gamma );
					for (int n=0;n<nfd;n++) sdata[xp][yp][zp+1][n] += weight*smoothed[n];
					sden[xp][yp][zp+1] += weight;
					weight = (float)( alpha*beta*(1.0f-gamma) );
					for (int n=0;n<nfd;n++) sdata[xp+1][yp+1][zp][n] += weight*smoothed[n];
					sden[xp+1][yp+1][zp] += weight;
					weight = (float)( (1.0f-alpha)*beta*gamma );
					for (int n=0;n<nfd;n++) sdata[xp][yp+1][zp+1][n] += weight*smoothed[n];
					sden[xp][yp+1][zp+1] += weight;
					weight = (float)( alpha*(1.0f-beta)*gamma );
					for (int n=0;n<nfd;n++) sdata[xp+1][yp][zp+1][n] += weight*smoothed[n];
					sden[xp+1][yp][zp+1] += weight;
					weight = (float)( alpha*beta*gamma );
					for (int n=0;n<nfd;n++) sdata[xp+1][yp+1][zp+1][n] += weight*smoothed[n];
					sden[xp+1][yp+1][zp+1] += weight;
				}
			}	
		}
		data = null;
		mapping = null;
		
		// normalize
		if (!spaceParam.getValue().equals("anatomical") 
				&& !spaceParam.getValue().equals("subsampled_anatomical")
					&& mapImage.getImageData()!=null) {
			for (int x=0;x<nfx*overstep;x++) for (int y=0;y<nfy*overstep;y++) for (int z=0;z<nfz*overstep;z++) {
				if (sden[x][y][z]>0) {
					for (int n=0;n<nfd;n++) sdata[x][y][z][n] /= sden[x][y][z];
				}
			}
		}

		
		// output mask
		byte[][][] outmask;
		if (spaceParam.getValue().equals("anatomical") 
			|| spaceParam.getValue().equals("subsampled_anatomical") 
			|| mapImage.getImageData()==null) {
			outmask = new byte[nax/substep][nay/substep][naz/substep];
			for (int x=0;x<nax;x+=substep) for (int y=0;y<nay;y+=substep) for (int z=0;z<naz;z+=substep) {
				if (msk[x][y][z]) outmask[x/substep][y/substep][z/substep] = 1;
			}
		} else {
			outmask = new byte[nfx*overstep][nfy*overstep][nfz*overstep];	
			for (int x=0;x<nfx*overstep;x++) for (int y=0;y<nfy*overstep;y++) for (int z=0;z<nfz*overstep;z++) {
				if (sden[x][y][z]>0) outmask[x][y][z] = 1;
			}
		}
		
		System.out.println("output..");
		
		// output
		sdataImg = new ImageDataFloat(sdata);		
		sdataImg.setHeader(dataHeader);
		sdataImg.setName(dataName + "_smooth");
		smoothdataImage.setValue(sdataImg);
		
		ImageDataByte smaskImg = new ImageDataByte(outmask);		
		smaskImg.setHeader(dataHeader);
		smaskImg.setName(dataName + "_smask");
		smoothmaskImage.setValue(smaskImg);
		
		// clean up
		sdataImg = null;
		sdata = null;
		regionSmoothing.finalize();
		
	}


}
