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

import edu.jhu.ece.iacl.algorithms.PrinceGroupAuthors;
import edu.jhu.ece.iacl.algorithms.ace.AnatomicallyConsistentEnhancement;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class JistCortexFullCRUISE extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume initImage;
	private ParamVolume wmImage;
	private ParamVolume gmImage;
	private ParamVolume csfImage;
	private ParamVolume vdImage;
	
	//private ParamFloat 	edgeParam;
	private	ParamFloat 	balloonParam;
	private ParamFloat 	curvParam;
	private ParamInteger 	iterationParam;
	//private ParamOption 	gwImageParam;
	//private ParamOption 	cgImageParam;
	private ParamOption 	topologyParam;
	private ParamBoolean	normalizeParam;
	private ParamBoolean	pvwmParam;
	//private ParamFloat 	offsetParam;
	private ParamFloat 		wmdistParam;
	
	//private static final String[] opts = {"Iso image", "T1 map", "Probabilities", "Iso+Proba", "T1map+Proba", "All"};
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume cortexImage;
	private ParamVolume gwbImage;
	private ParamVolume cgbImage;
	private ParamVolume aceImage;
	private ParamVolume avgImage;
	private ParamVolume pvwmImage;
	private ParamVolume thickImage;
	private ParamVolume pcsfImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(initImage = new ParamVolume("Initial WM Segmentation Image"));
		imageParams.add(wmImage = new ParamVolume("Filled WM Probability Image"));
		imageParams.add(gmImage = new ParamVolume("GM Probability Image"));
		imageParams.add(csfImage = new ParamVolume("CSF and BG Probability Image"));
		imageParams.add(vdImage = new ParamVolume("Veins and Dura Probability Image (opt)"));
		vdImage.setMandatory(false);
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		//mainParams.add(edgeParam = new ParamFloat("Edge weight", -1E10f, 1E10f, 0.1f25));
		mainParams.add(balloonParam = new ParamFloat("Data weight (balloon force)", -1E10f, 1E10f, 0.9f));
		mainParams.add(curvParam = new ParamFloat("Regularization weight (curvature)", -1E10f, 1E10f, 0.1f));
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
			
		//mainParams.add(gwImageParam = new ParamOption("Image used for WM / GM", opts));
		//gwImageParam.setValue("Iso+Proba");

		//mainParams.add(cgImageParam = new ParamOption("Image used for GM / CSF", opts));
		//cgImageParam.setValue("Iso+Proba");

		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		mainParams.add(normalizeParam = new ParamBoolean("Normalize probabilities", true));
		mainParams.add(pvwmParam = new ParamBoolean("Correct for WM-GM partial voluming", true));
		//mainParams.add(offsetParam = new ParamFloat("WM/GM offset", -1, 1, 0.0f));
		mainParams.add(wmdistParam = new ParamFloat("WM drop-off distance", 0.1f, 100.0f, 1.0f));
		
		inputParams.add(mainParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing");
		inputParams.setLabel("Full CRUISE Cortex Extraction");
		inputParams.setName("FullCruiseCortexExtraction");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Xiao Han", "","http://www.iacl.ece.jhu.edu/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new Citation("X. Han, D.L. Pham, D. Tosun, M.E. Rettmann, C. Xu, and J. L. Prince, "
								+"CRUISE: Cortical Reconstruction Using Implicit Surface Evolution, "
								+"NeuroImage, vol. 23, pp. 997--1012, 2004."));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Segments the cortex from a whole brain segmented data set with the CRUISE method \n"
								+"(includes partial voluming corrections and ACE sulcal enhancement).");
		
		info.setVersion("3.1.2");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(cortexImage = new ParamVolume("Cortex Mask",VoxelType.UBYTE));
		outputParams.add(gwbImage = new ParamVolume("WM-GM Levelset",VoxelType.FLOAT));
		outputParams.add(cgbImage = new ParamVolume("GM-CSF Levelset",VoxelType.FLOAT));
		outputParams.add(avgImage = new ParamVolume("Central Levelset",VoxelType.FLOAT));
		outputParams.add(thickImage = new ParamVolume("Cortical thickness",VoxelType.FLOAT));
		outputParams.add(pvwmImage = new ParamVolume("Cerebral WM probability",VoxelType.FLOAT));
		pvwmImage.setMandatory(false);
		outputParams.add(aceImage = new ParamVolume("Cortical GM probability",VoxelType.FLOAT));
		outputParams.add(pcsfImage = new ParamVolume("Sulcal CSF probability",VoxelType.FLOAT));
		
		outputParams.setName("cortical images");
		outputParams.setLabel("cortical images");
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
				
		ImageDataFloat csfImg = new ImageDataFloat(csfImage.getImageData());
		float[][][]csf = csfImg.toArray3d();
		csfImg = null;
		
		float[][][] vd = null;
		if (vdImage.getImageData()!=null) {
			ImageDataFloat vdImg = new ImageDataFloat(vdImage.getImageData());
			vd = vdImg.toArray3d();
			vdImg = null;
		}
		
		// 0. Create combined probabilities
		if (normalizeParam.getValue().booleanValue()) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				float pwm = 0.5f + 0.5f*(wm[x][y][z]-Numerics.max(gm[x][y][z],csf[x][y][z]));
				float pgm = 0.5f + 0.5f*(gm[x][y][z]-Numerics.max(wm[x][y][z],csf[x][y][z]));
				float pcsf = 0.5f + 0.5f*(csf[x][y][z]-Numerics.max(gm[x][y][z],wm[x][y][z]));
				
				wm[x][y][z] = pwm;
				gm[x][y][z] = pgm;
				csf[x][y][z] = pcsf;
			}
		}
		
		
		// 1. Run partial volume correction on WM membership (opt)
		float[][][] wmpv = null;;
		//float offset = offsetParam.getValue().floatValue();
		float offset = 0.0f;
		if (pvwmParam.getValue().booleanValue()) {
			wmpv = new float[nx][ny][nz];
			
			float avg = 0.0f;
			float sig = 0.0f;
			float den = 0.0f;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(wm,x,y,z,2)) {
				// check for zero-valued neighbors as well
				wmpv[x][y][z] = 0.0f;
				float best = 0.0f;
				float dmax = 13;
				for (int d=0;d<dmax;d++) {
					float val = ImageFilters.planeScore(wm, x,y,z,d);
					if (val*val>best*best) best = val;
				}
				wmpv[x][y][z] = best;
				avg += best;
				den++;
			}
			avg /= den;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(wm,x,y,z,2)) {
				sig += (wmpv[x][y][z]-avg)*(wmpv[x][y][z]-avg);
				
				if (wmpv[x][y][z]<0) wmpv[x][y][z] = 0.0f;
			}
			sig /= (den-1.0f);
			
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(wm,x,y,z,2)) {
				if (wmpv[x][y][z]>0) {
					wmpv[x][y][z] = 1.0f - (float)Math.exp( -0.5f*(wmpv[x][y][z]-avg)*(wmpv[x][y][z]-avg)/sig );
					if (wmpv[x][y][z]>wm[x][y][z]) {
						wm[x][y][z] = wmpv[x][y][z];
						// modify also the other regions?
						gm[x][y][z] *= wm[x][y][z]/wmpv[x][y][z];
						csf[x][y][z] *= wm[x][y][z]/wmpv[x][y][z];
					}
					/*
					// cube and multiply by offset (arbitrary!)
					wmpv[x][y][z] = wmpv[x][y][z]*wmpv[x][y][z]*wmpv[x][y][z]*offset;
					*/
				}
			}
		
			// intermediate result: PVWM- modified WM probability
			ImageDataFloat pvwmData = new ImageDataFloat(wm);		
			pvwmData.setHeader(gmImage.getImageData().getHeader());
			pvwmData.setName(gmImage.getImageData().getName()+"_pwm");
			pvwmImage.setValue(pvwmData);
			pvwmData = null;
			//wmpv = null;
		}
		
		// 2. Run ACE on CSF

		// pre-process: scale into [0,254]
		// and set the lowest value to zero (so that sum=1)
		float[][][] wmace = new float[nx][ny][nz];
		float[][][] gmace = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			/*  no normalization anymore?
			float sum = wm[x][y][z]+gm[x][y][z]+csf[x][y][z];
			if (sum>0.0f001) {
				wm[x][y][z] *= 254/sum;
				gm[x][y][z] *= 254/sum;
				csf[x][y][z] *= 254/sum;
			}
			*/
			wm[x][y][z] *= 254;
			gm[x][y][z] *= 254;
			csf[x][y][z] *= 254;
			
			/*
			if (wm[x][y][z]<gm[x][y][z] && wm[x][y][z]<csf[x][y][z]) wmace[x][y][z] = 0.0f;
			else wmace[x][y][z] = Numerics.bounded(wm[x][y][z], 0.0f, 254.0f);
			
			if (gm[x][y][z]<wm[x][y][z] && gm[x][y][z]<csf[x][y][z]) gmace[x][y][z] = 0.0f;
			else gmace[x][y][z] = Numerics.bounded(gm[x][y][z], 0.0f, 254.0f);
			*/
			// ACE: maximum separation between WM, CSF? better not to normalize...
			wmace[x][y][z] = Numerics.bounded(wm[x][y][z], 0.0f, 254.0f);
			gmace[x][y][z] = Numerics.bounded(gm[x][y][z], 0.0f, 254.0f);
			/* not good for ACE, but useful later on
			wmace[x][y][z] = Numerics.bounded(254.0f*wm[x][y][z]/(wm[x][y][z] + gm[x][y][z] + csf[x][y][z]), 0.0f, 254.0f);
			float csface = Numerics.bounded(254.0f*csf[x][y][z]/(wm[x][y][z] + gm[x][y][z] + csf[x][y][z]), 0.0f, 254.0f);
			gmace[x][y][z] = Numerics.bounded(254.0f-Numerics.max(wmace[x][y][z],csface), 0.0f, 254.0f);
			*/
		}
		
		// get WM mask
		ImageDataUByte initImg = new ImageDataUByte(initImage.getImageData());
		byte[][][] bytebuffer = initImg.toArray3d();
		byte[] init = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			init[xyz] = bytebuffer[x][y][z];
		}
		initImg = null;
		bytebuffer = null;
		// compute distance factor
		float wmdist = wmdistParam.getValue().floatValue();
		float[] factor = new float[nxyz];
		float[] update = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (init[xyz]>0) factor[xyz] = 1.0f;
		}
		float scaling = 1.0f/(1.0f+Numerics.square(1.0f/wmdist));
		for (int t=0;t<2*wmdist;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (factor[xyz]==0) {
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
						if (factor[xyz+i+j*nx+l*nx*ny]>0) update[xyz] = Numerics.max(update[xyz], scaling*factor[xyz+i+j*nx+l*nx*ny]);
					}
				}
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (update[xyz]>0) factor[xyz] = update[xyz];
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			wmace[x][y][z] *= factor[xyz];
		}

		// main algorithm
		AnatomicallyConsistentEnhancement ace = new AnatomicallyConsistentEnhancement();
		monitor.observe(ace);
		ImageDataFloat wmimg = new ImageDataFloat(wmace);
		ImageDataFloat gmimg = new ImageDataFloat(gmace);
		ace.solve(wmimg, gmimg, 126.9f);
		wmimg = null; gmimg = null;
		
		// change and rescale
		ImageDataFloat aceGmImg = (ImageDataFloat)ace.getEnhancedGM();
		float[][][] acegm = aceGmImg.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			csf[x][y][z] += gmace[x][y][z]-acegm[x][y][z];
			
			wm[x][y][z] /= 254;
			acegm[x][y][z] /= 254;
			csf[x][y][z] /= 254;
		}
		ace = null;
		aceGmImg = null;
		wmace = null;
		gmace = null;
		gm = null;
		
		// 3. run CRUISE
		float[] pwm = new float[nxyz];
		float[] pgm = new float[nxyz];
		float[] pcsf = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			// WM offset for segmentation of WM/GM interface
			float wmr = (float)FastMath.sin(Numerics.abs(pwm[xyz]-0.5f)*FastMath.PI);
			
			pwm[xyz] = Numerics.bounded(wm[x][y][z],0,1);
			//pgm[xyz] = Numerics.bounded(acegm[x][y][z]-offset,0,1);
			// offset: decreased in regions where pWM is close to 0.5f
			//float wmr = Numerics.square((pwm[xyz]-0.5f)/0.5f);
			
			pgm[xyz] = Numerics.bounded(acegm[x][y][z]-wmr*offset,0,1);
			acegm[x][y][z] = Numerics.bounded(acegm[x][y][z], 0,1);
			/*
			if (wmpv!=null) pgm[xyz] = Numerics.bounded(acegm[x][y][z]-wmpv[x][y][z],0,1);
			else pgm[xyz] = Numerics.bounded(acegm[x][y][z],0,1);
			*/
			pcsf[xyz] = Numerics.bounded(csf[x][y][z],0,1);
			/*
			// normalize?
			float sum = pgm[xyz]+pwm[xyz]+pcsf[xyz];
			if (sum>0.001) {
				pgm[xyz] /= sum;
				pwm[xyz] /= sum;
				pcsf[xyz] /= sum;
			}
			*/
		}
		//wm = null;
		//csf = null;
		
		// intermediate result: ACE- and PVWM- modified GM probability
		ImageDataFloat acegmData = new ImageDataFloat(acegm);		
		acegmData.setHeader(gmImage.getImageData().getHeader());
		acegmData.setName(gmImage.getImageData().getName()+"_pgm");
		aceImage.setValue(acegmData);
		acegmData = null;
		//acegm = null;
		
		// main algorithm
		boolean useProbas = true;
		boolean gwbflag = true;
		
		// levelset evolution: wm /gm		
		Mp2rageCortexGdm gdm = new Mp2rageCortexGdm(init, null, nx, ny, nz, rx, ry, rz,
														pwm, pgm, pcsf, null,
														null, null,
														0.0f, balloonParam.getValue().floatValue(), 
														curvParam.getValue().floatValue(), 0.0f,
														topologyParam.getValue(),null,useProbas, gwbflag);
		
		if (iterationParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			
			gdm.evolveNarrowBand(iterationParam.getValue().intValue(), 0.001f);
		}
		
		// output
		byte[][][] seg = new byte[nx][ny][nz];
		float[][][] gwb = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[x][y][z] = gdm.getSegmentation()[xyz];
			gwb[x][y][z] = gdm.getLevelSet()[xyz];
		}
		
		// second step

		// use the wm result as the filler mask => guarantees the boundaries never cross
		init = gdm.getSegmentation();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			pgm[xyz] = Numerics.bounded(acegm[x][y][z],0,1);
			if (vd!=null) csf[x][y][z] = Numerics.max(csf[x][y][z],vd[x][y][z]);
			csf[x][y][z] = Numerics.bounded(csf[x][y][z],0,1);
			pcsf[xyz] = csf[x][y][z];
			// remove WM proba far from current boundary (to help against dura mater problems)
			if (gwb[x][y][z]>0) pwm[xyz] = Numerics.bounded(wm[x][y][z],0,1)/(1.0f+Numerics.square(gwb[x][y][z]/wmdist));
			else pwm[xyz] = Numerics.bounded(wm[x][y][z],0,1);
			
			// merge WM and GM
			pgm[xyz] = Numerics.bounded(pgm[xyz]+pwm[xyz],0,1);
			/*
			// normalize
			float sum = pgm[xyz]+pwm[xyz]+pcsf[xyz];
			if (sum>0.001) {
				pgm[xyz] /= sum;
				pwm[xyz] /= sum;
				pcsf[xyz] /= sum;
			}
			*/
		}
		boolean cgbflag = !gwbflag;
		
		// intermediate result: ACE- and Vessel/Dura- modified CSF probability
		ImageDataFloat pcsfData = new ImageDataFloat(csf);		
		pcsfData.setHeader(gmImage.getImageData().getHeader());
		pcsfData.setName(gmImage.getImageData().getName()+"_pcsf");
		pcsfImage.setValue(pcsfData);
		pcsfData = null;
		//acegm = null;
		
		// levelset evolution: gm / csf
		gdm = new Mp2rageCortexGdm(init, null, nx, ny, nz, rx, ry, rz,
									pwm, pgm, pcsf, null,
									null, null,
									0.0f, balloonParam.getValue().floatValue(), 
									curvParam.getValue().floatValue(), 0.0f,
									topologyParam.getValue(),null, useProbas, cgbflag);
		
		if (iterationParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			
			gdm.evolveNarrowBand(iterationParam.getValue().intValue(), 0.001f);
		}

		float[][][] cgb = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[x][y][z] = (byte)(gdm.getSegmentation()[xyz]+seg[x][y][z]);
			cgb[x][y][z] = gdm.getLevelSet()[xyz];
		}
		gdm.finalize(); gdm = null;
		
		// estimate simple cortical thickness
		float[][][] thick = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (gwb[x][y][z]>=0 && cgb[x][y][z]<=0) {
				thick[x][y][z] = gwb[x][y][z]-cgb[x][y][z];
			} else {
				thick[x][y][z] = 0.0f;
			}
		}
		
		// 4. generate central surface
		boolean[] bgmask = new boolean[nxyz];
		float[] gwb2 = new float[nxyz];
		float[] ctr2 = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			// only do the masking for the outside
			//bgmask[xyz] = (gwb[x][y][z]>=-5.0f && cgb[x][y][z]<=5.0f);
			bgmask[xyz] = (cgb[x][y][z]<=5.0f);
			gwb2[xyz] = gwb[x][y][z];
			ctr2[xyz] = 0.5f*(gwb[x][y][z]+cgb[x][y][z]);
		}
		
		// main algorithm
		SmoothGdm gdm2 = new SmoothGdm(gwb2, ctr2, nx, ny, nz, rx, ry, rz,
										bgmask, balloonParam.getValue().floatValue(), 
										curvParam.getValue().floatValue(), topologyParam.getValue(), null);

		BasicInfo.displayMessage("average estimation...\n");
		gdm2.evolveNarrowBand(iterationParam.getValue().intValue(), 0.001f);
		float[][][] avg = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			avg[x][y][z] = gdm2.getLevelSet()[xyz];
		}
		gdm2.finalize(); gdm2 = null;
		
		
		ImageDataUByte cortexData = new ImageDataUByte(seg);		
		cortexData.setHeader(gmImage.getImageData().getHeader());
		cortexData.setName(gmImage.getImageData().getName()+"_cortex");
		cortexImage.setValue(cortexData);
		cortexData = null;
		seg = null;
		
		ImageDataFloat gwbData = new ImageDataFloat(gwb);		
		gwbData.setHeader(gmImage.getImageData().getHeader());
		gwbData.setName(gmImage.getImageData().getName()+"_gwb");
		gwbImage.setValue(gwbData);
		gwbData = null;
		gwb = null;
		
		ImageDataFloat cgbData = new ImageDataFloat(cgb);		
		cgbData.setHeader(gmImage.getImageData().getHeader());
		cgbData.setName(gmImage.getImageData().getName()+"_cgb");
		cgbImage.setValue(cgbData);
		cgbData = null;
		cgb = null;
		
		ImageDataFloat avgData = new ImageDataFloat(avg);		
		avgData.setHeader(gmImage.getImageData().getHeader());
		avgData.setName(gmImage.getImageData().getName()+"_avg");
		avgImage.setValue(avgData);
		avgData = null;
		avg = null;
		
		ImageDataFloat thickData = new ImageDataFloat(thick);		
		thickData.setHeader(gmImage.getImageData().getHeader());
		thickData.setName(gmImage.getImageData().getName()+"_thick");
		thickImage.setValue(thickData);
		thickData = null;
		thick = null;
		
	}

	boolean zeroNeighbor(float[][][] image, int x, int y, int z, int d) {
		for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int l=-d;l<=d;l++) {
			if (image[x+i][y+j][z+l]!=image[x][y][z] && i*i+j*j+l*l<=2*d*d) return false;
		}
		return true;
	}
	

}
