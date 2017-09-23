package de.mpg.cbs.core.cortex;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class CortexOptimCRUISE {

	private byte[] initImage;
	private float[] wmImage;
	private float[] gmImage;
	private float[] csfImage;
	private float[] vdImage;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private	float 		balloonParam = 0.4f;
	private float 		curvParam = 0.1f;
	private int 		iterationParam = 500;
	private String 		topologyParam = "wcs";
	private boolean		normalizeParam = true;
	private boolean		pvwmParam = true;
	private float	 	wmdistParam = 1.0f;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	private String	lutdir = null;
	
	private byte[] cortexImage;
	private float[] gwbImage;
	private float[] cgbImage;
	private float[] aceImage;
	private float[] avgImage;
	private float[] pwmImage;
	private float[] thickImage;
	private float[] pcsfImage;
	
	public final void setInitialWMSegmentationImage(byte[] val) { initImage = val; }
	public final void importInitialWMSegmentationImage(int[] val) { 
	    initImage = new byte[val.length];
	    for (int xyz=0;xyz<val.length;xyz++) initImage[xyz] = (byte)val[xyz];
	}
	public final void setFilledWMProbabilityImage(float[] val) { wmImage = val; }
	public final void setGMProbabilityImage(float[] val) { gmImage = val; }
	public final void setCSFandBGProbabilityImage(float[] val) { csfImage = val; }
	public final void setVeinsAndDuraProbabilityImage(float[] val) { vdImage = val; }
	
	public final void setDataWeight(float val) { balloonParam = val; }
	public final void setRegularizationWeight(float val) { curvParam = val; }
	public final void setMaxIterations(int val) { iterationParam = val; }
	public final void setNormalizeProbabilities(boolean val) { normalizeParam = val; }
	public final void setCorrectForWMGMpartialVoluming(boolean val) { pvwmParam = val; }
	public final void setWMdropoffDistance(float val) { wmdistParam = val; }
		
	public final void setTopology(String val) { topologyParam = val; }
	public final void setTopologyLUTdirectory(String val) { lutdir = val; }
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public  final String getPackage() { return "CBS Tools"; }
	public  final String getCategory() { return "Cortex Processing"; }
	public  final String getLabel() { return "CRUISE Cortex Extraction"; }
	public  final String getName() { return "CRUISECortexExtraction"; }
	public  final String getVersion() { return "3.1.2"; }

	public  final String[] getAlgorithmAuthors() {return new String[] {"Xiao Han","Pierre-Louis Bazin"}; }
	public  final String getCitation() { return "X. Han, D.L. Pham, D. Tosun, M.E. Rettmann, C. Xu, and J. L. Prince, "
												+"CRUISE: Cortical Reconstruction Using Implicit Surface Evolution, "
												+"NeuroImage, vol. 23, pp. 997--1012, 2004."; }
	public  final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public  final String getDescription() { return "Segments the cortex from a whole brain segmented data set with the CRUISE method \n"
													+"(includes customized partial voluming corrections and ACE sulcal enhancement)."; }
	public  final String getLongDescription() { return getDescription(); }

	// create outputs
	public final byte[] getCortexMask() { return cortexImage; }
	public final float[] getWMGMLevelset() { return gwbImage; }
	public final float[] getGMCSFLevelset() { return cgbImage; }
	public final float[] getCentralLevelset() { return avgImage; }
	public final float[] getCorticalThickness() { return thickImage; }
	public final float[] getCerebralWMprobability() { return pwmImage; }
	public final float[] getCorticalGMprobability() { return aceImage; }
	public final float[] getSulcalCSFprobability() { return pcsfImage; }
		

	public void execute(){
		
		// import the image data into 1D arrays
		float[] wm = wmImage;
		float[] gm = gmImage;
		float[] csf = csfImage;
		float[] vd = vdImage;
		byte[] init = initImage;
		
		// 0. Create combined probabilities
		if (normalizeParam) {
			for (int xyz=0;xyz<nxyz;xyz++) {
				float pwm = 0.5f + 0.5f*(wm[xyz]-Numerics.max(gm[xyz],csf[xyz]));
				float pgm = 0.5f + 0.5f*(gm[xyz]-Numerics.max(wm[xyz],csf[xyz]));
				float pcsf = 0.5f + 0.5f*(csf[xyz]-Numerics.max(gm[xyz],wm[xyz]));
				
				wm[xyz] = pwm;
				gm[xyz] = pgm;
				csf[xyz] = pcsf;
			}
		}
		
		
		// 1. Run partial volume correction on WM membership (opt)
		float[] wmpv = null;;
		//float offset = offsetParam.getValue().floatValue();
		float offset = 0.0f;
		if (pvwmParam) {
			wmpv = new float[nxyz];
			boolean[] positive = new boolean[nxyz];
			for (int xyz=0;xyz<nxyz;xyz++) positive[xyz] = (wm[xyz]>0);
			wmpv = ImageFilters.surfaceEdgeFilter(wm, positive, nx,ny,nz);
			float avg = 0.0f;
			float sig = 0.0f;
			float den = 0.0f;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (wmpv[xyz]!=0) {
					avg += wmpv[xyz];
					den++;
				}
			}
			avg /= den;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (wmpv[xyz]!=0) {
					sig += (wmpv[xyz]-avg)*(wmpv[xyz]-avg);
				}
				if (wmpv[xyz]<0) wmpv[xyz] = 0.0f;
			}
			sig /= (den-1.0f);
			
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (wmpv[xyz]>0) {
					wmpv[xyz] = 1.0f - (float)Math.exp( -0.5f*(wmpv[xyz]-avg)*(wmpv[xyz]-avg)/sig );
					if (wmpv[xyz]>wm[xyz]) {
						wm[xyz] = wmpv[xyz];
						// modify also the other regions?
						gm[xyz] *= wm[xyz]/wmpv[xyz];
						csf[xyz] *= wm[xyz]/wmpv[xyz];
					}
				}
			}
		}
		
		// 2. Run ACE on CSF (new version)
		
		// compute distance factor for reducing impact of WM far from initial boundary / 0.5 probability (for dura mostly)
		boolean[] wmmask = new boolean[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (init[xyz]>0 || wm[xyz]>=0.5f) wmmask[xyz] = true;
			else wmmask[xyz] = false;
		}
		wmmask = ObjectLabeling.largestObject(wmmask, nx,ny,nz, 26);
		float wmdist = wmdistParam;
		float[] factor = new float[nxyz];
		float[] update = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (wmmask[xyz]) factor[xyz] = 1.0f;
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
		float[] wmace = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			wmace[xyz] = wm[xyz]*factor[xyz];
			//wm[xyz] *= factor[xyz];
		}

		// main algorithm
		AnatomicallyConsistentEnhancement ace = new AnatomicallyConsistentEnhancement(init,wmace,gm,csf,nx,ny,nz,topologyParam,lutdir);
		ace.computeSulcalRidges(0.8f);
		ace.thinSulcalSkeleton(true, 10);
		
		for (int xyz=0;xyz<nxyz;xyz++) {
			gm[xyz] = ace.getACEprobaGM()[xyz];
			csf[xyz] = ace.getACEprobaCSF()[xyz];
		}
		ace = null;
		
		// 3. run CRUISE
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			wm[xyz] = Numerics.bounded(wm[xyz],0,1);
			gm[xyz] = Numerics.bounded(gm[xyz],0,1);
			csf[xyz] = Numerics.bounded(csf[xyz],0,1);
		}
		
		// main algorithm
		boolean useProbas = true;
		boolean gwbflag = true;
		
		// levelset evolution: wm /gm		
		Mp2rageCortexGdm gdm = new Mp2rageCortexGdm(init, null, nx, ny, nz, rx, ry, rz,
														wm, gm, csf, null,
														null, null,
														0.0f, balloonParam, 
														curvParam, 0.0f,
														topologyParam,lutdir,useProbas, gwbflag);
		
		if (iterationParam>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			
			gdm.evolveNarrowBand(iterationParam, 0.001f);
		}
		
		// output
		byte[] seg = new byte[nxyz];
		float[] gwb = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[xyz] = gdm.getSegmentation()[xyz];
			gwb[xyz] = gdm.getLevelSet()[xyz];
		}
		
		// second step

		// use the wm result as the filler mask => guarantees the boundaries never cross
		init = gdm.getSegmentation();
		float[] wmgm = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (vd!=null) csf[xyz] = Numerics.bounded(Numerics.max(csf[xyz],vd[xyz]),0,1);
			// remove WM proba far from current boundary (to help against dura mater problems)
			if (gwb[xyz]>0) wm[xyz] = Numerics.bounded(wm[xyz],0,1)/(1.0f+Numerics.square(gwb[xyz]/wmdist));
			
			// merge WM and GM
			wmgm[xyz] = Numerics.bounded(gm[xyz]+wm[xyz],0,1);
		}
		boolean cgbflag = !gwbflag;
		
		// levelset evolution: gm / csf
		gdm = new Mp2rageCortexGdm(init, null, nx, ny, nz, rx, ry, rz,
									wm, wmgm, csf, null,
									null, null,
									0.0f, balloonParam, 
									curvParam, 0.0f,
									topologyParam,lutdir,useProbas, cgbflag);
		
		if (iterationParam>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			
			gdm.evolveNarrowBand(iterationParam, 0.001f);
		}

		float[] cgb = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[xyz] = (byte)(gdm.getSegmentation()[xyz]+seg[xyz]);
			cgb[xyz] = gdm.getLevelSet()[xyz];
		}
		gdm.finalize(); gdm = null;
		
		// estimate simple cortical thickness
		float[] thick = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (gwb[xyz]>=0 && cgb[xyz]<=0) {
				thick[xyz] = gwb[xyz]-cgb[xyz];
			} else {
				thick[xyz] = 0.0f;
			}
		}
		
		// 4. generate central surface? only for display
		boolean[] bgmask = new boolean[nxyz];
		float[] gwb2 = new float[nxyz];
		float[] ctr2 = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			// only do the masking for the outside
			//bgmask[xyz] = (gwb[x][y][z]>=-5.0f && cgb[x][y][z]<=5.0f);
			bgmask[xyz] = (cgb[xyz]<=5.0f);
			gwb2[xyz] = gwb[xyz];
			ctr2[xyz] = 0.5f*(gwb[xyz]+cgb[xyz]);
		}
		
		// main algorithm
		SmoothGdm gdm2 = new SmoothGdm(gwb2, ctr2, nx, ny, nz, rx, ry, rz,
										bgmask, balloonParam, 
										curvParam, topologyParam, lutdir);

		BasicInfo.displayMessage("average estimation...\n");
		gdm2.evolveNarrowBand(iterationParam, 0.001f);
		float[] avg = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			avg[xyz] = gdm2.getLevelSet()[xyz];
		}
		gdm2.finalize(); gdm2 = null;
		
		// outputs
		cortexImage = seg;
		gwbImage = gwb;
		cgbImage = cgb;
		avgImage = avg;
		thickImage = thick;
		
		pwmImage = wm;
		aceImage = gm;
		pcsfImage = csf;
		
		return;
	}

	boolean zeroNeighbor(float[] image, int x, int y, int z, int d) {
		int xyz = x+nx*y+nx*ny*z;
		for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int l=-d;l<=d;l++) {
			if (image[xyz+i+nx*j+nx*ny*l]!=image[xyz] && i*i+j*j+l*l<=2*d*d) return false;
		}
		return true;
	}
	
}
