package de.mpg.cbs.core.brain;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class BrainMp2rageSkullStripping {

	// jist containers
	private float[] inv2Image;
	private float[] t1mapImage;
	private float[] isoImage;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private		static final byte	HNORM = 101;
	private		static final byte	EXP = 102;
	
	private boolean	skip0Param = false;
	private String	lutdir = null;
	
	private byte[] brainmaskImage;
	private float[] t1maskImage;
	private float[] isomaskImage;
	private float[] inv2maskImage;
	
	// create inputs
	public final void setSecondInversionImage(float[] val) { inv2Image = val; }
	public final void setT1MapImage(float[] val) { t1mapImage = val; }
	public final void setT1weightedImage(float[] val) { isoImage = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setSkipZeroValues(boolean val) { skip0Param = val; }
	public final void setTopologyLUTdirectory(String val) { lutdir = val; }

	// to be used for JIST definitions, generic info / help
	public static final String getPackage() { return "CBS Tools"; }
	public static final String getCategory() { return "Brain Processing"; }
	public static final String getLabel() { return "MP2RAGE Skull Stripping"; }
	public static final String getName() { return "Mp2rageSkullStripping"; }

	public final String[] getAlgorithmAuthors() {return new String[] {"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Estimate a brain mask for a MP2RAGE dataset. At least a T1-weighted or a T1 map image is required."; }
		
	public static final String getVersion() { return "3.1"; }

	// create outputs
	public final byte[] getBrainMaskImage() { return brainmaskImage; }
	public final float[] getMaskedT1MapImage() { return t1maskImage; }
	public final float[] getMaskedT1weightedImage() { return isomaskImage; }
	public final float[] getMaskedSecondInversionImage() { return inv2maskImage; }
		
	public void execute(){
		
		// main algorithm

		// assume  an exponential + outlier distributions for the inv2 image
		
		// use default parameters: 
		float maxdiff = 0.0001f;
		int itermax = 20;
		
		// min max for t1 map fixing
		float max = -1e10f;
		float min = 1e10f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (inv2Image[xyz]>max) max = inv2Image[xyz];
			if (inv2Image[xyz]<min) min = inv2Image[xyz];
		}
		// re-normalize
		for (int xyz=0;xyz<nxyz;xyz++) {
			inv2Image[xyz] = (inv2Image[xyz]-min)/(max-min);
		}
		System.out.println("image range: "+min+", "+max);

		double mean = 0.0;
		double den = 0.0;
		for (int xyz=0;xyz<nxyz;xyz++) if (!skip0Param || inv2Image[xyz]!=0) {
			mean += inv2Image[xyz];
			den++;
		}
		mean /= den;
		System.out.println("mean parameters: "+mean);
		
		double pi = Math.PI;
		
		byte model = EXP;
		//if (bgParam.getValue().equals("half-normal")) model = HNORM;
		
		// re-normalized probability map
		float[] proba = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (model==HNORM)
				// half-normal model
				proba[xyz] = (float)(2.0/(mean*pi)*Math.exp(-inv2Image[xyz]*inv2Image[xyz]/(pi*mean*mean)));
			else if (model==EXP)
				// exponential model
				proba[xyz] = (float)(Math.exp(-inv2Image[xyz]/mean)/mean);
			// bg vs. outlier test : proba is p(xyz is background)
			proba[xyz] = proba[xyz]/(1.0f+proba[xyz]);
		}
		// loop
		double diff = 1.0;
		for (int t=0;t<itermax && diff>maxdiff;t++) {
			System.out.println("iteration "+(t+1));
			diff = mean;
			mean = 0.0;
			den = 0.0;
			for (int xyz=0;xyz<nxyz;xyz++) if (!skip0Param || inv2Image[xyz]!=0) {
				mean += proba[xyz]*inv2Image[xyz];
				den += proba[xyz];
			}
			mean /= den;
			System.out.println("mean parameters: "+mean);
			diff = Numerics.abs(diff-mean);
			System.out.println("diff parameters: "+diff);
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (model==HNORM)
					// half-normal model
					proba[xyz] = (float)(2.0/(mean*pi)*Math.exp(-inv2Image[xyz]*inv2Image[xyz]/(pi*mean*mean)));
				else if (model==EXP)
					// exponential model
					proba[xyz] = (float)(Math.exp(-inv2Image[xyz]/mean)/mean);
				// bg vs. outlier test : proba is p(xyz is background)
				proba[xyz] = proba[xyz]/(1.0f+proba[xyz]);
			}
		}

		BasicInfo.displayMessage("background-based skull stripping");
		
		// start from the bg mask
		MinMaxFiltering minmax = new MinMaxFiltering(proba, nx,ny,nz, rx,ry,rz);
		
		float[] brain = minmax.growRegion(new float[]{0.0f}, new float[]{0.9f}, new float[]{0.9f}, 16, 10, false);
		
		BinaryTopology topo = null;
		Gdm3d gdm = null;

		// final refinement using T1, iso images + GDM
		if (t1mapImage!=null || isoImage!=null) {
			float t1max = 0.0f, isomax = 0.0f, inv2max = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (t1mapImage!=null && t1mapImage[xyz]>t1max) t1max = t1mapImage[xyz];
				if (isoImage!=null && isoImage[xyz]>isomax) isomax = isoImage[xyz];
				if (inv2Image[xyz]>inv2max) inv2max = inv2Image[xyz];
			}
			// keep CSF inside + regions away from boundary (=> limit number of iterations)
			int[] mask = new int[nxyz];
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (brain[xyz]>0) {
					mask[xyz] = 1;
				} else {
					mask[xyz] = 0;
				}
			}
			
			float[] balloon = new float[nxyz];
			for (int xyz=0;xyz<nxyz;xyz++) {
				float force = -0.1f;
				if (t1mapImage!=null) force = Numerics.max(force, Numerics.bounded( (t1mapImage[xyz]/t1max-0.75f)/0.25f, -1.0f, 1.0f));
				if (isoImage!=null) force = Numerics.max(force, 1.0f-Numerics.bounded( (isoImage[xyz]/isomax-0.25f)/0.25f, -1.0f, 1.0f));
				balloon[xyz] = force;
			}
			
			// topology correction for the mask?
			topo = new BinaryTopology(mask, nx, ny, nz, rx, ry, rz, "wcs", lutdir);
			
			topo.outsideSphericalTopology();
			
			int[] toposeg = topo.exportIntSegmentation();
			gdm = new Gdm3d(toposeg, nx, ny, nz, rx, ry, rz, null, balloon, 0.0f, 0.1f, 0.9f, "wcs", lutdir);
					
			gdm.evolveNarrowBand(100, 0.001f);

			brain = gdm.exportSegmentation();
		}
		
		// generate outputs
		byte[] brainmask = new byte[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (brain[xyz]>0) brainmask[xyz] = 1;
			else brainmask[xyz] = 0;
		}
		brainmaskImage = brainmask;

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (brainmaskImage[xyz]==0) inv2Image[xyz] = 0.0f;
		}
		inv2maskImage = inv2Image;
	
		if (t1mapImage!=null) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (brainmaskImage[xyz]==0) t1mapImage[xyz] = 0.0f;
			}
			t1maskImage = t1mapImage;
		}
		
		if (isoImage!=null) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (brainmaskImage[xyz]==0) isoImage[xyz] = 0.0f;
			}
			isomaskImage = isoImage;
		}
	}
}
