package de.mpg.cbs.core.brain;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class BrainMp2rageDuraEstimation {

	// parameters
	public		static final String[]	outTypes = {"dura_region","boundary","dura_prior","bg_prior", "intens_prior"};
	
	private float[] inv2Image;
	private int[] maskImage;
	private float[] resultImage;

	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private String outParam = "dura_region";
	private float distParam = 5.0f;
	
	// create inputs
	public final void setSecondInversionImage(float[] val) { inv2Image = val; }
	public final void setSkullStrippingMask(int[] val) { maskImage = val; }

	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

    public final void setDistanceToBackground_mm(float val) { distParam = val; }
	public final void setOutputType(String val) { outParam = val; }
	
	
	// to be used for JIST definitions, generic info / help
	public static final String getPackage() { return "CBS Tools"; }
	public static final String getCategory() { return "Brain Processing"; }
	public static final String getLabel() { return "MP2RAGE Dura Estimatio"; }
	public static final String getName() { return "Mp2rageDuraEstimation"; }

	public final String[] getAlgorithmAuthors() {return new String[] {"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Filters a MP2RAGE brain image to obtain a probability map of dura matter."; }
		
	public static final String getVersion() { return "3.1"; }

	// create outputs
	public final float[] getDuraImage() { return resultImage; }

	public final void execute() {
		
		float[][][] inv2img = new float[nx][ny][nz];
		byte[][][] maskimg = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    inv2img[x][y][z] = inv2Image[x+nx*y+nx*ny*z];
		    maskimg[x][y][z] = (byte)maskImage[x+nx*y+nx*ny*z];
		}
		inv2Image = null;
		maskImage = null;
		    
		// main algorithm
		
		// 1. Estimate regions of potential partial volume
		float[][][] pvmap = new float[nx][ny][nz];
		float pvavg = 0.0f;
		float pvsig = 0.0f;
		float pvden = 0.0f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			if (!zeroNeighbor(inv2img,x,y,z,2)) {
				// check for zero-valued neighbors as well
				pvmap[x][y][z] = 0.0f;
				float best = 0.0f;
				float dmax = 13;
				for (int d=0;d<dmax;d++) {
					float val = planeScore(inv2img, x,y,z,d);
					if (val*val>best*best) best = val;
				}
				pvmap[x][y][z] = best;
				pvavg += best;
				pvden++;
			}
		}
		pvavg /= pvden;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			//if (!zeroNeighbor(inv2img,x,y,z,2)) {
			pvsig += Numerics.square(pvmap[x][y][z]-pvavg);
		}
		pvsig /= (pvden-1.0f);
		// probability
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			if (pvmap[x][y][z]>pvavg) {
				pvmap[x][y][z] = 1.0f - (float)Math.exp( -0.5f*(pvmap[x][y][z]-pvavg)*(pvmap[x][y][z]-pvavg)/pvsig );
			} else {
				pvmap[x][y][z] = 0.0f;
			}
		}

		// 2. Expand background region
		float scale = distParam/rx;
		float[][][] bgmap = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			bgmap[x][y][z] = maskimg[x][y][z];
		}
		float[][] kernel = ImageFilters.separableGaussianKernel(scale, scale, scale);
		int kx = (kernel[0].length-1)/2;
		bgmap = ImageFilters.separableConvolution(bgmap, nx, ny, nz, kernel, kx, kx, kx);
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			bgmap[x][y][z] = 2.0f*(1.0f-bgmap[x][y][z])*maskimg[x][y][z];
		}

		// 3. Find regions of low intensity
		// use default parameters: 
		float maxdiff = 0.0001f;
		int itermax = 20;
		
		// min max for t1 map fixing
		float max = -1e10f;
		float min = 1e10f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (inv2img[x][y][z]>max) max = inv2img[x][y][z];
			if (inv2img[x][y][z]<min) min = inv2img[x][y][z];
		}
		// re-normalize
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			inv2img[x][y][z] = (inv2img[x][y][z]-min)/(max-min);
		}
		System.out.println("image range: "+min+", "+max);

		double mean = 0.0f;
		double den = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (inv2img[x][y][z]!=0) {
				mean += inv2img[x][y][z];
				den++;
			}
		}
		mean /= den;
		System.out.println("mean parameters: "+mean);
		
		// re-normalized probability map
		float[][][] intensmap = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			intensmap[x][y][z] = (float)(FastMath.exp(-inv2img[x][y][z]/mean)/mean);
			intensmap[x][y][z] = intensmap[x][y][z]/(1.0f+intensmap[x][y][z]);
		}
		// loop
		double diff = 1.0f;
		for (int t=0;t<itermax && diff>maxdiff;t++) {
			System.out.println("iteration "+(t+1));
			diff = mean;
			mean = 0.0f;
			den = 0.0f;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (inv2img[x][y][z]!=0) {
					mean += intensmap[x][y][z]*inv2img[x][y][z];
					den += intensmap[x][y][z];
				}
			}
			mean /= den;
			System.out.println("mean parameters: "+mean);
			diff = Numerics.abs(diff-mean);
			System.out.println("diff parameters: "+diff);
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				intensmap[x][y][z] = (float)(FastMath.exp(-(inv2img[x][y][z])/mean)/mean);
				intensmap[x][y][z] = intensmap[x][y][z]/(1.0f+intensmap[x][y][z]);
			}
		}

		// 4. Combine information sources, expand PV in regions of higher BG
		float[] result = new float[nxyz];
		if (outParam.equals("dura_prior")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				result[x+nx*y+nx*ny*z] = pvmap[x][y][z];
			}
		} else 
		if (outParam.equals("bg_prior")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				result[x+nx*y+nx*ny*z] = Numerics.min(1.0f, bgmap[x][y][z]);
			}
		} else 
		if (outParam.equals("intens_prior")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				result[x+nx*y+nx*ny*z] = intensmap[x][y][z];
			}
		} else 
		if (outParam.equals("boundary")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				//result[x][y][z] = (float)Math.sqrt(bgmap[x][y][z]*Numerics.max(pvmap[x][y][z],intensmap[x][y][z]));
				result[x+nx*y+nx*ny*z] = Numerics.min(1.0f, bgmap[x][y][z])*Numerics.max(pvmap[x][y][z],intensmap[x][y][z]);
			}
		} else 
		if (outParam.equals("dura_region")) {
			// take the product
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				//result[x][y][z] = (float)Math.sqrt(bgmap[x][y][z]*Numerics.max(pvmap[x][y][z],intensmap[x][y][z]));
				result[x+nx*y+nx*ny*z] = Numerics.min(1.0f, bgmap[x][y][z])*Numerics.max(pvmap[x][y][z],intensmap[x][y][z]);
			}
			// propagate maximum values toward background
			float[] tmp = new float[nxyz];
			for (int t=0;t<2*scale;t++) {
				for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				    int xyz = x+nx*y+nx*ny*z;
					tmp[xyz] = result[xyz];
					if (bgmap[x+1][y][z]<bgmap[x][y][z] && result[xyz+1]>tmp[xyz]) tmp[xyz] = result[xyz+1];
					if (bgmap[x-1][y][z]<bgmap[x][y][z] && result[xyz-1]>tmp[xyz]) tmp[xyz] = result[xyz-1];
					if (bgmap[x][y+1][z]<bgmap[x][y][z] && result[xyz+nx]>tmp[xyz]) tmp[xyz] = result[xyz+nx];
					if (bgmap[x][y-1][z]<bgmap[x][y][z] && result[xyz-nx]>tmp[xyz]) tmp[xyz] = result[xyz-nx];
					if (bgmap[x][y][z+1]<bgmap[x][y][z] && result[xyz+nx*ny]>tmp[xyz]) tmp[xyz] = result[xyz+nx*ny];
					if (bgmap[x][y][z-1]<bgmap[x][y][z] && result[xyz-nx*ny]>tmp[xyz]) tmp[xyz] = result[xyz-nx*ny];
				}
				for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
					result[x+nx*y+nx*ny*z] = tmp[x+nx*y+nx*ny*z];
				}
			}
			/*
			// square the result to dampen lower values?
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				result[x][y][z] = result[x][y][z]*result[x][y][z];
			}
			*/
		}
		
		resultImage = result;		
		
		return;
	}

	boolean zeroNeighbor(float[][][] image, int x, int y, int z, int d) {
		for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int l=-d;l<=d;l++) {
			if (image[x+i][y+j][z+l]!=image[x][y][z] && i*i+j*j+l*l<=2*d*d) return false;
		}
		return true;
	}
	
	float planeScore(float[][][] image, int x, int y, int z, int d) {
		float val = 0.0f;
		if (d==0) {		
			val =(2.0f*image[x][y][z]-image[x-1][y][z]-image[x+1][y][z]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z]-image[x+1][y-1][z]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z]-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y][z-1]-image[x+1][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y][z+1]-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z-1]-image[x-1][y-1][z-1]-image[x+1][y-1][z-1]
				 +2.0f*image[x][y-1][z+1]-image[x-1][y-1][z+1]-image[x+1][y-1][z+1]
				 +2.0f*image[x][y+1][z-1]-image[x-1][y+1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x][y+1][z+1]-image[x-1][y+1][z+1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==1) {
			val =(2.0f*image[x][y][z]-image[x][y-1][z]-image[x][y+1][z]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z]-image[x-1][y+1][z]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z]-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]-image[x][y-1][z-1]-image[x][y+1][z-1]
				 +2.0f*image[x][y][z+1]-image[x][y-1][z+1]-image[x][y+1][z+1]
				 +2.0f*image[x-1][y][z-1]-image[x-1][y-1][z-1]-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y][z+1]-image[x-1][y-1][z+1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y][z-1]-image[x+1][y-1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y][z+1]-image[x+1][y-1][z+1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==2) { 			
			val =(2.0f*image[x][y][z]-image[x][y][z-1]-image[x][y][z+1]
				 +2.0f*image[x-1][y][z]-image[x-1][y][z-1]-image[x-1][y][z+1]
				 +2.0f*image[x+1][y][z]-image[x+1][y][z-1]-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z]-image[x][y-1][z-1]-image[x][y-1][z+1]
				 +2.0f*image[x][y+1][z]-image[x][y+1][z-1]-image[x][y+1][z+1]
				 +2.0f*image[x-1][y-1][z]-image[x-1][y-1][z-1]-image[x-1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z]-image[x-1][y+1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z]-image[x+1][y-1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y+1][z]-image[x+1][y+1][z-1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==3) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z]-image[x+1][y+1][z]
				 +2.0f*image[x-1][y+1][z]-image[x-2][y][z]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z]-image[x][y-2][z]-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y-1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-1]-image[x][y+2][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y-1][z+1]-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+1]-image[x][y+2][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/18.0f;
		} else if (d==4) { 			
			val =(2.0f*image[x][y][z]-image[x][y-1][z-1]-image[x][y+1][z+1]
				 +2.0f*image[x][y+1][z-1]-image[x][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x][y-1][z+1]-image[x][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-1][y][z-2]-image[x-1][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-1][y-2][z]-image[x-1][y][z+2]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x+1][y][z-2]-image[x+1][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x+1][y-2][z]-image[x+1][y][z+2])/18.0f;
		} else if (d==5) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y][z-1]-image[x+1][y][z+1]
				 +2.0f*image[x+1][y][z-1]-image[x][y][z-2]-image[x+2][y][z]
				 +2.0f*image[x-1][y][z+1]-image[x-2][y][z]-image[x][y][z+2]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-1][z-2]-image[x+2][y-1][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-1][z]-image[x][y-1][z+2]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+1][z-2]-image[x+2][y+1][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+1][z]-image[x][y+1][z+2])/18.0f;
		} else if (d==6) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y+1][z]-image[x+1][y-1][z]
				 +2.0f*image[x-1][y-1][z]-image[x-2][y][z]-image[x][y-2][z]
				 +2.0f*image[x+1][y+1][z]-image[x][y-2][z]-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y+1][z-1]-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-1]-image[x][y-2][z-1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y+1][z+1]-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z+1]-image[x][y-2][z+1]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/18.0f;
		} else if (d==7) { 			
			val =(2.0f*image[x][y][z]-image[x][y-1][z+1]-image[x][y+1][z-1]
				 +2.0f*image[x][y-1][z-1]-image[x][y-2][z]-image[x][y][z-2]
				 +2.0f*image[x][y+1][z+1]-image[x][y][z+2]-image[x][y+2][z]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z+1]-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-1][y-2][z]-image[x-1][y][z-2]
				 +2.0f*image[x-1][y+1][z+1]-image[x-1][y][z+2]-image[x-1][y+2][z]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x+1][y-2][z]-image[x+1][y][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+1][y][z+2]-image[x+1][y+2][z])/18.0f;
		} else if (d==8) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y][z+1]-image[x+1][y][z-1]
				 +2.0f*image[x-1][y][z-1]-image[x-2][y][z]-image[x][y][z-2]
				 +2.0f*image[x+1][y][z+1]-image[x][y][z+2]-image[x+2][y][z]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z+1]-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-1][z]-image[x][y-1][z-2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-1][z+2]-image[x+2][y-1][z]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y+1][z]-image[x][y+1][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+1][z+2]-image[x+2][y+1][z])/18.0f;
		} else if (d==9) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-2]-image[x+2][y][z]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y][z-2]-image[x+2][y+2][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z]-image[x][y+2][z+2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z]-image[x+2][y][z+2])/14.0f;
		} else if (d==10) { 			
			val =(2.0f*image[x][y][z]-image[x+1][y-1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x+2][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x+2][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x][y-2][z-2]-image[x-2][y][z]
				 +2.0f*image[x-1][y+1][z-1]-image[x][y][z-2]-image[x-2][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x][y-2][z]-image[x-2][y][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+2][y][z]-image[x][y+2][z+2])/14.0f;
		} else if (d==11) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y+1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+2][z]-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+2][z-2]-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-2]-image[x][y-2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y][z-2]-image[x+2][y-2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z]-image[x][y-2][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+2][z]-image[x+2][y][z+2])/14.0f;
		} else if (d==12) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+2]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+2]-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-2][z]-image[x][y][z-2]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z]-image[x+2][y][z-2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z]-image[x][y+2][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y][z+2]-image[x+2][y+2][z])/14.0f;
		}
		return 0.5f*val;
	}

}
