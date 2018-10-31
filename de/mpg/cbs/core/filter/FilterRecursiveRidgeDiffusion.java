package de.mpg.cbs.core.filter;

import de.mpg.cbs.utilities.Ngb;
import de.mpg.cbs.utilities.Numerics;
import de.mpg.cbs.utilities.BasicInfo;
import de.mpg.cbs.libraries.ImageFilters;
import de.mpg.cbs.libraries.ImageStatistics;
import de.mpg.cbs.libraries.ObjectExtraction;
import de.mpg.cbs.libraries.ObjectLabeling;
import de.mpg.cbs.structures.BinaryHeap2D;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.special.Erf;
//import org.apache.commons.math3.analysis.*;

import java.util.BitSet;

public class FilterRecursiveRidgeDiffusion {
	
	// jist containers
	private float[] inputImage;
	private String brightParam = "bright";
	public static final String[] brightTypes = {"bright","dark","both"};
	
	private float[] surfaceImage;
	private String	orientationParam = "undefined";
	public static final String[] orientationTypes = {"parallel","orthogonal","undefined"};
	private float	angleParam = 1.0f;
	
	private float[] locationImage;
	
	private String filterParam = "2D";
	public static final String[] filterTypes = {"2D","1D","0D"};
	private int minscaleParam = 0;
	private int maxscaleParam = 3;
	
	private String propagationParam = "diffusion";
	public static final String[] propagationTypes = {"diffusion","belief","SUR","flooding","growing","labeling","none"};
	private float difffactorParam = 1.0f;
	private float simscaleParam = 0.1f;
	private int ngbParam = 4;
	private int iterParam = 100;
	private float maxdiffParam = 0.001f;
	
	private float[] pvImage;
	private float[] filterImage;
	private float[] probaImage;
	private float[] propagImage;
	private int[] scaleImage;
	private float[] directionImage;
	private float[] correctImage;
	private float[] sizeImage;
	
	// global variables
	private int nx, ny, nz, nc, nxyz;
	private float rx, ry, rz;
	
	// numerical quantities
	private static final	float	INVSQRT2 = (float)(1.0/FastMath.sqrt(2.0));
	private static final	float	INVSQRT3 = (float)(1.0/FastMath.sqrt(3.0));
	private static final	float	SQRT2 = (float)FastMath.sqrt(2.0);
	private static final	float	SQRT3 = (float)FastMath.sqrt(3.0);
	private static final	float	SQRT2PI = (float)FastMath.sqrt(2.0*(float)Math.PI);
	private static final	float	PI2 = (float)(Math.PI/2.0);
	private static final	float   	L2N2=2.0f*(float)(FastMath.sqrt(2.0*(float)(FastMath.log(2.0))));
	
	// direction labeling		
	public	static	final	byte	X = 0;
	public	static	final	byte	Y = 1;
	public	static	final	byte	Z = 2;
	public	static 	final 	byte 	XpY = 3;
	public	static 	final 	byte 	YpZ = 4;
	public	static 	final 	byte 	ZpX = 5;
	public	static 	final 	byte 	XmY = 6;
	public	static 	final 	byte 	YmZ = 7;
	public	static 	final 	byte 	ZmX = 8;
	public	static 	final 	byte 	XpYpZ = 9;
	public	static 	final 	byte 	XmYmZ = 10;
	public	static 	final 	byte 	XmYpZ = 11;
	public	static 	final 	byte 	XpYmZ = 12;
	public	static	final	byte	mX = 13;
	public	static	final	byte	mY = 14;
	public	static	final	byte	mZ = 15;
	public	static 	final 	byte 	mXpY = 16;
	public	static 	final 	byte 	mYpZ = 17;
	public	static 	final 	byte 	mZpX = 18;
	public	static 	final 	byte 	mXmY = 19;
	public	static 	final 	byte 	mYmZ = 20;
	public	static 	final 	byte 	mZmX = 21;
	public	static 	final 	byte 	mXpYpZ = 22;
	public	static 	final 	byte 	mXmYmZ = 23;
	public	static 	final 	byte 	mXmYpZ = 24;
	public	static 	final 	byte 	mXpYmZ = 25;
	
	public	static 	final 	byte 	NC = 13;
	public	static 	final 	byte 	NC2 = 26;
	
	private static final boolean debug=true;
	
	//set inputs
	public final void setInputImage(float[] val) { inputImage = val; }
	public final void setRidgeIntensities(String val) { brightParam = val; }
	public final void setRidgeFilter(String val) { filterParam = val; }
	public final void setSurfaceLevelSet(float[] val) { surfaceImage = val; }
	public final void setOrientationToSurface(String val) { orientationParam = val; }
	public final void setAngularFactor(float val) { angleParam = val; }
		
	public final void setLocationPrior(float[] val) { locationImage = val; }

	public final void setMinimumScale(int val) { minscaleParam = val; }
	public final void setMaximumScale(int val) { maxscaleParam = val; }

	public final void setPropagationModel(String val) { propagationParam = val; }
	public final void setDiffusionFactor(float val) { difffactorParam = val; }
	public final void setSimilarityScale(float val) { simscaleParam = val; }
	public final void setNeighborhoodSize(int val) { ngbParam = val; }
	public final void setMaxIterations(int val) { iterParam = val; }
	public final void setMaxDifference(float val) { maxdiffParam = val; }
		
	// set generic inputs	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	//set JIST definitions
	//to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Filter.devel"; }
	public final String getLabel() { return "Recursive Ridge Diffusion"; }
	public final String getName() { return "RecursiveRidgeDiffusion"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Extracts planar of tubular structures across multiple scales, with an optional directional bias."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.2"; };
			
	// create outputs
	public final float[] getRidgePartialVolumeImage() { return pvImage;}
	public final float[] getFilterResponseImage() { return filterImage;}
	public final float[] getProbabilityResponseImage() { return probaImage;}
	public final float[] getPropagatedResponseImage() { return propagImage;}
	public final int[] getDetectionScaleImage() { return scaleImage;}
	public final float[] getRidgeDirectionImage() { return directionImage;}
	public final float[] getDirectionalCorrectionImage() { return correctImage;}
	public final float[] getRidgeSizeImage() { return sizeImage;}

	public void execute(){
		BasicInfo.displayMessage("recursive ridge diffusion:\n");
		
		// import the inputImage data into 1D arrays: already done
		BasicInfo.displayMessage("...load data\n");
		
		// find surfaceImage dimension, if multiple surfaces are used
		nc = 0;
		if (surfaceImage==null) orientationParam = "undefined";
		
		if (orientationParam.equals("parallel") || orientationParam.equals("orthogonal")) {
			BasicInfo.displayMessage("...load surfaceImage orientationParam(s)\n");
			nc = Numerics.round(surfaceImage.length/nxyz);
		}
		
		
		boolean[] mask = new boolean[nxyz];
		float minI = 1e9f;
		float maxI = -1e9f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			// mask
			if (inputImage[id]==0) mask[id] = false;
			else mask[id] = true;
			// remove border from computations
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[id] = false;
			// normalize
			if (inputImage[id]>maxI) maxI = inputImage[id];
			if (inputImage[id]<minI) minI = inputImage[id];
		}
		
		// normalize, invert inputImage if looking for dark features
		BasicInfo.displayMessage("...normalize intensities (detection: "+brightParam+")\n");
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (brightParam.equals("bright"))
				inputImage[xyz] = (inputImage[xyz]-minI)/(maxI-minI);
			else if (brightParam.equals("dark"))
				inputImage[xyz] = (maxI-inputImage[xyz])/(maxI-minI);
			else 
				inputImage[xyz] = (inputImage[xyz]-minI)/(maxI-minI);
		}
		boolean unidirectional = true;
		if (brightParam.equals("both")) unidirectional = false;
		
		// Compute filter at different scales
		// new filter response from raw inputImage		
		float[] maxresponse = new float[nxyz];
		byte[] maxdirection = new byte[nxyz];
		int[] maxscale = new int[nxyz];
		
		if (minscaleParam==0) {
            BasicInfo.displayMessage("...first filter response\n");
            
            if (filterParam.equals("0D")) directionFromRecursiveRidgeFilter0D(inputImage, mask, maxresponse, maxdirection, unidirectional);
            else if (filterParam.equals("1D")) directionFromRecursiveRidgeFilter1D(inputImage, mask, maxresponse, maxdirection, unidirectional);
            else if (filterParam.equals("2D")) directionFromRecursiveRidgeFilter2D(inputImage, mask, maxresponse, maxdirection, unidirectional);
        }
        
		for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
			if (maxresponse[xyz]>0) {
				maxscale[xyz] = 0;
			} else {
				maxscale[xyz] = -1;
			}
		}
		for (int i=Numerics.max(minscaleParam-1,0);i<maxscaleParam;i++) {
			float scale = 1.0f+i;
			float[] smoothed = new float[nxyz];

			BasicInfo.displayMessage("...filter response at scale "+scale+"\n");
		
			// Gaussian Kernel
			float[][] G = ImageFilters.separableGaussianKernel(scale/L2N2,scale/L2N2,scale/L2N2);
			int gx = (G[0].length-1)/2;
			int gy = (G[1].length-1)/2;
			int gz = (G[2].length-1)/2;
				
			// smoothed inputImage
			smoothed = ImageFilters.separableConvolution(inputImage,nx,ny,nz,G,gx,gy,gz); 

			byte[] direction = new byte[nxyz];
			float[] response = new float[nxyz];
			if (filterParam.equals("0D")) directionFromRecursiveRidgeFilter0D(smoothed, mask, response, direction, unidirectional);
			else if (filterParam.equals("1D")) directionFromRecursiveRidgeFilter1D(smoothed, mask, response, direction, unidirectional);
			else if (filterParam.equals("2D")) directionFromRecursiveRidgeFilter2D(smoothed, mask, response, direction, unidirectional);
			
			//Combine scales: keep maximum response
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				if (response[xyz]>maxresponse[xyz]) {
					maxresponse[xyz] = response[xyz];
					maxdirection[xyz] = direction[xyz];
					maxscale[xyz] = (1+i);
				}
			}
		}
		
		// Equalize histogram (Exp Median)
		BasicInfo.displayMessage("...normalization into probabilities\n");
		float[] proba = new float[nxyz];
		probabilityFromRecursiveRidgeFilter(maxresponse, proba);	
		
		// Attenuate response according to orientationParam
		float[] correction = new float[nxyz];
		if (orientationParam.equals("parallel") || orientationParam.equals("orthogonal")) {
			BasicInfo.displayMessage("...orientationParam response filtering\n");
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				float[] vec = directionVector(maxdirection[xyz]);
				// find the closest surfaceImages; weighted mean
				int closest = 0;
				float dist = Numerics.abs(surfaceImage[xyz]);
				for (int c=1;c<nc;c++) {
					if (Numerics.abs(surfaceImage[xyz+c*nxyz])<dist) {
						closest = c;
						dist = Numerics.abs(surfaceImage[xyz+c*nxyz]);
					}
				}
				int second = closest;
				if (closest<nc-1 && surfaceImage[xyz+closest*nxyz]*surfaceImage[xyz+(closest+1)*nxyz]<0) second = closest+1; 
				if (closest>0 && surfaceImage[xyz+closest*nxyz]*surfaceImage[xyz+(closest-1)*nxyz]>0) second = closest-1; 
				
				float w1 = 1.0f, w2 = 0.0f;
				if (second!=closest) {
					w1 = surfaceImage[xyz+second*nxyz]/(surfaceImage[xyz+second*nxyz]-surfaceImage[xyz+closest*nxyz]);
					w2 = surfaceImage[xyz+closest*nxyz]/(surfaceImage[xyz+closest*nxyz]-surfaceImage[xyz+second*nxyz]);
				}
				double[] grad = new double[3];
				// central differences
				grad[X] = w1*0.5f*(surfaceImage[xyz+closest*nxyz+1]-surfaceImage[xyz+closest*nxyz-1]);
				grad[Y] = w1*0.5f*(surfaceImage[xyz+closest*nxyz+nx]-surfaceImage[xyz+closest*nxyz-nx]);
				grad[Z] = w1*0.5f*(surfaceImage[xyz+closest*nxyz+nx*ny]-surfaceImage[xyz+closest*nxyz-nx*ny]);
				if (second!=closest) {
					grad[X] += w2*0.5f*(surfaceImage[xyz+second*nxyz+1]-surfaceImage[xyz+second*nxyz-1]);
					grad[Y] += w2*0.5f*(surfaceImage[xyz+second*nxyz+nx]-surfaceImage[xyz+second*nxyz-nx]);
					grad[Z] += w2*0.5f*(surfaceImage[xyz+second*nxyz+nx*ny]-surfaceImage[xyz+second*nxyz-nx*ny]);
				}					
				double norm = grad[X]*grad[X]+grad[Y]*grad[Y]+grad[Z]*grad[Z];
				if (norm>0) {
					norm = FastMath.sqrt(norm);
					grad[X] /= norm;
					grad[Y] /= norm;
					grad[Z] /= norm;
				}
				// angle projected back into [0,1], 0: parallel, 1:orthogonal
				double angle = 2.0*FastMath.acos(grad[X]*vec[X]+grad[Y]*vec[Y]+grad[Z]*vec[Z])/FastMath.PI;
				angle = Numerics.min(angle,2.0-angle);
				if (orientationParam.equals("parallel") && filterParam.equals("1D")) {
					correction[xyz] = (float)FastMath.pow(angle, angleParam);
				} else if (orientationParam.equals("parallel") && filterParam.equals("2D")) {
					correction[xyz] = (float)FastMath.pow(1.0-angle, angleParam);
				} else if (orientationParam.equals("orthogonal") && filterParam.equals("1D")) {
					correction[xyz] = (float)FastMath.pow(1.0-angle, angleParam);
				} else if (orientationParam.equals("orthogonal") && filterParam.equals("2D")) {
					correction[xyz] = (float)FastMath.pow(angle, angleParam);
				} else {
					correction[xyz] = 1.0f;
				}
				//maxresponse[xyz] *= correction[xyz];
				proba[xyz] *= correction[xyz];
			}
		}
		
		// mask response with locationImage prior before equalization
		if (locationImage!=null) {
			//for (int xyz=0;xyz<nxyz;xyz++) maxresponse[xyz] *= locationImage[xyz];
			for (int xyz=0;xyz<nxyz;xyz++) proba[xyz] *= locationImage[xyz];
		}

		/* to be done before the other prior selection, in order to maintain the noise / response structure? 
		// Equalize histogram (Exp Median)
		BasicInfo.displayMessage("...normalization into probabilities\n");
		float[] proba = new float[nxyz];
		probabilityFromRecursiveRidgeFilter(maxresponse, proba);	
		
		float[] proba = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			proba[xyz] = maxresponse[xyz];
		}
		*/
		
		
		// rescale the probability response
		float pmax = ImageStatistics.robustMaximum(proba, 0.000001f, 6, nx, ny, nz);
		if (pmax>0) for (int xyz=0;xyz<nxyz;xyz++) proba[xyz] = Numerics.min(proba[xyz]/pmax,1.0f);
		
		// generate a direction vector
		float[][] direction = new float[3][nxyz];
		for (int xyz=0;xyz<nxyz;xyz++){
			float[] vec = directionVector(maxdirection[xyz]);
			direction[X][xyz] = vec[X];
			direction[Y][xyz] = vec[Y];
			direction[Z][xyz] = vec[Z];
		}
		
		// 3. diffuse the data to neighboring structures with SUR
		BasicInfo.displayMessage("...diffusion ("+propagationParam+")\n");
		
		float[] propag = new float[nxyz];
		if (filterParam.equals("0D")) {
		     propag = spatialExpansion0D(proba, maxscale, inputImage);
		} else
		if (propagationParam.equals("SUR")) {
			if (filterParam.equals("1D"))  propag = spatialUncertaintyRelaxation1D(proba, maxdirection, ngbParam, simscaleParam, difffactorParam, iterParam);
			else if (filterParam.equals("2D"))  propag = spatialUncertaintyRelaxation2D(proba, maxdirection, ngbParam, simscaleParam, difffactorParam, iterParam);
		} else if  (propagationParam.equals("diffusion")) {
			if (filterParam.equals("1D"))  propag = probabilisticDiffusion1D(proba, maxdirection, ngbParam, maxdiffParam, simscaleParam, difffactorParam, iterParam);
			else if (filterParam.equals("2D"))  propag = probabilisticDiffusion2D(proba, maxdirection, ngbParam, maxdiffParam, simscaleParam, difffactorParam, iterParam);
		} else if  (propagationParam.equals("belief")) {
			if (filterParam.equals("1D"))  propag = beliefPropagation1D(proba, maxdirection, ngbParam, simscaleParam, iterParam, maxdiffParam);
			else if (filterParam.equals("2D"))  propag = beliefPropagation2D(proba, maxdirection, ngbParam, simscaleParam, iterParam, maxdiffParam);
		} else if  (propagationParam.equals("flooding")) {
			if (filterParam.equals("1D"))  propag = probabilisticDiffusion1D(proba, maxdirection, ngbParam, maxdiffParam, simscaleParam, difffactorParam, iterParam);
			else if (filterParam.equals("2D"))  propag = probabilisticFlooding2D(proba, maxdirection, ngbParam, maxdiffParam, simscaleParam, difffactorParam, iterParam);
		} else if  (propagationParam.equals("growing")) {
			if (filterParam.equals("1D"))  propag = regionGrowing1D(proba, maxdirection, ngbParam, maxdiffParam, simscaleParam, difffactorParam, iterParam);
			else if (filterParam.equals("2D"))  propag = regionGrowing2D(proba, maxdirection, ngbParam, maxdiffParam, simscaleParam, difffactorParam, iterParam);
		} else if  (propagationParam.equals("labeling")) {
			if (filterParam.equals("1D"))  propag = regionLabeling1D(proba, maxdirection, ngbParam, maxdiffParam, simscaleParam, difffactorParam, iterParam);
			else if (filterParam.equals("2D"))  propag = regionLabeling2D(proba, maxdirection, ngbParam, maxdiffParam, simscaleParam, difffactorParam, iterParam);
		} else {
		    propag = proba;
		}
		
		/*
		probabilityFromRecursiveRidgeFilter(propag, proba);	
		
		// rescale the probability response
		float pmax = ImageStatistics.robustMaximum(proba, 0.0001f, 5, nx, ny, nz);
		if (pmax>0) for (int xyz=0;xyz<nxyz;xyz++) proba[xyz] = Numerics.min(proba[xyz]/pmax,1.0f);
		*/
		
		
		// 4. Measure size of connected region
		boolean[] detected = ObjectExtraction.objectFromImage(propag, nx, ny, nz, 0.5f, ObjectExtraction.SUPEQUAL);
		int[] components = ObjectLabeling.connected26Object3D(detected, nx, ny, nz);
		int Ncomp = ObjectLabeling.countLabels(components, nx, ny, nz);
		float[] length = new float[Ncomp];
		for (int xyz=0;xyz<nxyz;xyz++) if (components[xyz]>0) {
			length[components[xyz]-1] += propag[xyz];
		}
		float[] lengthmap = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) if (components[xyz]>0) {
			lengthmap[xyz] = length[components[xyz]-1];
		}
		
		// 5. Get the partial volume estimates : note they may be incorrect after propagation
		float[] pvol = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (propag[xyz]>0.0f) {
				// when the diameter is smaller than the voxel, linear approx
				pvol[xyz] = Numerics.bounded(propag[xyz], 0.0f, 1.0f);
				float radius = Numerics.max(0.0f, maxscale[xyz]/2.0f); // scale is diameter
				int nsc = Numerics.ceil(radius); 
				for (int i=-nsc;i<=nsc;i++) for (int j=-nsc;j<=nsc;j++) for (int k=-nsc;k<=nsc;k++) {
					if (x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+k>=0 && z+k<nz) {
						float dist = (float)FastMath.sqrt(i*i+j*j+k*k);
						// update the neighbor value for diameters above the voxel size
						float ratio = Numerics.bounded(radius+0.5f-dist, 0.0f, 1.0f);
						//pvol[xyz+i+j*nx+k*nx*ny] = Numerics.max(pvol[xyz+i+j*nx+k*nx*ny], Numerics.bounded(scaleMax[xyz]/2.0f+0.5f-dist, 0.0f, 1.0f));
						pvol[xyz+i+j*nx+k*nx*ny] = Numerics.max(pvol[xyz+i+j*nx+k*nx*ny], ratio*pvol[xyz]);
					}
				}
			}
		}
		
		// Output
		BasicInfo.displayMessage("...output inputImages\n");
		filterImage = maxresponse;
		probaImage = proba;
		propagImage = propag;
		scaleImage = maxscale;
		pvImage = pvol;
		
		float[] result4 = new float[nxyz*3];
		for (int xyz=0;xyz<nxyz;xyz++) {
			result4[xyz+X*nxyz] = direction[X][xyz];
			result4[xyz+Y*nxyz] = direction[Y][xyz];
			result4[xyz+Z*nxyz] = direction[Z][xyz];
		}					
		directionImage = result4;
		direction = null;
		
		correctImage = correction;
		sizeImage = lengthmap;
		
		return;
	}
	
	private final void directionFromRecursiveRidgeFilter0D(float[] img, boolean[] mask, float[] filter, byte[] direction, boolean unidirectional) {
			
			// get the tubular filter response
			float[][][] planescore = new float[nx][ny][nz];
			byte[][][] planedir = new byte[nx][ny][nz];
			byte[][][] linedir = new byte[nx][ny][nz];
			float[][][] inputImage = new float[nx][ny][nz];
			//float[] filter = new float[nx*ny*nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				inputImage[x][y][z] = img[xyz];
			}
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				filter[xyz] = 0.0f;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					// check for zero-valued neighbors as well
					minmaxplaneScore(inputImage, planescore, planedir, x,y,z, 13);
					// sign issue: remove all that have different sign, keep global sign
					float linescore = minmaxlineScore(inputImage, planedir, linedir, x,y,z, 4);
					if (planescore[x][y][z]*linescore>0) {
					    // now measure along the line
					    byte[] dl = directionNeighbor(linedir[x][y][z]);
					    float pointscore = Numerics.minmag(inputImage[x][y][z]-inputImage[x+dl[X]][y+dl[Y]][z+dl[Z]],inputImage[x][y][z]-inputImage[x-dl[X]][y-dl[Y]][z-dl[Z]]);
					    if (pointscore*linescore>0) {
                            filter[xyz] = (float)FastMath.cbrt(planescore[x][y][z]*linescore*pointscore);
                            direction[xyz] = -1;
                        } else {
                            filter[xyz] = 0.0f;
                            direction[xyz] = -1;
                        }
						if(filter[xyz]<0) if (unidirectional) { filter[xyz]=0; direction[xyz] = -1; } else filter[xyz]*=-1.0f;
					} else {
						filter[xyz] = 0.0f;
						direction[xyz] = -1;
					}
				}
			}
			planescore = null;
			planedir = null;
			linedir = null;
			inputImage = null;
			return;
	}
	private final void directionFromRecursiveRidgeFilter1D(float[] img, boolean[] mask, float[] filter, byte[] direction, boolean unidirectional) {
			
			// get the tubular filter response
			float[][][] planescore = new float[nx][ny][nz];
			byte[][][] planedir = new byte[nx][ny][nz];
			byte[][][] linedir = new byte[nx][ny][nz];
			float[][][] inputImage = new float[nx][ny][nz];
			//float[] filter = new float[nx*ny*nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				inputImage[x][y][z] = img[xyz];
			}
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				filter[xyz] = 0.0f;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					// check for zero-valued neighbors as well
					minmaxplaneScore(inputImage, planescore, planedir, x,y,z, 13);
					// sign issue: remove all that have different sign, keep global sign
					float linescore = minmaxlineScore(inputImage, planedir, linedir, x,y,z, 4);
					if (planescore[x][y][z]*linescore>0) {
						filter[xyz] = Numerics.sign(linescore)*Numerics.sqrt(planescore[x][y][z]*linescore);
						direction[xyz] = linedir[x][y][z];
						if(filter[xyz]<0) if (unidirectional) { filter[xyz]=0; direction[xyz] = -1; } else filter[xyz]*=-1.0f;
					} else {
						filter[xyz] = 0.0f;
						direction[xyz] = -1;
					}
				}
			}
			planescore = null;
			planedir = null;
			linedir = null;
			inputImage = null;
			return;
	}
	private final void directionFromRecursiveRidgeFilter2D(float[] img, boolean[] mask, float[] filter,byte[] direction, boolean unidirectional) {
			
			// get the tubular filter response
			float[][][] planescore = new float[nx][ny][nz];
			byte[][][] planedir = new byte[nx][ny][nz];
			float[][][] inputImage = new float[nx][ny][nz];
			//float[] filter = new float[nx*ny*nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				inputImage[x][y][z] = img[xyz];
			}
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				filter[xyz] = 0.0f;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					// check for zero-valued neighbors as well
					minmaxplaneScore(inputImage, planescore, planedir, x,y,z, 13);
					filter[xyz] = planescore[x][y][z];
					direction[xyz] = planedir[x][y][z];
					if (filter[xyz]<0) if (unidirectional) { filter[xyz]=0; direction[xyz] = -1; } else filter[xyz]*=-1.0f;
				}
			}
			planescore = null;
			planedir = null;
			inputImage = null;
			return;
	}
	private final void probabilityFromRecursiveRidgeFilter( float[] filter, float[] shape) {
		// normalization: best is the iterParamative robust exponential (others are removed)
		int nb = 0;
		double min = 1e9;
		double max = 0.0;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
				// keep only the proper sign
				if (filter[xyz]<=0) filter[xyz] = 0.0f;
				else {
					// fit exp only to non-zero data
					nb++;
					min = Numerics.min(filter[xyz], min);
					max = Numerics.max(filter[xyz], max);
				}
		}
		
		// robust measures? pb is the exponential is not steep enough
		double[] response = new double[nb];
		int n=0;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
				if (filter[xyz]>0) {
					response[n] = filter[xyz];
					n++;
				}
		}
		Percentile measure = new Percentile();
		double median = measure.evaluate(response, 50.0);
		double beta = median/FastMath.log(2.0);
		
		BasicInfo.displayMessage("exponential parameter estimates: median "+median+", beta "+beta+",\n");
		
		// model the filter response as something more interesting, e.g. log-normal (removing the bg samples)
		double[] weights = new double[nb];
		for (int b=0;b<nb;b++) { 
			weights[b] = (1.0-FastMath.exp( -response[b]/beta));
			response[b] = FastMath.log(response[b]);
		}
		
		double fmean = ImageStatistics.weightedPercentile(response,weights,50.0,nb);
		
		// stdev: 50% +/- 1/2*erf(1/sqrt(2)) (~0.341344746..)
		double dev = 100.0*0.5*Erf.erf(1.0/FastMath.sqrt(2.0));
		double fdev = 0.5*(ImageStatistics.weightedPercentile(response,weights,50.0+dev,nb) - ImageStatistics.weightedPercentile(response,weights,50.0-dev,nb));
		
		BasicInfo.displayMessage("Log-normal parameter estimates: mean = "+FastMath.exp(fmean)+", stdev = "+FastMath.exp(fdev)+",\n");
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			int xyz = x + nx*y + nx*ny*z;
			if (filter[xyz]>0) {
				double pe = FastMath.exp( -filter[xyz]/beta)/beta;
				double plg = FastMath.exp(-Numerics.square(FastMath.log(filter[xyz])-fmean)/(2.0*fdev*fdev))/FastMath.sqrt(2.0*FastMath.PI*fdev*fdev);
				shape[xyz] = (float)(plg/(plg+pe));
				//shape[xyz] = (float)(1.0-pe);
			}
		}
		return;
	}
	boolean zeroNeighbor(float[] inputImage, boolean[] mask, int x, int y, int z, int d) {
			for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int l=-d;l<=d;l++) {
				if (inputImage[x+i+nx*(y+j)+nx*ny*(z+l)]!=inputImage[x+nx*y+nx*ny*z] && i*i+j*j+l*l<=2*d*d) return false;
			}
			return true;
		}
	void minmaxplaneScore(float[][][] inputImage, float[][][] plane, byte[][][] dir, int x, int y, int z, int dmax) {
		float maxgrad = 0.0f; 
		float minval = 0.0f; 
		float sign = 0.0f;
		byte direction = -1;
		for (byte d=0;d<dmax;d++) {
			float val1 = 0.0f, val2 = 0.0f;
			if (d==0) {		
				val1=(inputImage[x][y][z]		-inputImage[x-1][y][z]
					 +inputImage[x][y-1][z]		-inputImage[x-1][y-1][z]
					 +inputImage[x][y+1][z]		-inputImage[x-1][y+1][z]
					 +inputImage[x][y][z-1]		-inputImage[x-1][y][z-1]
					 +inputImage[x][y][z+1]		-inputImage[x-1][y][z+1]
					 +inputImage[x][y-1][z-1]	-inputImage[x-1][y-1][z-1]
					 +inputImage[x][y-1][z+1]	-inputImage[x-1][y-1][z+1]
					 +inputImage[x][y+1][z-1]	-inputImage[x-1][y+1][z-1]
					 +inputImage[x][y+1][z+1]	-inputImage[x-1][y+1][z+1])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x+1][y][z]
					 +inputImage[x][y-1][z]		-inputImage[x+1][y-1][z]
					 +inputImage[x][y+1][z]		-inputImage[x+1][y+1][z]
					 +inputImage[x][y][z-1]		-inputImage[x+1][y][z-1]
					 +inputImage[x][y][z+1]		-inputImage[x+1][y][z+1]
					 +inputImage[x][y-1][z-1]	-inputImage[x+1][y-1][z-1]
					 +inputImage[x][y-1][z+1]	-inputImage[x+1][y-1][z+1]
					 +inputImage[x][y+1][z-1]	-inputImage[x+1][y+1][z-1]
					 +inputImage[x][y+1][z+1]	-inputImage[x+1][y+1][z+1])/9.0f;
			} else if (d==1) {
				val1=(inputImage[x][y][z]		-inputImage[x][y-1][z]
					 +inputImage[x-1][y][z]		-inputImage[x-1][y-1][z]
					 +inputImage[x+1][y][z]		-inputImage[x+1][y-1][z]	
					 +inputImage[x][y][z-1]		-inputImage[x][y-1][z-1]	
					 +inputImage[x][y][z+1]		-inputImage[x][y-1][z+1]
					 +inputImage[x-1][y][z-1]	-inputImage[x-1][y-1][z-1]
					 +inputImage[x-1][y][z+1]	-inputImage[x-1][y-1][z+1]
					 +inputImage[x+1][y][z-1]	-inputImage[x+1][y-1][z-1]
					 +inputImage[x+1][y][z+1]	-inputImage[x+1][y-1][z+1])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x][y+1][z]
					 +inputImage[x-1][y][z]		-inputImage[x-1][y+1][z]
					 +inputImage[x+1][y][z]		-inputImage[x+1][y+1][z]
					 +inputImage[x][y][z-1]		-inputImage[x][y+1][z-1]
					 +inputImage[x][y][z+1]		-inputImage[x][y+1][z+1]
					 +inputImage[x-1][y][z-1]	-inputImage[x-1][y+1][z-1]
					 +inputImage[x-1][y][z+1]	-inputImage[x-1][y+1][z+1]
					 +inputImage[x+1][y][z-1]	-inputImage[x+1][y+1][z-1]
					 +inputImage[x+1][y][z+1]	-inputImage[x+1][y+1][z+1])/9.0f;
			} else if (d==2) { 			
				val1=(inputImage[x][y][z]		-inputImage[x][y][z-1]
					 +inputImage[x-1][y][z]		-inputImage[x-1][y][z-1]
					 +inputImage[x+1][y][z]		-inputImage[x+1][y][z-1]
					 +inputImage[x][y-1][z]		-inputImage[x][y-1][z-1]
					 +inputImage[x][y+1][z]		-inputImage[x][y+1][z-1]	
					 +inputImage[x-1][y-1][z]	-inputImage[x-1][y-1][z-1]
					 +inputImage[x-1][y+1][z]	-inputImage[x-1][y+1][z-1]
					 +inputImage[x+1][y-1][z]	-inputImage[x+1][y-1][z-1]
					 +inputImage[x+1][y+1][z]	-inputImage[x+1][y+1][z-1])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x][y][z+1]
					 +inputImage[x-1][y][z]		-inputImage[x-1][y][z+1]
					 +inputImage[x+1][y][z]		-inputImage[x+1][y][z+1]
					 +inputImage[x][y-1][z]		-inputImage[x][y-1][z+1]
					 +inputImage[x][y+1][z]		-inputImage[x][y+1][z+1]
					 +inputImage[x-1][y-1][z]	-inputImage[x-1][y-1][z+1]
					 +inputImage[x-1][y+1][z]	-inputImage[x-1][y+1][z+1]
					 +inputImage[x+1][y-1][z]	-inputImage[x+1][y-1][z+1]
					 +inputImage[x+1][y+1][z]	-inputImage[x+1][y+1][z+1])/9.0f;
			} else if (d==3) { 			
				val1=(inputImage[x][y][z]		-inputImage[x-1][y-1][z]
					 +inputImage[x-1][y+1][z]	-inputImage[x-2][y][z]
					 +inputImage[x+1][y-1][z]	-inputImage[x][y-2][z]
					 +inputImage[x][y][z-1]		-inputImage[x-1][y-1][z-1]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x-2][y][z-1]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x][y-2][z-1]
					 +inputImage[x][y][z+1]		-inputImage[x-1][y-1][z+1]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x-2][y][z+1]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x][y-2][z+1])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x+1][y+1][z]
					 +inputImage[x-1][y+1][z]	-inputImage[x][y+2][z]
					 +inputImage[x+1][y-1][z]	-inputImage[x+2][y][z]
					 +inputImage[x][y][z-1]		-inputImage[x+1][y+1][z-1]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x][y+2][z-1]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x+2][y][z-1]
					 +inputImage[x][y][z+1]		-inputImage[x+1][y+1][z+1]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x][y+2][z+1]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x+2][y][z+1])/9.0f;
			} else if (d==4) { 			
				val1=(inputImage[x][y][z]		-inputImage[x][y-1][z-1]
					 +inputImage[x][y+1][z-1]	-inputImage[x][y][z-2]
					 +inputImage[x][y-1][z+1]	-inputImage[x][y-2][z]
					 +inputImage[x-1][y][z]		-inputImage[x-1][y-1][z-1]
					 +inputImage[x-1][y+1][z-1]-inputImage[x-1][y][z-2]
					 +inputImage[x-1][y-1][z+1]-inputImage[x-1][y-2][z]
					 +inputImage[x+1][y][z]		-inputImage[x+1][y-1][z-1]
					 +inputImage[x+1][y+1][z-1]-inputImage[x+1][y][z-2]
					 +inputImage[x+1][y-1][z+1]-inputImage[x+1][y-2][z])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x][y+1][z+1]
					 +inputImage[x][y+1][z-1]	-inputImage[x][y+2][z]
					 +inputImage[x][y-1][z+1]	-inputImage[x][y][z+2]
					 +inputImage[x-1][y][z]		-inputImage[x-1][y+1][z+1]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x-1][y+2][z]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x-1][y][z+2]
					 +inputImage[x+1][y][z]		-inputImage[x+1][y+1][z+1]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x+1][y+2][z]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x+1][y][z+2])/9.0f;
			} else if (d==5) { 			
				val1=(inputImage[x][y][z]		-inputImage[x-1][y][z-1]
					 +inputImage[x+1][y][z-1]	-inputImage[x][y][z-2]
					 +inputImage[x-1][y][z+1]	-inputImage[x-2][y][z]
					 +inputImage[x][y-1][z]		-inputImage[x-1][y-1][z-1]
					 +inputImage[x+1][y-1][z-1]-inputImage[x][y-1][z-2]
					 +inputImage[x-1][y-1][z+1]-inputImage[x-2][y-1][z]
					 +inputImage[x][y+1][z]		-inputImage[x-1][y+1][z-1]
					 +inputImage[x+1][y+1][z-1]-inputImage[x][y+1][z-2]
					 +inputImage[x-1][y+1][z+1]-inputImage[x-2][y+1][z])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x+1][y][z+1]
					 +inputImage[x+1][y][z-1]	-inputImage[x+2][y][z]
					 +inputImage[x-1][y][z+1]	-inputImage[x][y][z+2]
					 +inputImage[x][y-1][z]		-inputImage[x+1][y-1][z+1]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x+2][y-1][z]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x][y-1][z+2]
					 +inputImage[x][y+1][z]		-inputImage[x+1][y+1][z+1]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x+2][y+1][z]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x][y+1][z+2])/9.0f;
			} else if (d==6) { 			
				val1=(inputImage[x][y][z]		-inputImage[x-1][y+1][z]
					 +inputImage[x-1][y-1][z]	-inputImage[x-2][y][z]
					 +inputImage[x+1][y+1][z]	-inputImage[x][y-2][z]
					 +inputImage[x][y][z-1]		-inputImage[x-1][y+1][z-1]
					 +inputImage[x-1][y-1][z-1]-inputImage[x-2][y][z-1]
					 +inputImage[x+1][y+1][z-1]-inputImage[x][y-2][z-1]
					 +inputImage[x][y][z+1]		-inputImage[x-1][y+1][z+1]
					 +inputImage[x-1][y-1][z+1]-inputImage[x-2][y][z+1]
					 +inputImage[x+1][y+1][z+1]-inputImage[x][y-2][z+1])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x+1][y-1][z]
					 +inputImage[x-1][y-1][z]	-inputImage[x][y-2][z]
					 +inputImage[x+1][y+1][z]	-inputImage[x+2][y][z]
					 +inputImage[x][y][z-1]		-inputImage[x+1][y-1][z-1]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x][y-2][z-1]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x+2][y][z-1]
					 +inputImage[x][y][z+1]		-inputImage[x+1][y-1][z+1]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x][y-2][z+1]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x+2][y][z+1])/9.0f;
			} else if (d==7) { 			
				val1=(inputImage[x][y][z]		-inputImage[x][y-1][z+1]
					 +inputImage[x][y-1][z-1]	-inputImage[x][y-2][z]
					 +inputImage[x][y+1][z+1]	-inputImage[x][y][z+2]
					 +inputImage[x-1][y][z]		-inputImage[x-1][y-1][z+1]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x-1][y-2][z]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x-1][y][z+2]
					 +inputImage[x+1][y][z]		-inputImage[x+1][y-1][z+1]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x+1][y-2][z]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x+1][y][z+2])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x][y+1][z-1]
					 +inputImage[x][y-1][z-1]	-inputImage[x][y][z-2]
					 +inputImage[x][y+1][z+1]	-inputImage[x][y+2][z]
					 +inputImage[x-1][y][z]		-inputImage[x-1][y+1][z-1]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x-1][y][z-2]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x-1][y+2][z]
					 +inputImage[x+1][y][z]		-inputImage[x+1][y+1][z-1]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x+1][y][z-2]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x+1][y+2][z])/9.0f;
			} else if (d==8) { 			
				val1=(inputImage[x][y][z]		-inputImage[x-1][y][z+1]
					 +inputImage[x-1][y][z-1]	-inputImage[x-2][y][z]
					 +inputImage[x+1][y][z+1]	-inputImage[x][y][z+2]
					 +inputImage[x][y-1][z]		-inputImage[x-1][y-1][z+1]
					 +inputImage[x-1][y-1][z-1]-inputImage[x-2][y-1][z]
					 +inputImage[x+1][y-1][z+1]-inputImage[x][y-1][z+2]
					 +inputImage[x][y+1][z]		-inputImage[x-1][y+1][z+1]
					 +inputImage[x-1][y+1][z-1]-inputImage[x-2][y+1][z]
					 +inputImage[x+1][y+1][z+1]-inputImage[x][y+1][z+2])/9.0f;
				val2=(inputImage[x][y][z]		-inputImage[x+1][y][z-1]
					 +inputImage[x-1][y][z-1]	-inputImage[x][y][z-2]
					 +inputImage[x+1][y][z+1]	-inputImage[x+2][y][z]
					 +inputImage[x][y-1][z]		-inputImage[x+1][y-1][z-1]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x][y-1][z-2]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x+2][y-1][z]
					 +inputImage[x][y+1][z]		-inputImage[x+1][y+1][z-1]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x][y+1][z-2]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x+2][y+1][z])/9.0f;
			} else if (d==9) { 			
				val1=(inputImage[x][y][z]		-inputImage[x-1][y-1][z-1]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x-2][y-2][z]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x-2][y][z-2]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x][y-2][z-2]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x][y][z-2]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x-2][y][z]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x][y-2][z])/7.0f;
				val2=(inputImage[x][y][z]		-inputImage[x+1][y+1][z+1]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x][y][z+2]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x][y+2][z]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x+2][y][z]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x+2][y+2][z]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x][y+2][z+2]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x+2][y][z+2])/7.0f;
			} else if (d==10) { 			
				val1=(inputImage[x][y][z]		-inputImage[x+1][y-1][z-1]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x+2][y-2][z]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x+2][y][z-2]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x][y-2][z-2]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x][y][z-2]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x][y-2][z]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x+2][y][z])/7.0f;
				val2=(inputImage[x][y][z]		-inputImage[x-1][y+1][z+1]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x][y][z+2]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x][y+2][z]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x-2][y][z]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x-2][y+2][z]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x-2][y][z+2]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x][y+2][z+2])/7.0f;
			} else if (d==11) { 			
				val1=(inputImage[x][y][z]		-inputImage[x-1][y+1][z-1]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x-2][y+2][z]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x][y+2][z-2]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x-2][y][z-2]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x][y][z-2]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x-2][y][z]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x][y+2][z])/7.0f;
				val2=(inputImage[x][y][z]		-inputImage[x+1][y-1][z+1]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x][y][z+2]
					 +inputImage[x+1][y+1][z-1]	-inputImage[x+2][y][z]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x][y-2][z]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x+2][y-2][z]
					 +inputImage[x-1][y-1][z+1]	-inputImage[x][y-2][z+2]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x+2][y][z+2])/7.0f;
			} else if (d==12) { 			
				val1=(inputImage[x][y][z]		-inputImage[x-1][y-1][z+1]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x-2][y][z+2]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x][y-2][z+2]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x-2][y-2][z]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x][y-2][z]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x-2][y][z]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x][y][z+2])/7.0f;
				val2=(inputImage[x][y][z]		-inputImage[x+1][y+1][z-1]
					 +inputImage[x-1][y+1][z+1]	-inputImage[x][y+2][z]
					 +inputImage[x+1][y-1][z+1]	-inputImage[x+2][y][z]
					 +inputImage[x-1][y-1][z-1]	-inputImage[x][y][z-2]
					 +inputImage[x+1][y-1][z-1]	-inputImage[x+2][y][z-2]
					 +inputImage[x-1][y+1][z-1]	-inputImage[x][y+2][z-2]
					 +inputImage[x+1][y+1][z+1]	-inputImage[x+2][y+2][z])/7.0f;
			}
			// find the strongest gradient direction, then estimate the corresponding filter response
			if (val1*val1+val2*val2>maxgrad) {
				maxgrad = val1*val1+val2*val2;
				if (val1*val1<val2*val2) minval = val1;
				else minval = val2;
				sign = val1*val2;
				direction = d;
			}
		}
		if (sign>0) {
			plane[x][y][z] = minval;
			dir[x][y][z] = direction;
		} else {
			plane[x][y][z] = 0.0f;
			dir[x][y][z] = -1;
		}
		return;
		
	}
	
	float minmaxlineScore(float[][][] inputImage, byte[][][] planedir, byte[][][] linedir, int x, int y, int z, int dmax) {
		float maxgrad = 0.0f; 
		float minval = 0.0f; 
		float sign = 0.0f;
		byte direction = -1;
		
		float val1 = 0.0f, val2 = 0.0f;
		if (planedir[x][y][z]==0) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z]	+inputImage[x][y+1][z]
							-inputImage[x][y][z-1]	-inputImage[x][y-1][z-1]	-inputImage[x][y+1][z-1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z]	+inputImage[x][y+1][z]
							-inputImage[x][y][z+1]	-inputImage[x][y-1][z+1]	-inputImage[x][y+1][z+1];
					direction = Y;
				} else if (d==1) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y][z-1]	+inputImage[x][y][z+1]
							-inputImage[x][y-1][z]	-inputImage[x][y-1][z-1]	-inputImage[x][y-1][z+1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y][z-1]	+inputImage[x][y][z+1]
							-inputImage[x][y+1][z]	-inputImage[x][y+1][z-1]	-inputImage[x][y+1][z+1];
					direction = Z;		
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z-1]	+inputImage[x][y+1][z+1]
							-inputImage[x][y-1][z+1]	-inputImage[x][y-2][z]	-inputImage[x][y][z+2];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z-1]	+inputImage[x][y+1][z+1]
							-inputImage[x][y+1][z-1]	-inputImage[x][y][z-2]	-inputImage[x][y+2][z];
					direction = YpZ;		
				} else if (d==3) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z+1]	+inputImage[x][y+1][z-1]
							-inputImage[x][y-1][z-1]	-inputImage[x][y-2][z]	-inputImage[x][y][z-2];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y+1][z-1]	+inputImage[x][y-1][z+1]
							-inputImage[x][y+1][z+1]	-inputImage[x][y+2][z]	-inputImage[x][y][z+2];
					direction = YmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==1) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y][z-1]	+inputImage[x][y][z+1]
							-inputImage[x-1][y][z]	-inputImage[x-1][y][z-1]	-inputImage[x-1][y][z+1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y][z-1]	+inputImage[x][y][z+1]
							-inputImage[x+1][y][z]	-inputImage[x+1][y][z-1]	-inputImage[x+1][y][z+1];
					direction = Z;		
				} else if (d==1) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z]	+inputImage[x+1][y][z]
							-inputImage[x][y][z-1]	-inputImage[x-1][y][z-1]	-inputImage[x+1][y][z-1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z]	+inputImage[x+1][y][z]
							-inputImage[x][y][z+1]	-inputImage[x-1][y][z+1]	-inputImage[x+1][y][z+1];
					direction = X;		
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z-1]	+inputImage[x+1][y][z+1]
							-inputImage[x-1][y][z+1]	-inputImage[x-2][y][z]	-inputImage[x][y][z+2];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z-1]	+inputImage[x+1][y][z+1]
							-inputImage[x+1][y][z-1]	-inputImage[x][y][z-2]	-inputImage[x+2][y][z];
					direction = ZpX;		
				} else if (d==3) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z+1]	+inputImage[x+1][y][z-1]
							-inputImage[x-1][y][z-1]	-inputImage[x-2][y][z]	-inputImage[x][y][z-2];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z+1]	+inputImage[x+1][y][z-1]
							-inputImage[x+1][y][z+1]	-inputImage[x][y][z+2]	-inputImage[x+2][y][z];
					direction = ZmX;		
				}	
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==2) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z]	+inputImage[x+1][y][z]
							-inputImage[x][y-1][z]	-inputImage[x-1][y-1][z]	-inputImage[x+1][y-1][z];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z]	+inputImage[x+1][y][z]
							-inputImage[x][y+1][z]	-inputImage[x-1][y+1][z]	-inputImage[x+1][y+1][z];
					direction = X;		
				} else if (d==1) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z]	+inputImage[x][y+1][z]
							-inputImage[x-1][y][z]	-inputImage[x-1][y-1][z]	-inputImage[x-1][y+1][z];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z]	+inputImage[x][y+1][z]
							-inputImage[x+1][y][z]	-inputImage[x+1][y-1][z]	-inputImage[x+1][y+1][z];
					direction = Y;		
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y-1][z]	+inputImage[x+1][y+1][z]
							-inputImage[x-1][y+1][z]	-inputImage[x-2][y][z]	-inputImage[x][y+2][z];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y-1][z]	+inputImage[x+1][y+1][z]
							-inputImage[x+1][y-1][z]	-inputImage[x][y-2][z]	-inputImage[x+2][y][z];
					direction = XpY;		
				} else if (d==3) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y+1][z]	+inputImage[x+1][y-1][z]
							-inputImage[x-1][y-1][z]	-inputImage[x-2][y][z]	-inputImage[x][y-2][z];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y+1][z]	+inputImage[x+1][y-1][z]
							-inputImage[x+1][y+1][z]	-inputImage[x][y+2][z]	-inputImage[x+2][y][z];
					direction = XmY;		
				}	
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==3) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y+1][z]		+inputImage[x+1][y-1][z]
							-inputImage[x][y][z-1]	-inputImage[x-1][y+1][z-1]	-inputImage[x+1][y-1][z-1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y+1][z]		+inputImage[x+1][y-1][z]
							-inputImage[x][y][z+1]	-inputImage[x-1][y+1][z+1]	-inputImage[x+1][y-1][z+1];
					direction = XmY;		
				} else if (d==1) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y][z-1]		+inputImage[x][y][z+1]
							-inputImage[x-1][y+1][z]	-inputImage[x-1][y+1][z-1]	-inputImage[x-1][y+1][z+1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y][z-1]		+inputImage[x][y][z+1]
							-inputImage[x+1][y-1][z]	-inputImage[x+1][y-1][z-1]	-inputImage[x+1][y-1][z+1];
					direction = Z;		
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x-1][y+1][z-1]	+inputImage[x+1][y-1][z+1]
							-inputImage[x-1][y+1][z+1]	-inputImage[x-2][y+2][z]		-inputImage[x][y][z+2];
					val2 = 	 inputImage[x][y][z]			+inputImage[x-1][y+1][z-1]	+inputImage[x+1][y-1][z+1]
							-inputImage[x+1][y-1][z-1]	-inputImage[x][y][z-2]		-inputImage[x+2][y-2][z];
					direction = XmYpZ;		
				} else if (d==3) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x-1][y+1][z+1]	+inputImage[x+1][y-1][z-1]
							-inputImage[x-1][y+1][z-1]	-inputImage[x-2][y+2][z]		-inputImage[x][y][z-2];
					val2 = 	 inputImage[x][y][z]			+inputImage[x-1][y+1][z+1]	+inputImage[x+1][y-1][z-1]
							-inputImage[x+1][y-1][z+1]	-inputImage[x][y][z+2]		-inputImage[x+2][y-2][z];
					direction = XmYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==4) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y+1][z-1]		+inputImage[x][y-1][z+1]
							-inputImage[x-1][y][z]	-inputImage[x-1][y+1][z-1]	-inputImage[x-1][y-1][z+1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y+1][z-1]		+inputImage[x][y-1][z+1]
							-inputImage[x+1][y][z]	-inputImage[x+1][y+1][z-1]	-inputImage[x+1][y-1][z+1];
					direction = YmZ;		
				} else if (d==1) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z]		+inputImage[x+1][y][z]
							-inputImage[x][y+1][z-1]	-inputImage[x-1][y+1][z-1]	-inputImage[x+1][y+1][z-1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y][z]		+inputImage[x+1][y][z]
							-inputImage[x][y-1][z+1]	-inputImage[x-1][y-1][z+1]	-inputImage[x+1][y-1][z+1];
					direction = X;		
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x-1][y+1][z-1]	+inputImage[x+1][y-1][z+1]
							-inputImage[x-1][y-1][z+1]	-inputImage[x-2][y][z]		-inputImage[x][y-2][z+2];
					val2 = 	 inputImage[x][y][z]			+inputImage[x-1][y+1][z-1]	+inputImage[x+1][y-1][z+1]
							-inputImage[x+1][y+1][z-1]	-inputImage[x][y+2][z-2]		-inputImage[x+2][y][z];
					direction = XmYpZ;		
				} else if (d==3) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x+1][y+1][z-1]	+inputImage[x-1][y-1][z+1]
							-inputImage[x+1][y-1][z+1]	-inputImage[x+2][y][z]		-inputImage[x][y-2][z+2];
					val2 = 	 inputImage[x][y][z]			+inputImage[x+1][y+1][z-1]	+inputImage[x-1][y-1][z+1]
							-inputImage[x-1][y+1][z-1]	-inputImage[x][y+2][z-2]		-inputImage[x-2][y][z];
					direction = XpYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==5) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x+1][y][z-1]		+inputImage[x-1][y][z+1]
							-inputImage[x][y-1][z]	-inputImage[x+1][y-1][z-1]	-inputImage[x-1][y-1][z+1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x+1][y][z-1]		+inputImage[x-1][y][z+1]
							-inputImage[x][y+1][z]	-inputImage[x+1][y+1][z-1]	-inputImage[x-1][y+1][z+1];
					direction = ZmX;			
				} else if (d==1) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z]		+inputImage[x][y+1][z]
							-inputImage[x+1][y][z-1]	-inputImage[x+1][y-1][z-1]	-inputImage[x+1][y+1][z-1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y-1][z]		+inputImage[x][y+1][z]
							-inputImage[x-1][y][z+1]	-inputImage[x-1][y-1][z+1]	-inputImage[x-1][y+1][z+1];
					direction = Y;		
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x+1][y-1][z-1]	+inputImage[x-1][y+1][z+1]
							-inputImage[x-1][y-1][z+1]	-inputImage[x][y-2][z]		-inputImage[x-2][y][z+2];
					val2 = 	 inputImage[x][y][z]			+inputImage[x+1][y-1][z-1]	+inputImage[x-1][y+1][z+1]
							-inputImage[x+1][y+1][z-1]	-inputImage[x+2][y][z-2]		-inputImage[x][y+2][z];
					direction = XmYmZ;		
				} else if (d==3) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x+1][y+1][z-1]	+inputImage[x-1][y-1][z+1]
							-inputImage[x-1][y+1][z+1]	-inputImage[x][y+2][z]		-inputImage[x-2][y][z+2];
					val2 = 	 inputImage[x][y][z]			+inputImage[x+1][y+1][z-1]	+inputImage[x-1][y-1][z+1]
							-inputImage[x+1][y-1][z-1]	-inputImage[x+2][y][z-2]		-inputImage[x][y-2][z];
					direction = XpYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==6) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x-1][y-1][z]		+inputImage[x+1][y+1][z]
							-inputImage[x][y][z-1]	-inputImage[x-1][y-1][z-1]	-inputImage[x+1][y+1][z-1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x-1][y-1][z]		+inputImage[x+1][y+1][z]
							-inputImage[x][y][z+1]	-inputImage[x-1][y-1][z+1]	-inputImage[x+1][y+1][z+1];
					direction = XpY;		
				} else if (d==1) {
					val1 = 	 inputImage[x][y][z]		+inputImage[x][y][z-1]		+inputImage[x][y][z+1]
							-inputImage[x-1][y-1][z]	-inputImage[x-1][y-1][z-1]	-inputImage[x-1][y-1][z+1];
					val2 = 	 inputImage[x][y][z]		+inputImage[x][y][z-1]		+inputImage[x][y][z+1]
							-inputImage[x+1][y+1][z]	-inputImage[x+1][y+1][z-1]	-inputImage[x+1][y+1][z+1];
					direction = Z;		
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x-1][y-1][z-1]	+inputImage[x+1][y+1][z+1]
							-inputImage[x-1][y-1][z+1]	-inputImage[x-2][y-2][z]		-inputImage[x][y][z+2];
					val2 = 	 inputImage[x][y][z]			+inputImage[x-1][y-1][z-1]	+inputImage[x+1][y+1][z+1]
							-inputImage[x+1][y+1][z-1]	-inputImage[x][y][z-2]		-inputImage[x+2][y+2][z];
					direction = XpYpZ;		
				} else if (d==3) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x-1][y-1][z+1]	+inputImage[x+1][y+1][z-1]
							-inputImage[x-1][y-1][z-1]	-inputImage[x-2][y-2][z]		-inputImage[x][y][z-2];
					val2 = 	 inputImage[x][y][z]			+inputImage[x-1][y-1][z+1]	+inputImage[x+1][y+1][z-1]
							-inputImage[x+1][y+1][z+1]	-inputImage[x][y][z+2]		-inputImage[x+2][y+2][z];
					direction = XpYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==7) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 inputImage[x][y][z]		+inputImage[x][y-1][z-1]		+inputImage[x][y+1][z+1]
							-inputImage[x-1][y][z]	-inputImage[x-1][y-1][z-1]	-inputImage[x-1][y+1][z+1];
					val2 =	 inputImage[x][y][z]		+inputImage[x][y-1][z-1]		+inputImage[x][y+1][z+1]
							-inputImage[x+1][y][z]	-inputImage[x+1][y-1][z-1]	-inputImage[x+1][y+1][z+1];
					direction = YpZ;		
				} else if (d==1) {
					val1 =	 inputImage[x][y][z]		+inputImage[x-1][y][z]		+inputImage[x+1][y][z]
							-inputImage[x][y-1][z-1]	-inputImage[x-1][y-1][z-1]	-inputImage[x+1][y-1][z-1];
					val2 =	 inputImage[x][y][z]		+inputImage[x-1][y][z]		+inputImage[x+1][y][z]
							-inputImage[x][y+1][z+1]	-inputImage[x-1][y+1][z+1]	-inputImage[x+1][y+1][z+1];
					direction = X;		
				} else if (d==2) {
					val1 =   inputImage[x][y][z]			+inputImage[x-1][y-1][z-1]	+inputImage[x+1][y+1][z+1]
							-inputImage[x+1][y-1][z-1]	-inputImage[x][y-2][z-2]		-inputImage[x+2][y][z];
					val2 =	 inputImage[x][y][z]			+inputImage[x-1][y-1][z-1]	+inputImage[x+1][y+1][z+1]
							-inputImage[x-1][y+1][z+1]	-inputImage[x-2][y][z]		-inputImage[x][y+2][z+2];
					direction = XpYpZ;		
				} else if (d==3) {
					val1 =   inputImage[x][y][z]			+inputImage[x+1][y-1][z-1]	+inputImage[x-1][y+1][z+1]
							-inputImage[x-1][y-1][z-1]	-inputImage[x][y-2][z-2]		-inputImage[x-2][y][z];
					val2 =	 inputImage[x][y][z]			+inputImage[x+1][y-1][z-1]	+inputImage[x-1][y+1][z+1]
							-inputImage[x+1][y+1][z+1]	-inputImage[x+2][y][z]		-inputImage[x][y+2][z+2];
					direction = XmYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==8) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 inputImage[x][y][z]		+inputImage[x-1][y][z-1]		+inputImage[x+1][y][z+1]
							-inputImage[x][y-1][z]	-inputImage[x-1][y-1][z-1]	-inputImage[x+1][y-1][z+1];
					val2 =	 inputImage[x][y][z]		+inputImage[x-1][y][z-1]		+inputImage[x+1][y][z+1]
							-inputImage[x][y+1][z]	-inputImage[x-1][y+1][z-1]	-inputImage[x+1][y+1][z+1];
					direction = ZpX;		
				} else if (d==1) {
					val1 =	 inputImage[x][y][z]		+inputImage[x][y-1][z]		+inputImage[x][y+1][z]
							-inputImage[x-1][y][z-1]	-inputImage[x-1][y-1][z-1]	-inputImage[x-1][y+1][z-1];
					val2 =	 inputImage[x][y][z]		+inputImage[x][y-1][z]		+inputImage[x][y+1][z]
							-inputImage[x+1][y][z+1]	-inputImage[x+1][y-1][z+1]	-inputImage[x+1][y+1][z+1];
					direction = Y;		
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x-1][y-1][z-1]	+inputImage[x+1][y+1][z+1]
							-inputImage[x-1][y+1][z-1]	-inputImage[x-2][y][z-2]		-inputImage[x][y+2][z];
					val2 = 	 inputImage[x][y][z]			+inputImage[x-1][y-1][z-1]	+inputImage[x+1][y+1][z+1]
							-inputImage[x+1][y-1][z+1]	-inputImage[x][y-2][z]		-inputImage[x+2][y][z+2];
					direction = XpYpZ;		
				} else if (d==3) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x-1][y+1][z-1]	+inputImage[x+1][y-1][z+1]
							-inputImage[x-1][y-1][z-1]	-inputImage[x-2][y][z-2]		-inputImage[x][y-2][z];
					val2 = 	 inputImage[x][y][z]			+inputImage[x-1][y+1][z-1]	+inputImage[x+1][y-1][z+1]
							-inputImage[x+1][y+1][z+1]	-inputImage[x][y+2][z]		-inputImage[x+2][y][z+2];
					direction = XmYpZ;		
				}					
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==9) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 inputImage[x][y][z]			+inputImage[x-1][y-1][z+1]	+inputImage[x+1][y+1][z-1]
											  -1.5f*(inputImage[x-1][y+1][z-1]	+inputImage[x-1][y+1][z+1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x-1][y-1][z+1]	+inputImage[x+1][y+1][z-1]
											  -1.5f*(inputImage[x+1][y-1][z+1]	+inputImage[x+1][y-1][z-1]);
					direction = XpYmZ;						  
				} else if (d==1) {
					val1 =	 inputImage[x][y][z]			+inputImage[x-1][y+1][z-1]	+inputImage[x+1][y-1][z+1]
											  -1.5f*(inputImage[x-1][y+1][z+1]	+inputImage[x-1][y-1][z+1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x-1][y+1][z-1]	+inputImage[x+1][y-1][z+1]
											  -1.5f*(inputImage[x+1][y-1][z-1]	+inputImage[x+1][y+1][z-1]);
					direction = XmYpZ;						  
				} else if (d==2) {
					val1 = 	 inputImage[x][y][z]			+inputImage[x+1][y-1][z-1]	+inputImage[x-1][y+1][z+1]
											  -1.5f*(inputImage[x-1][y-1][z+1]	+inputImage[x+1][y-1][z+1]);
					val2 = 	 inputImage[x][y][z]			+inputImage[x+1][y-1][z-1]	+inputImage[x-1][y+1][z+1]
											  -1.5f*(inputImage[x+1][y+1][z-1]	+inputImage[x-1][y+1][z-1]);
					direction = XmYmZ;						  
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==10) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 inputImage[x][y][z]			+inputImage[x+1][y-1][z+1]	+inputImage[x-1][y+1][z-1]
											  -1.5f*(inputImage[x+1][y+1][z-1]	+inputImage[x+1][y+1][z+1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x+1][y-1][z+1]	+inputImage[x-1][y+1][z-1]
											  -1.5f*(inputImage[x-1][y-1][z+1]	+inputImage[x-1][y-1][z-1]);
					direction = XmYpZ;						  
				} else 	if (d==1) {
					val1 =	 inputImage[x][y][z]			+inputImage[x+1][y+1][z-1]	+inputImage[x-1][y-1][z+1]
											  -1.5f*(inputImage[x-1][y-1][z-1]	+inputImage[x-1][y+1][z-1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x+1][y+1][z-1]	+inputImage[x-1][y-1][z+1]
											  -1.5f*(inputImage[x+1][y+1][z+1]	+inputImage[x+1][y-1][z+1]);
					direction = XpYmZ;						  
				} else 	if (d==2) {
					val1 =	 inputImage[x][y][z]			+inputImage[x+1][y+1][z+1]	+inputImage[x-1][y-1][z-1]
											  -1.5f*(inputImage[x+1][y-1][z+1]	+inputImage[x-1][y-1][z+1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x+1][y+1][z+1]	+inputImage[x-1][y-1][z-1]
											  -1.5f*(inputImage[x-1][y+1][z-1]	+inputImage[x+1][y+1][z-1]);
					direction = XpYpZ;						  
				}		
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}			
		} else if (planedir[x][y][z]==11) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 inputImage[x][y][z]			+inputImage[x-1][y+1][z+1]	+inputImage[x+1][y-1][z-1]
											  -1.5f*(inputImage[x+1][y+1][z-1]	+inputImage[x+1][y+1][z+1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x-1][y+1][z+1]	+inputImage[x+1][y-1][z-1]
											  -1.5f*(inputImage[x-1][y-1][z+1]	+inputImage[x-1][y-1][z-1]);
					direction = XmYmZ;						  
				} else if (d==1) {
					val1 =	 inputImage[x][y][z]			+inputImage[x+1][y+1][z-1]	+inputImage[x-1][y-1][z+1]
											  -1.5f*(inputImage[x-1][y-1][z-1]	+inputImage[x+1][y-1][z-1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x+1][y+1][z-1]	+inputImage[x-1][y-1][z+1]
											  -1.5f*(inputImage[x+1][y+1][z+1]	+inputImage[x-1][y+1][z+1]);
					direction = XpYmZ;						  
				} else 	if (d==2) {
					val1 =	 inputImage[x][y][z]			+inputImage[x+1][y+1][z+1]	+inputImage[x-1][y-1][z-1]
											  -1.5f*(inputImage[x-1][y+1][z+1]	+inputImage[x-1][y-1][z+1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x+1][y+1][z+1]	+inputImage[x-1][y-1][z-1]
											  -1.5f*(inputImage[x+1][y-1][z-1]	+inputImage[x+1][y+1][z-1]);
					direction = XpYpZ;						  
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==12) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 inputImage[x][y][z]			+inputImage[x-1][y+1][z+1]	+inputImage[x+1][y-1][z-1]
											  -1.5f*(inputImage[x+1][y-1][z+1]	+inputImage[x+1][y+1][z+1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x-1][y+1][z+1]	+inputImage[x+1][y-1][z-1]
											  -1.5f*(inputImage[x-1][y+1][z-1]	+inputImage[x-1][y-1][z-1]);
					direction = XmYmZ;						  
				} else if (d==1) {
					val1 =	 inputImage[x][y][z]			+inputImage[x+1][y-1][z+1]	+inputImage[x-1][y+1][z-1]
											  -1.5f*(inputImage[x-1][y-1][z-1]	+inputImage[x+1][y-1][z-1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x+1][y-1][z+1]	+inputImage[x-1][y+1][z-1]
											  -1.5f*(inputImage[x+1][y+1][z+1]	+inputImage[x-1][y+1][z+1]);
					direction = XmYpZ;						  
				} else if (d==2) {
					val1 =	 inputImage[x][y][z]			+inputImage[x-1][y-1][z-1]	+inputImage[x+1][y+1][z+1]
											  -1.5f*(inputImage[x-1][y+1][z+1]	+inputImage[x-1][y+1][z-1]);
					val2 =	 inputImage[x][y][z]			+inputImage[x-1][y-1][z-1]	+inputImage[x+1][y+1][z+1]
											  -1.5f*(inputImage[x-1][y+1][z+1]	+inputImage[x-1][y+1][z-1]);
					direction = XpYpZ;						  
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		}
		if (sign>0) return minval;
		else return 0.0f;
	}
	private final float[] directionVector(int d) {
		if (d==X) return new float[]{1.0f, 0.0f, 0.0f};
		else if (d==Y) return new float[]{0.0f, 1.0f, 0.0f};
		else if (d==Z) return new float[]{0.0f, 0.0f, 1.0f};
		else if (d==XpY) return new float[]{INVSQRT2, INVSQRT2, 0.0f};
		else if (d==YpZ) return new float[]{0.0f, INVSQRT2, INVSQRT2};
		else if (d==ZpX) return new float[]{INVSQRT2, 0.0f, INVSQRT2};
		else if (d==XmY) return new float[]{INVSQRT2, -INVSQRT2, 0.0f};
		else if (d==YmZ) return new float[]{0.0f, INVSQRT2, -INVSQRT2};
		else if (d==ZmX) return new float[]{-INVSQRT2, 0.0f, INVSQRT2};
		else if (d==XpYpZ) return new float[]{INVSQRT3, INVSQRT3, INVSQRT3};
		else if (d==XmYmZ) return new float[]{INVSQRT3, -INVSQRT3, -INVSQRT3};
		else if (d==XmYpZ) return new float[]{INVSQRT3, -INVSQRT3, INVSQRT3};
		else if (d==XpYmZ) return new float[]{INVSQRT3, INVSQRT3, -INVSQRT3};
		else return new float[]{0.0f, 0.0f, 0.0f};
	}
	private final byte[] directionNeighbor(int d) {
		if (d==X) return new byte[]{1, 0, 0};
		else if (d==Y) return new byte[]{0, 1, 0};
		else if (d==Z) return new byte[]{0, 0, 1};
		else if (d==XpY) return new byte[]{1, 1, 0};
		else if (d==YpZ) return new byte[]{0, 1, 1};
		else if (d==ZpX) return new byte[]{1, 0, 1};
		else if (d==XmY) return new byte[]{1, -1, 0};
		else if (d==YmZ) return new byte[]{0, 1, -1};
		else if (d==ZmX) return new byte[]{-1, 0, 1};
		else if (d==XpYpZ) return new byte[]{1, 1, 1};
		else if (d==XmYmZ) return new byte[]{1, -1, -1};
		else if (d==XmYpZ) return new byte[]{1, -1, 1};
		else if (d==XpYmZ) return new byte[]{1, 1, -1};
		else return new byte[]{0, 0, 0};
	}
	
	
	private final int neighborIndex(byte d, int id) {
		int idn=id;
		
			 if (d==X) 	idn+=1; 		
		else if (d==mX)	idn-=1;
		else if (d==Y) 	idn+=nx;
		else if (d==mY)	idn-=nx;
		else if (d==Z) 	idn+=nx*ny;
		else if (d==mZ)	idn-=nx*ny;
		else if (d==XpY) 	idn+=1+nx;
		else if (d==mXpY) 	idn-=1+nx;
		else if (d==YpZ) 	idn+=nx+nx*nx;
		else if (d==mYpZ)	idn-=nx+nx*ny;
		else if (d==ZpX) 	idn+=1+nx*ny;	
		else if (d==mZpX)	idn-=1+nx*ny;
		else if (d==XmY) 	idn+=1-nx;	
		else if (d==mXmY)	idn-=1-nx;
		else if (d==YmZ) 	idn+=nx-nx*ny;
		else if (d==mYmZ)	idn-=nx-nx*ny;
		else if (d==ZmX) 	idn+=1-nx*ny;
		else if (d==mZmX)	idn-=1-nx*ny;
		else if (d==XpYpZ) 		idn+=1+nx+nx*ny;
		else if (d==mXpYpZ)		idn-=1+nx+nx*ny;
		else if (d==XmYmZ) 		idn+=1-nx-nx*ny; 
		else if (d==mXmYmZ)		idn-=1-nx-nx*ny;
		else if (d==XmYpZ) 		idn+=1-nx+nx*ny;
		else if (d==mXmYpZ)		idn-=1-nx+nx*ny;
		else if (d==XpYmZ) 		idn+=1+nx-nx*ny; 
		else if (d==mXpYmZ)		idn-=1+nx-nx*ny;

		return idn;
	}
	private final float directionProduct(int dir, int id, float[][] imdir) {
		float dv=0.0f;
			
			 if (dir==X) dv = imdir[0][id];
		else if (dir==Y) dv = imdir[1][id];
		else if (dir==Z) dv = imdir[2][id];
		else if (dir==XpY) dv = imdir[0][id]/SQRT2+imdir[1][id]/SQRT2;
		else if (dir==YpZ) dv = imdir[1][id]/SQRT2+imdir[2][id]/SQRT2;
		else if (dir==ZpX) dv = imdir[2][id]/SQRT2+imdir[0][id]/SQRT2;
		else if (dir==XmY) dv = imdir[0][id]/SQRT2-imdir[1][id]/SQRT2;
		else if (dir==YmZ) dv = imdir[1][id]/SQRT2-imdir[2][id]/SQRT2;
		else if (dir==ZmX) dv = imdir[2][id]/SQRT2-imdir[0][id]/SQRT2;
		else if (dir==XpYpZ) dv = imdir[0][id]/SQRT3+imdir[1][id]/SQRT3+imdir[2][id]/SQRT3;
		else if (dir==XmYpZ) dv = imdir[0][id]/SQRT3-imdir[1][id]/SQRT3+imdir[2][id]/SQRT3;
		else if (dir==XpYmZ) dv = imdir[0][id]/SQRT3+imdir[1][id]/SQRT3-imdir[2][id]/SQRT3;
		else if (dir==XmYmZ) dv = imdir[0][id]/SQRT3-imdir[1][id]/SQRT3-imdir[2][id]/SQRT3;
		else if (dir==mX) dv = -imdir[0][id];
		else if (dir==mY) dv = -imdir[1][id];
		else if (dir==mZ) dv = -imdir[2][id];
		else if (dir==mXpY) dv = -(imdir[0][id]/SQRT2+imdir[1][id]/SQRT2);
		else if (dir==mYpZ) dv = -(imdir[1][id]/SQRT2+imdir[2][id]/SQRT2);
		else if (dir==mZpX) dv = -(imdir[2][id]/SQRT2+imdir[0][id]/SQRT2);
		else if (dir==mXmY) dv = -(imdir[0][id]/SQRT2-imdir[1][id]/SQRT2);
		else if (dir==mYmZ) dv = -(imdir[1][id]/SQRT2-imdir[2][id]/SQRT2);
		else if (dir==mZmX) dv = -(imdir[2][id]/SQRT2-imdir[0][id]/SQRT2);
		else if (dir==mXpYpZ) dv = -(imdir[0][id]/SQRT3+imdir[1][id]/SQRT3+imdir[2][id]/SQRT3);
		else if (dir==mXmYpZ) dv = -(imdir[0][id]/SQRT3-imdir[1][id]/SQRT3+imdir[2][id]/SQRT3);
		else if (dir==mXpYmZ) dv = -(imdir[0][id]/SQRT3+imdir[1][id]/SQRT3-imdir[2][id]/SQRT3);
		else if (dir==mXmYmZ) dv = -(imdir[0][id]/SQRT3-imdir[1][id]/SQRT3-imdir[2][id]/SQRT3);

		return dv;
	}	
	
	private final float[] spatialExpansion0D(float[] proba, int[] scale, float[] intens) {
		float[] propag = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			// expand along the same gaussian kernel as used in detection
			if (proba[xyz]>0.0f && scale[xyz]>0) {
				int dim = Numerics.ceil(Math.max(3.0f*scale[xyz]-0.5f,0.0f)); 
				for (int i=-dim;i<=dim;i++) for (int j=-dim;j<=dim;j++) for (int k=-dim;k<=dim;k++) {
					if (x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+k>=0 && z+k<nz) {
					    float ratio = (float)FastMath.exp( - 0.5f*(i*i+j*j+k*k)/(scale[xyz]*scale[xyz]) );
					    float imgratio = intens[xyz+i+j*nx+k*nx*ny]*intens[xyz];
						propag[xyz+i+j*nx+k*nx*ny] = Numerics.max(propag[xyz+i+j*nx+k*nx*ny], imgratio*ratio*proba[xyz]);
					}
				}
			}
			// by default, unless increased previously
			propag[xyz] = Numerics.max(propag[xyz],proba[xyz]);
		}
	    return propag;
	}

	private final float[] spatialUncertaintyRelaxation1D(float[] proba, byte[] dir, int ngbParam, float scale, float factor, int iterParam) {
		// run SUR with weighting against normal directions
		//mix with the neighbors?
		float[] certainty = new float[nx*ny*nz];   	
		float[] newproba = new float[nx*ny*nz];
		float[] prevproba = new float[nx*ny*nz];
		float[] ngbweight = new float[26];
		float[][] imgweight = new float[nx*ny*nz][26];  
		float[][] parallelweight = new float[26][26];
		for (int d1=0;d1<26;d1++) for (int d2=0;d2<26;d2++) {
			float[] dir1 = directionVector(d1);
			float[] dir2 = directionVector(d2);
			parallelweight[d1][d2] = (float)FastMath.pow(2.0f*FastMath.asin(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
		}
		
		// SOR-scheme?
		//float sorfactor = 1.95f;
		float sorfactor = 1.0f;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyzi = x+nx*y+nx*ny*z;
			for (byte dj=0;dj<26;dj++) {
				int xyzj = Ngb.neighborIndex(dj, xyzi, nx, ny, nz);
				imgweight[xyzi][dj] = sorfactor*diffusionWeightFunction(proba, dir, xyzi, xyzj, dj, parallelweight, scale)/ngbParam;
			}
		}
		for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) {
			newproba[xyzi] = proba[xyzi];
			prevproba[xyzi] = proba[xyzi];
		}
		
		float maxdiffParam = 1.0f;
		float meandiff = 1.0f;
		float t0diff = 1.0f;
		for (int t=0;t<iterParam && meandiff>0.001f*t0diff;t++) {
		//for (int t=0;t<iterParam && maxdiffParam>0.1f;t++) {
			BasicInfo.displayMessage("iterParam "+(t+1));
			maxdiffParam = 0.0f;
			meandiff = 0.0f;
			float ndiff = 0.0f;
			float nflip = 0.0f;
			float nproc = 0.0f;
							
			// mask: remove regions 0,1 away from ]0,1[ regions
			boolean[] diffmask = new boolean[nxyz];
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (prevproba[xyz]>0 && prevproba[xyz]<1) diffmask[xyz] = true;
				else diffmask[xyz] = false;
				for (int i=-1;i<=1 && !diffmask[xyz];i++) for (int j=-1;j<=1 && !diffmask[xyz];j++) for (int k=-1;k<=1 && !diffmask[xyz];k++) {
					if (prevproba[xyz+i+j*nx+k*nx*ny]>0 && prevproba[xyz+i+j*nx+k*nx*ny]<1) diffmask[xyz] = true;
				}
			}
			
			// main loop: label-per-label
			//BasicInfo.displayMessage("propagate gain for label "+n+"\n");
			//BasicInfo.displayMessage(".");

			// get the gain ; normalize
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				certainty[xyzi] = Numerics.abs(2.0f*prevproba[xyzi]-1.0f);
			}
			// propagate the values : diffusion
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				nproc++;
						
				float den = certainty[xyzi];
				float num = den*prevproba[xyzi];
				float prev = prevproba[xyzi];
					
				for (byte j=0;j<26;j++) {
					int xyzj = Ngb.neighborIndex(j, xyzi, nx, ny, nz);
					ngbweight[j] = certainty[xyzj]*imgweight[xyzi][j];
				}
				byte[] rank = Numerics.argmax(ngbweight, ngbParam);
				for (int l=0;l<ngbParam;l++) {
					int xyzl = Ngb.neighborIndex(rank[l], xyzi, nx, ny, nz);
					num += ngbweight[rank[l]]*prevproba[xyzl];
					den += ngbweight[rank[l]];
				}
				if (den>1e-9f) num /= den;
					
				newproba[xyzi] = num;
						
				meandiff += Numerics.abs(num-prev);
				ndiff++;
				maxdiffParam = Numerics.max(maxdiffParam, Numerics.abs(num-prev));
				if (prev<0.5f && num>0.5f) nflip++;
				if (prev>0.5f && num<0.5f) nflip++;
			}
			if (ndiff>0) meandiff /= ndiff;
			if (t==0) t0diff = meandiff;
			
			// make a hard copy
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				prevproba[xyzi] = newproba[xyzi];
			}
			BasicInfo.displayMessage("mean diff. "+meandiff+", max diff. "+maxdiffParam+"\n");
			//BasicInfo.displayMessage("n processed "+nproc+", n flipped "+nflip+"\n");
			
			//BasicInfo.displayMessage("n resorted"+nresort+"\n");
			
		}
		return newproba;
	}
	
	private final float[] spatialUncertaintyRelaxation2D(float[] proba, byte[] dir, int ngbParam, float scale, float factor, int iterParam) {
		// run SUR with weighting against normal directions
		//mix with the neighbors?
		float[] certainty = new float[nx*ny*nz];   	
		float[] newproba = new float[nx*ny*nz];
		float[] prevproba = new float[nx*ny*nz];
		float[] ngbweight = new float[26];
		float[][] imgweight = new float[nx*ny*nz][26];  
		float[][] orthogonalweight = new float[26][26];
		for (int d1=0;d1<26;d1++) for (int d2=0;d2<26;d2++) {
			float[] dir1 = directionVector(d1);
			float[] dir2 = directionVector(d2);
			orthogonalweight[d1][d2] = (float)FastMath.pow(2.0f*FastMath.acos(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
		}
		
		// SOR-scheme?
		//float sorfactor = 1.95f;
		float sorfactor = 1.0f;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyzi = x+nx*y+nx*ny*z;
			for (byte dj=0;dj<26;dj++) {
				int xyzj = Ngb.neighborIndex(dj, xyzi, nx, ny, nz);
				imgweight[xyzi][dj] = sorfactor*diffusionWeightFunction(proba, dir, xyzi, xyzj, dj, orthogonalweight, scale)/ngbParam;
			}
		}
		for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) {
			newproba[xyzi] = proba[xyzi];
			prevproba[xyzi] = proba[xyzi];
		}
		
		float maxdiffParam = 1.0f;
		float meandiff = 1.0f;
		float t0diff = 1.0f;
		for (int t=0;t<iterParam && meandiff>0.001f*t0diff;t++) {
		//for (int t=0;t<iterParam && maxdiffParam>0.1f;t++) {
			BasicInfo.displayMessage("iterParam "+(t+1));
			maxdiffParam = 0.0f;
			meandiff = 0.0f;
			float ndiff = 0.0f;
			float nflip = 0.0f;
			float nproc = 0.0f;
							
			// mask: remove regions 0,1 away from ]0,1[ regions
			boolean[] diffmask = new boolean[nxyz];
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (prevproba[xyz]>0 && prevproba[xyz]<1) diffmask[xyz] = true;
				else diffmask[xyz] = false;
				for (int i=-1;i<=1 && !diffmask[xyz];i++) for (int j=-1;j<=1 && !diffmask[xyz];j++) for (int k=-1;k<=1 && !diffmask[xyz];k++) {
					if (prevproba[xyz+i+j*nx+k*nx*ny]>0 && prevproba[xyz+i+j*nx+k*nx*ny]<1) diffmask[xyz] = true;
				}
			}
			
			// main loop: label-per-label
			//BasicInfo.displayMessage("propagate gain for label "+n+"\n");
			//BasicInfo.displayMessage(".");

			// get the gain ; normalize
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				certainty[xyzi] = Numerics.abs(2.0f*prevproba[xyzi]-1.0f);
			}
			// propagate the values : diffusion
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				nproc++;
						
				float den = certainty[xyzi];
				float num = den*prevproba[xyzi];
				float prev = prevproba[xyzi];
					
				for (byte j=0;j<26;j++) {
					int xyzj = Ngb.neighborIndex(j, xyzi, nx, ny, nz);
					ngbweight[j] = certainty[xyzj]*imgweight[xyzi][j];
				}
				byte[] rank = Numerics.argmax(ngbweight, ngbParam);
				for (int l=0;l<ngbParam;l++) {
					int xyzl = Ngb.neighborIndex(rank[l], xyzi, nx, ny, nz);
					num += ngbweight[rank[l]]*prevproba[xyzl];
					den += ngbweight[rank[l]];
				}
				if (den>1e-9f) num /= den;
					
				newproba[xyzi] = num;
						
				meandiff += Numerics.abs(num-prev);
				ndiff++;
				maxdiffParam = Numerics.max(maxdiffParam, Numerics.abs(num-prev));
				if (prev<0.5f && num>0.5f) nflip++;
				if (prev>0.5f && num<0.5f) nflip++;
			}
			if (ndiff>0) meandiff /= ndiff;
			if (t==0) t0diff = meandiff;
			
			// make a hard copy
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				prevproba[xyzi] = newproba[xyzi];
			}
			BasicInfo.displayMessage("mean diff. "+meandiff+", max diff. "+maxdiffParam+"\n");
			//BasicInfo.displayMessage("n processed "+nproc+", n flipped "+nflip+"\n");
			
			//BasicInfo.displayMessage("n resorted"+nresort+"\n");
			
		}
		return newproba;
	}
	
	private final float diffusionWeightFunction(float[] proba, byte[] dir, int xyz, int ngb, byte dngb, float[][] corr, float scale) {
		float diff = Numerics.abs((proba[xyz] - proba[ngb])/scale);
		float weight = 1.0f;
		if (dir[xyz]>-1 && dngb>-1) weight = corr[dir[xyz]][dngb];
		return weight/(1.0f+Numerics.square(diff/scale));
	}
	
	
	private final float[] probabilisticDiffusion1D(float[] proba, byte[] dir, int ngbParam, float maxdiffParam, float angle, float factor, int iterParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    	estimateSimpleDiffusionSimilarity1D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		// run the diffusion process
		float[] diffused = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
		}

		factor /= (float)ngbParam;
		for (int t=0;t<iterParam;t++) {
			BasicInfo.displayMessage("iterParamation "+(t+1)+": ");
			float diff = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				float prev = diffused[xyz];
				diffused[xyz] = proba[xyz];
				// weight with distance to 0 or 1
				for (int n=0;n<ngbParam;n++) {
					float ngb = diffused[neighborIndex(neighbor[n][xyz],xyz)];
					// remap neighbors?
					ngb = (float)FastMath.exp(ngb)-1.0f;
				
					// integration over the whole vessel (log version is more stable (??) )
					diffused[xyz] += factor*similarity[n][xyz]*ngb;
				}
				diffused[xyz] = (float)FastMath.log(1.0f + diffused[xyz]);
				
				if (diffused[xyz]+prev>0) {
					if (Numerics.abs(prev-diffused[xyz])/(prev+diffused[xyz])>diff) diff = Numerics.abs(prev-diffused[xyz])/(prev+diffused[xyz]);
				}
			}
			
			BasicInfo.displayMessage("diff "+diff+"\n");
			if (diff<maxdiffParam) t=iterParam;
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (iterParam<2) diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0), 0.0f, 1.0f);
			else diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0)/(1.0f+2.0f*factor), 0.0f, 1.0f);
		}
		return diffused;
	}

	private final float[] probabilisticFlooding1D(float[] proba, byte[] dir, int ngbParam, float maxdiffParam, float angle, float factor, int iterParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    		estimateSimpleDiffusionSimilarity1D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		// run the diffusion process
		float[] diffused = new float[nxyz];
		float[] previous = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			//diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
			diffused[xyz] = proba[xyz];
		}

		for (int t=0;t<iterParam;t++) {
			BasicInfo.displayMessage("iterParamation "+(t+1)+": ");
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				previous[xyz] = diffused[xyz];
			}
			float max = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				// weight with distance to 0 or 1
				for (int n=0;n<ngbParam;n++) {
					float ngb = previous[neighborIndex(neighbor[n][xyz],xyz)];
					// remap neighbors?
					//ngb = (float)FastMath.exp(ngb)-1.0f;
				
					// integration over the whole vessel (log version is more stable (??) )
					diffused[xyz] += factor*similarity[n][xyz]*ngb;
				}
				//diffused[xyz] = (float)FastMath.log(1.0f + diffused[xyz]);
				if (diffused[xyz]>max) max = diffused[xyz];
			}
			// renormalization (global)
			float diff = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				diffused[xyz] /= max;
				if (Numerics.abs(previous[xyz]-diffused[xyz])>diff) diff = Numerics.abs(previous[xyz]-diffused[xyz]);
			}
			
			BasicInfo.displayMessage("diff "+diff+"\n");
			if (diff<maxdiffParam) t=iterParam;
		}
		/*
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (iterParam==0) diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0), 0.0f, 1.0f);
			else diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0)/(1.0f+2.0f*factor), 0.0f, 1.0f);
		}
		*/
		return diffused;
	}

	private final float[] regionGrowing1D(float[] proba, byte[] dir, int ngbParam, float maxdiffParam, float angle, float factor, int iterParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    		estimateSimpleDiffusionSimilarity1D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		// sample all possible paths?
		float[] count = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			count[xyz] = 0.0f;
		}

		BinaryHeap2D heap = new BinaryHeap2D(nx+ny+nz, nx+ny+nz, BinaryHeap2D.MINTREE);
		BitSet used = new BitSet(nxyz);
		for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>factor) {
			BasicInfo.displayMessage(".");
			heap.reset();
			used.clear();
			heap.addValue(0, xyz, (byte)1);
			while (heap.isNotEmpty()) {
				int pos = heap.getFirstId();
				float val = heap.getFirst();
				heap.removeFirst();
				count[pos]++;
				used.set(pos);
				
				for (int n=0;n<ngbParam;n++) if (similarity[n][pos]>=0.5f ) {
					int ngb = neighborIndex(neighbor[n][pos],pos);
					if (!used.get(ngb)) {
						if (proba[ngb]>factor) {
							heap.addValue(val+1, ngb, (byte)1);
						}
					}
				}
			}
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			count[id] = (float)FastMath.log(1.0f+count[id]);
		}
		return count;
	}

    private final void estimateSimpleDiffusionSimilarity1D(byte[] dir, float[] proba, int ngbParam, byte[][] neighbor, float[][] similarity, float factor) {
    	
    	float[][] parallelweight = new float[26][26];
		for (int d1=0;d1<26;d1++) for (int d2=0;d2<26;d2++) {
			float[] dir1 = directionVector(d1);
			float[] dir2 = directionVector(d2);
			parallelweight[d1][d2] = (float)FastMath.pow(2.0f*FastMath.asin(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
		}
		float[] weight = new float[26];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int id = x+nx*y+nx*ny*z;
			if (proba[id]>0) {
				// find the N best aligned discrete directions
				for (byte d=0;d<NC2;d++) {
					int idn = neighborIndex(d,id);
					
					if (proba[idn]>0) {
						if (ngbParam==2) weight[d] = parallelweight[d][dir[id]];
						else weight[d] = parallelweight[d][dir[id]]*proba[idn];
					} else {
						weight[d] = 0.0f;
					}
				}
				byte[] ngb = Numerics.argmax(weight, ngbParam);
				for (int n=0;n<ngbParam;n++) {
					neighbor[n][id] = ngb[n];
					int idn = neighborIndex(ngb[n],id);
					if (proba[idn]>0) {
						similarity[n][id] = parallelweight[dir[id]][dir[idn]];
					} else {
						similarity[n][id] = 0.0f;
					}
				}
			}
		}
		
		return;
    }
    
	private final float[] probabilisticDiffusion2D(float[] proba, byte[] dir, int ngbParam, float maxdiffParam, float angle, float factor, int iterParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    	estimateSimpleDiffusionSimilarity2D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		// run the diffusion process
		float[] diffused = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
		}

		factor /= (float)ngbParam;
		for (int t=0;t<iterParam;t++) {
			BasicInfo.displayMessage("iterParamation "+(t+1)+": ");
			float diff = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				float prev = diffused[xyz];
				diffused[xyz] = proba[xyz];
				// weight with distance to 0 or 1
				for (int n=0;n<ngbParam;n++) {
					float ngb = diffused[neighborIndex(neighbor[n][xyz],xyz)];
					// remap neighbors?
					ngb = (float)FastMath.exp(ngb)-1.0f;
				
					// integration over the whole vessel (log version is more stable (??) )
					diffused[xyz] += factor*similarity[n][xyz]*ngb;
				}
				diffused[xyz] = (float)FastMath.log(1.0f + diffused[xyz]);
				
				if (Numerics.abs(prev-diffused[xyz])>diff) diff = Numerics.abs(prev-diffused[xyz]);
			}
			
			BasicInfo.displayMessage("diff "+diff+"\n");
			if (diff<maxdiffParam) t=iterParam;
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (iterParam<2) diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0), 0.0f, 1.0f);
			else diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0)/(1.0f+2.0f*factor), 0.0f, 1.0f);
		}
		return diffused;
	}

 	private final float[] probabilisticFlooding2D(float[] proba, byte[] dir, int ngbParam, float maxdiffParam, float angle, float factor, int iterParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    		estimateSimpleDiffusionSimilarity2D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		// run the diffusion process
		float[] diffused = new float[nxyz];
		float[] previous = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			//diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
			diffused[xyz] = proba[xyz];
		}

		for (int t=0;t<iterParam;t++) {
			BasicInfo.displayMessage("iterParamation "+(t+1)+": ");
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				previous[xyz] = diffused[xyz];
			}
			float max = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				// weight with distance to 0 or 1
				for (int n=0;n<ngbParam;n++) {
					float ngb = previous[neighborIndex(neighbor[n][xyz],xyz)];
					// remap neighbors?
					//ngb = (float)FastMath.exp(ngb)-1.0f;
				
					// integration over the whole vessel (log version is more stable (??) )
					diffused[xyz] += factor*similarity[n][xyz]*ngb;
				}
				//diffused[xyz] = (float)FastMath.log(1.0f + diffused[xyz]);
				if (diffused[xyz]>max) max = diffused[xyz];
			}
			// renormalization (global)
			float diff = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				diffused[xyz] /= max;
				if (Numerics.abs(previous[xyz]-diffused[xyz])>diff) diff = Numerics.abs(previous[xyz]-diffused[xyz]);
			}
			
			BasicInfo.displayMessage("diff "+diff+"\n");
			if (diff<maxdiffParam) t=iterParam;
		}
		/*
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (iterParam==0) diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0), 0.0f, 1.0f);
			else diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0)/(1.0f+2.0f*factor), 0.0f, 1.0f);
		}
		*/
		return diffused;
	}

	private final float[] regionGrowing2D(float[] proba, byte[] dir, int ngbParam, float maxdiffParam, float angle, float factor, int iterParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    		estimateSimpleDiffusionSimilarity2D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		// sample all possible paths?
		float[] count = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			count[xyz] = 0.0f;
		}

		BinaryHeap2D heap = new BinaryHeap2D(nx+ny+nz, nx+ny+nz, BinaryHeap2D.MINTREE);
		BitSet used = new BitSet(nxyz);
		for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>factor) {
			BasicInfo.displayMessage(".");
			heap.reset();
			used.clear();
			heap.addValue(0, xyz, (byte)1);
			while (heap.isNotEmpty()) {
				int pos = heap.getFirstId();
				float val = heap.getFirst();
				heap.removeFirst();
				count[pos]++;
				used.set(pos);
				
				for (int n=0;n<ngbParam;n++) if (similarity[n][pos]>=0.5f ) {
					int ngb = neighborIndex(neighbor[n][pos],pos);
					if (!used.get(ngb)) {
						if (proba[ngb]>factor) {
							heap.addValue(val+1, ngb, (byte)1);
						}
					}
				}
			}
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			count[id] = (float)FastMath.log(1.0f+count[id]);
		}
		return count;
	}
	private final float[] regionLabeling2D(float[] proba, byte[] dir, int ngbParam, float maxdiffParam, float angle, float factor, int iterParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    		estimateSimpleDiffusionSimilarity2D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		int Nlabel = 0;
		float[]   label = new float[nx*ny*nz];
		int[]       lb = new int[nx+ny+nz];
		int lbMin;
		int Nlb;
		int[]   connect = new int[26];
		int Nconnect;
		int AddLabel;
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++)
			label[xyz] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					int xyz = x + nx*y + nx*ny*z;
					if (proba[xyz]>factor) {
						// object point: neighbors ?
						Nconnect = 0;
						for (int n=0;n<ngbParam;n++) if (similarity[n][xyz]>=0.5f ) {
							int ngb = neighborIndex(neighbor[n][xyz],xyz);
							if (label[ngb] > 0) {
								connect[Nconnect] = lb[ (int)label[ngb] ];
								Nconnect++;
							}
						}
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[xyz] = lbMin;
						} else {
							// new, unconnected region
							label[xyz] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (int k=1;k<Nlabel;k++) {
			int c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// copy on label inputImage
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			label[xyz] = lb2[ lb[ (int)label[xyz] ] ];
		}
		// clean up
		lb = null;
		lb2 = null;
		connect = null;
	   
		return label;
	}

	private final float[] regionLabeling1D(float[] proba, byte[] dir, int ngbParam, float maxdiffParam, float angle, float factor, int iterParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    		estimateSimpleDiffusionSimilarity1D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		int Nlabel = 0;
		float[]   label = new float[nx*ny*nz];
		int[]       lb = new int[nx+ny+nz];
		int lbMin;
		int Nlb;
		int[]   connect = new int[26];
		int Nconnect;
		int AddLabel;
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++)
			label[xyz] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					int xyz = x + nx*y + nx*ny*z;
					if (proba[xyz]>factor) {
						// object point: neighbors ?
						Nconnect = 0;
						for (int n=0;n<ngbParam;n++) if (similarity[n][xyz]>=0.5f ) {
							int ngb = neighborIndex(neighbor[n][xyz],xyz);
							if (label[ngb] > 0) {
								connect[Nconnect] = lb[ (int)label[ngb] ];
								Nconnect++;
							}
						}
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[xyz] = lbMin;
						} else {
							// new, unconnected region
							label[xyz] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (int k=1;k<Nlabel;k++) {
			int c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// copy on label inputImage
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			label[xyz] = lb2[ lb[ (int)label[xyz] ] ];
		}
		// clean up
		lb = null;
		lb2 = null;
		connect = null;
	   
		return label;
	}

	private final void estimateSimpleDiffusionSimilarity2D(byte[] dir, float[] proba, int ngbParam, byte[][] neighbor, float[][] similarity, float factor) {
    	
    	float[][] parallelweight = new float[26][26];
		float[][] orthogonalweight = new float[26][26];
		for (int d1=0;d1<26;d1++) for (int d2=0;d2<26;d2++) {
			float[] dir1 = directionVector(d1);
			float[] dir2 = directionVector(d2);
			parallelweight[d1][d2] = (float)FastMath.pow(2.0f*FastMath.asin(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
			orthogonalweight[d1][d2] = (float)FastMath.pow(2.0f*FastMath.acos(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
		}
		float[] weight = new float[26];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int id = x+nx*y+nx*ny*z;
			if (proba[id]>0) {
				// find the N best planar discrete directions
				for (byte d=0;d<NC2;d++) {
					int idn = neighborIndex(d,id);
					
					if (proba[idn]>0) {
						if (ngbParam==2) weight[d] = orthogonalweight[d][dir[id]];
						else weight[d] = orthogonalweight[d][dir[id]]*proba[idn];
					} else {
						weight[d] = 0.0f;
					}
				}
				byte[] ngb = Numerics.argmax(weight, ngbParam);
				for (int n=0;n<ngbParam;n++) {
					neighbor[n][id] = ngb[n];
					int idn = neighborIndex(ngb[n],id);
					if (proba[idn]>0) {
						// similarity comes from the normal direction
						similarity[n][id] = parallelweight[dir[id]][dir[idn]];
					} else {
						similarity[n][id] = 0.0f;
					}
				}
			}
		}
		
		return;
    }
    
    private final float[] beliefPropagation1D(float[] proba, byte[] dir, int ngbParam, float angle, int iterParam, float maxdiffParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    		estimateSimpleDiffusionSimilarity1D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		float[][] messagefg = new float[ngbParam][nxyz];
		float[][] messagebg = new float[ngbParam][nxyz];
		
		float[][] newmsgfg = new float[ngbParam][nxyz];
		float[][] newmsgbg = new float[ngbParam][nxyz];
		
		// init
		for (int xyz=0;xyz<nxyz;xyz++) for (int n=0;n<ngbParam;n++) {
			messagefg[n][xyz] = 1.0f;
			messagebg[n][xyz] = 1.0f;
			newmsgfg[n][xyz] = 1.0f;
			newmsgbg[n][xyz] = 1.0f;
		}
		
		// min value for probas, similarities
		float minp = 0.00f;
		float mins = 0.00f;
		
		// message passing
		for (int t=0;t<iterParam;t++) {
			BasicInfo.displayMessage("iterParamation "+(t+1)+": ");
			float diff = 0.0f;
			//float prev;
			int ngb;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbParam;n++) {
					float simfgfg = similarity[n][xyz];
					float simbgfg = 0.0f; //(1.0f-similarity1[xyz]); //similarity1[xyz];
					float simfgbg = (1.0f-similarity[n][xyz]); // good results
					float simbgbg = 1.0f; //(1.0f-similarity1[xyz]); //1.0f; //similarity1[xyz];
				
					ngb = neighborIndex(neighbor[n][xyz],xyz);
					if (proba[ngb]>0) {
						// neighbor message, fg label
						float mfgfg = Numerics.max(mins,simfgfg)*Numerics.max(minp,proba[ngb]);
						for (int m=0;m<ngbParam;m++) {
							if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagefg[m][ngb]>0) mfgfg *= messagefg[m][ngb];
						}
						
						float mfgbg = Numerics.max(mins,simfgbg)*Numerics.max(minp,(1.0f-proba[ngb]));
						for (int m=0;m<ngbParam;m++) {
							if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagebg[m][ngb]>0) mfgbg *= messagebg[m][ngb];
						}
						newmsgfg[n][xyz] = Numerics.max(mfgfg, mfgbg);
						diff = Numerics.max(diff, Numerics.abs(newmsgfg[n][xyz]-messagefg[n][xyz]));
					
						// neighbor message, bg label
						float mbgbg = Numerics.max(mins,simbgbg)*Numerics.max(minp,(1.0f-proba[ngb]));
						for (int m=0;m<ngbParam;m++) {
							if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagebg[m][ngb]>0) mbgbg *= messagebg[m][ngb];
						}
					
						float mbgfg = Numerics.max(mins,simbgfg)*Numerics.max(minp,proba[ngb]);
						for (int m=0;m<ngbParam;m++) {
							if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagefg[m][ngb]>0) mbgfg *= messagefg[m][ngb];
						}
						newmsgbg[n][xyz] = Numerics.max(mbgbg, mbgfg);
						diff = Numerics.max(diff, Numerics.abs(newmsgbg[n][xyz]-messagebg[n][xyz]));
					}
				}
			}
			// copy new messages onto older ones
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbParam;n++) {
					messagefg[n][xyz] = newmsgfg[n][xyz];
					messagebg[n][xyz] = newmsgbg[n][xyz];
				}
			}
			BasicInfo.displayMessage("diff "+diff+"\n");
			if (diff < maxdiffParam) t = iterParam;
		}
		
		// compute final belief
		float[] belief = new float[nxyz];
		float 	bgbelief;
		
		for (int xyz=0;xyz<nxyz;xyz++) {
			belief[xyz] = proba[xyz];
			
			if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbParam;n++) {
					if (messagefg[n][xyz]>0) belief[xyz] *= messagefg[n][xyz];
				}
			}
			
			bgbelief = (1.0f-proba[xyz]);
			if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbParam;n++) {
					if (messagebg[n][xyz]>0) bgbelief *= messagebg[n][xyz];
				}
			}
			// normalize
			//belief[xyz] = belief[xyz]/Numerics.max(1e-9f,belief[xyz]+bgbelief);
			belief[xyz] = belief[xyz]/(belief[xyz]+bgbelief);
		}
		return belief;
	}
	
    private final float[] beliefPropagation2D(float[] proba, byte[] dir, int ngbParam, float angle, int iterParam, float maxdiffParam) {
		// mask out inputImage boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (x<1 || y<1 || z<1 || x>=nx-1 || y>=ny-1 || z>=nz-1) proba[xyz] = 0.0f;
		}
		
		// build a similarity function from aligned neighbors
		float[][] similarity = new float[ngbParam][nxyz];
		byte[][] neighbor = new byte[ngbParam][nxyz];
    		estimateSimpleDiffusionSimilarity2D(dir, proba, ngbParam, neighbor, similarity, angle);
		
		float[][] messagefg = new float[ngbParam][nxyz];
		float[][] messagebg = new float[ngbParam][nxyz];
		
		float[][] newmsgfg = new float[ngbParam][nxyz];
		float[][] newmsgbg = new float[ngbParam][nxyz];
		
		// init
		for (int xyz=0;xyz<nxyz;xyz++) for (int n=0;n<ngbParam;n++) {
			messagefg[n][xyz] = 1.0f;
			messagebg[n][xyz] = 1.0f;
			newmsgfg[n][xyz] = 1.0f;
			newmsgbg[n][xyz] = 1.0f;
		}
		
		// min value for probas, similarities
		float minp = 0.00f;
		float mins = 0.00f;
		
		// message passing
		for (int t=0;t<iterParam;t++) {
			BasicInfo.displayMessage("iterParamation "+(t+1)+": ");
			float diff = 0.0f;
			//float prev;
			int ngb;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbParam;n++) {
					float simfgfg = similarity[n][xyz];
					float simbgfg = 0.0f; //(1.0f-similarity1[xyz]); //similarity1[xyz];
					float simfgbg = (1.0f-similarity[n][xyz]); // good results
					float simbgbg = 1.0f; //(1.0f-similarity1[xyz]); //1.0f; //similarity1[xyz];
				
					ngb = neighborIndex(neighbor[n][xyz],xyz);
					
					if (proba[ngb]>0) {
						// neighbor message, fg label
						float mfgfg = Numerics.max(mins,simfgfg)*Numerics.max(minp,proba[ngb]);
						for (int m=0;m<ngbParam;m++) {
							if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagefg[m][ngb]>0) mfgfg *= messagefg[m][ngb];
						}
						
						float mfgbg = Numerics.max(mins,simfgbg)*Numerics.max(minp,(1.0f-proba[ngb]));
						for (int m=0;m<ngbParam;m++) {
							if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagebg[m][ngb]>0) mfgbg *= messagebg[m][ngb];
						}
						newmsgfg[n][xyz] = Numerics.max(mfgfg, mfgbg);
						diff = Numerics.max(diff, Numerics.abs(newmsgfg[n][xyz]-messagefg[n][xyz]));
					
						// neighbor message, bg label
						float mbgbg = Numerics.max(mins,simbgbg)*Numerics.max(minp,(1.0f-proba[ngb]));
						for (int m=0;m<ngbParam;m++) {
							if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagebg[m][ngb]>0) mbgbg *= messagebg[m][ngb];
						}
					
						float mbgfg = Numerics.max(mins,simbgfg)*Numerics.max(minp,proba[ngb]);
						for (int m=0;m<ngbParam;m++) {
							if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagefg[m][ngb]>0) mbgfg *= messagefg[m][ngb];
						}
						newmsgbg[n][xyz] = Numerics.max(mbgbg, mbgfg);
						diff = Numerics.max(diff, Numerics.abs(newmsgbg[n][xyz]-messagebg[n][xyz]));
					}
				}
			}
			// copy new messages onto older ones
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbParam;n++) {
					messagefg[n][xyz] = newmsgfg[n][xyz];
					messagebg[n][xyz] = newmsgbg[n][xyz];
				}
			}
			BasicInfo.displayMessage("diff "+diff+"\n");
			if (diff < maxdiffParam) t = iterParam;
		}
		
		// compute final belief
		float[] belief = new float[nxyz];
		float 	bgbelief;
		
		for (int xyz=0;xyz<nxyz;xyz++) {
			belief[xyz] = proba[xyz];
			
			if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbParam;n++) {
					if (messagefg[n][xyz]>0) belief[xyz] *= messagefg[n][xyz];
				}
			}
			
			bgbelief = (1.0f-proba[xyz]);
			if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbParam;n++) {
					if (messagebg[n][xyz]>0) bgbelief *= messagebg[n][xyz];
				}
			}
			// normalize
			//belief[xyz] = belief[xyz]/Numerics.max(1e-9f,belief[xyz]+bgbelief);
			belief[xyz] = belief[xyz]/(belief[xyz]+bgbelief);
		}
		return belief;
	}
	

}
