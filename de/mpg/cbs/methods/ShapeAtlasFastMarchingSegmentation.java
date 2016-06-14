package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.jama.*;
import gov.nih.mipav.model.file.FileInfoBase;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm performs a simple alignment of the image with a topology template
 *	<p>
 *	The algorithm handles all the little things needed for image cropping, conversion
 * 	to an image buffer, padding, etc.
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class ShapeAtlasFastMarchingSegmentation {
		
	// data buffers
	private 	float[][]		image;  			// original images
	private 	int				nix,niy,niz;   		// image dimensions
	private 	float			rix,riy,riz;   		// image resolutions
	private		String[]		modality;			// image modality / contrast
	private		float[]			imrange;			// image intensity range (robust estimate)
	private		int				nc;					// number of channels
	
	// atlas parameters
	private		SimpleShapeAtlas	atlas;
	private 	int 				nobj;    	// number of shapes
	private 	int				nax,nay,naz;   			// atlas dimensions
	private 	float			rax,ray,raz;   			// atlas resolutions
	private		float[][]		intensity;		// atlas intensities

	// segmentation parameters
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	private 	BinaryHeap2D		heap;		   		// the binary tree used for the fast marching
	private		float[]				score;
	private		byte[]				segmentation;
	private		boolean[]			mask;
	private		byte[]				processed;
	
	// output results
	private		int				nsx,nsy,nsz;			// shape dimensions
	private 	float			rsx,rsy,rsz;   			// shape resolutions
	private		float			scaling;				// shape scaling
	private		byte			napprox;				// max number of retained functions
	private		float			pvsize;
	
	//private		float[][]		min, max;
	private		float[][]		gain;
	private		float[][]		bestgain;
	private		byte[][]		bestlabel;
	private		float[][]		field;
	private		float[]			curvweight;
	private		float			sigmaI, sigmaS;
	private		float[] 		sigmaN;
	private		float			factor;
	private		float			threshold1;
	private		float			threshold2;
	private		float			slope;
	private		int				ngbsize;
	private		float			alpha;
	
    static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// constants
	private static final	float   ISQRT2 = (float)(1.0/Math.sqrt(2.0f));
	private static final	float   ZERO = 1E-20f;
	private static final	float   INF = 1E20f;
		
	private	static	final	byte	EMPTY = -1;
			
	public	static	final	byte	pX = 0;
	public	static	final	byte	pY = 1;
	public	static	final	byte	pZ = 2;
	public	static	final	byte	mX = 3;
	public	static	final	byte	mY = 4;
	public	static	final	byte	mZ = 5;
	public	static 	final 	byte 	pXpY = 6;
	public	static 	final 	byte 	pYpZ = 7;
	public	static 	final 	byte 	pZpX = 8;
	public	static 	final 	byte 	pXmY = 9;
	public	static 	final 	byte 	pYmZ = 10;
	public	static 	final 	byte 	pZmX = 11;
	public	static 	final 	byte 	mXpY = 12;
	public	static 	final 	byte 	mYpZ = 13;
	public	static 	final 	byte 	mZpX = 14;
	public	static 	final 	byte 	mXmY = 15;
	public	static 	final 	byte 	mYmZ = 16;
	public	static 	final 	byte 	mZmX = 17;
	public	static 	final 	byte 	pXpYpZ = 18;
	public	static 	final 	byte 	pXmYmZ = 19;
	public	static 	final 	byte 	pXmYpZ = 20;
	public	static 	final 	byte 	pXpYmZ = 21;
	public	static 	final 	byte 	mXpYpZ = 22;
	public	static 	final 	byte 	mXmYmZ = 23;
	public	static 	final 	byte 	mXmYpZ = 24;
	public	static 	final 	byte 	mXpYmZ = 25;
	public	static 	final 	byte 	NGB = 26;

	private static final byte[] ngbx = {+1,  0,  0, -1,  0,  0, +1,  0, +1, +1,  0, -1, -1,  0, +1, -1,  0, -1, +1, +1, +1, +1, -1, -1, -1, -1};
	private static final byte[] ngby = { 0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1,  0, +1, -1, -1, +1, +1, -1, -1, +1};
	private static final byte[] ngbz = { 0,  0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1, +1, -1, +1, -1, +1, -1, +1, -1};
	
	private static final byte[] ngpx = {+1,  0,  0, +1,  0, +1, +1};
	private static final byte[] ngpy = { 0, +1,  0, +1, +1,  0, +1};
	private static final byte[] ngpz = { 0,  0, +1,  0, +1, +1, +1};
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	public	static 	final 	byte 	X = 0;
	public	static 	final 	byte 	Y = 1;
	public	static 	final 	byte 	Z = 2;

	/**
	 *  constructor
	 */
	public ShapeAtlasFastMarchingSegmentation(float[][] img_, String[] mod_, float[] rng_, int nc_,
										int nix_, int niy_, int niz_,
										float rix_, float riy_, float riz_,
										SimpleShapeAtlas atlas_,
										float scaling_,
										float p1_, float p2_, float p3_,
										float p4_, float p5_, float p6_,
										int n1_, String connectivityType_, float alpha_,
										int na_, float pv_) {
		image = img_;
		imrange = rng_;
		modality = mod_;
		nc = nc_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		
		scaling = scaling_;
		nsx = Numerics.ceil(scaling*nix);
		nsy = Numerics.ceil(scaling*niy);
		nsz = Numerics.ceil(scaling*niz);
		
		if (debug) System.out.println("scaled dimensions: "+nsx+" x "+nsy+" x "+nsz);
		
		
		rsx = rix/scaling;
		rsy = riy/scaling;
		rsz = riz/scaling;
		
		atlas = atlas_;
		nobj = atlas.getNumber();
		nax = atlas.getShapeDim()[0];
		nay = atlas.getShapeDim()[1];
		naz = atlas.getShapeDim()[2];
		rax = atlas.getShapeRes()[0];
		ray = atlas.getShapeRes()[1];
		raz = atlas.getShapeRes()[2];		
		intensity = atlas.getIntensityPriors(modality, nc);
		for (int c=0;c<nc;c++) {
			BasicInfo.displayMessage("channel "+c+" intensity priors: "+atlas.displayVector(intensity[c])+"\n");
		}
		
		sigmaI = p1_;
		sigmaN = new float[nc];
		for (int c=0;c<nc;c++) sigmaN[c] = p2_;
		sigmaS = p3_;
		
		threshold1 = p4_;
		threshold2 = p5_;
		
		slope	= p6_;
		
		ngbsize = n1_;
		
		alpha = alpha_;
		
		napprox = (byte)na_;
		
		pvsize = pv_;
		
		// 6-neighborhood: pre-compute the index offsets
		//xoff = new int[]{1, -1, 0, 0, 0, 0};
		//yoff = new int[]{0, 0, nsx, -nsx, 0, 0};
		//zoff = new int[]{0, 0, 0, 0, nsx*nsy, -nsx*nsy};

		try {
			segmentation = new byte[nsx*nsy*nsz];	
			//score = new float[nsx*nsy*nsz];	
			mask = new boolean[nsx*nsy*nsz];	
			// topology luts
			/*
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6")) lut = new CriticalPointLUT("critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26")) lut = new CriticalPointLUT("critical626LUT.raw.gz",200);
			else if (connectivityType_.equals("18/6")) lut = new CriticalPointLUT("critical186LUT.raw.gz",200);
			else if (connectivityType_.equals("6/18")) lut = new CriticalPointLUT("critical618LUT.raw.gz",200);
			else if (connectivityType_.equals("6/6")) lut = new CriticalPointLUT("critical66LUT.raw.gz",200);
			else if (connectivityType_.equals("wcs")) {
				lut = new CriticalPointLUT("criticalWCLUT.raw.gz",200);
				checkComposed=false;
			}
			else if (connectivityType_.equals("wco")) {
				lut = new CriticalPointLUT("critical66LUT.raw.gz",200);
				checkComposed=true;
			}
			else if (connectivityType_.equals("no")) {
				lut = null;
				checkTopology=false;
			}
			else {
				lut = null;
				checkTopology=false;
			}
			if (checkTopology) {
				if (!lut.loadCompressedPattern()) {
					finalize();
					System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
					BasicInfo.displayMessage("Problem loading the algorithm's LUT from: "+lut.getFilename()+"\n");
				} else {
					System.out.println("LUT loaded from: "+lut.getFilename());
				}
			}*/
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nsx; x++) for (int y=0; y<nsy; y++) for (int z = 0; z<nsz; z++) {
			if (x>1 && x<nsx-2 && y>1 && y<nsy-2 && z>1 && z<nsz-2) mask[x+nsx*y+nsx*nsy*z] = true;
			else mask[x+nsx*y+nsx*nsy*z] = false;
		}

		// init segmentation from atlas
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyz = x + nsx*y + nsx*nsy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x/scaling, y/scaling, z/scaling);
			
			int init = ImageInterpolation.nearestNeighborInterpolation(atlas.getTemplate(),(byte)1,XA[0],XA[1],XA[2],nax,nay,naz);
			byte nlb = EMPTY;
			for (byte n=0; n<nobj; n++) {
				if (atlas.getLabels()[n]==init) {
					nlb = n;
					continue;
				}
			}
			segmentation[xyz] = nlb;		
		}

		return;
	}

	final public void finalize() {
		gain = null;
		//min = null;
		//max = null;
		segmentation = null;
		System.gc();
	}

	public final float[][] getGain() { return gain; }

	public final float[][] getBestGain() { return bestgain; }

	public final byte[][] getBestLabel() { return bestlabel; }

	public final float[][] getField() { return field; }

	public final float[] getCurvatureWeight() { return curvweight; }

	public final byte[] getSegmentation() { return segmentation; }
	
	public final int nsx() { return nsx; }
	public final int nsy() { return nsy; }
	public final int nsz() { return nsz; }
	
	public final float rsx() { return rsx; }
	public final float rsy() { return rsy; }
	public final float rsz() { return rsz; }
	
	/** 
	 *  compute the average and std of image intensity at each atlas voxel
	 */
	/* 
    final public void computeImageMinMaxValues() {
    	min = new float[nc][nsx*nsy*nsz];
    	max = new float[nc][nsx*nsy*nsz];
    	
     	for (int c=0;c<nc;c++) for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
    		min[c][xyz] = INF;
    		max[c][xyz] = -INF;
    	}
    	
		float[] Xs = new float[3];
		int xyzs, xyzi;
		// scans the entire image only once (note: some atlas locations may not have values)
       for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
       //for (int x=ngbsize;x<nix-ngbsize;x++) for (int y=ngbsize;y<niy-ngbsize;y++) for (int z=ngbsize;z<niz-ngbsize;z++) {
	        xyzi = x + nix*y + nix*niy*z;
        	xyzs = Numerics.floor(x*scaling) + nsx*Numerics.floor(y*scaling) + nsx*nsy*Numerics.floor(z*scaling);
        	
        	for (int c=0;c<nc;c++) {
				min[c][xyzs] = Numerics.min(min[c][xyzs], image[c][xyzi]);
				max[c][xyzs] = Numerics.max(max[c][xyzs], image[c][xyzi]);
				// use a region
				//for (int i=-ngbsize;i<=ngbsize;i++) for (int j=-ngbsize;j<=ngbsize;j++) for (int l=-ngbsize;l<=ngbsize;l++) {
				for (int n=0;n<ngbsize;n++) { 
					min[c][xyzs] = Numerics.min(min[c][xyzs], image[c][xyzi+ngbx[n]+nix*ngby[n]+nix*niy*ngbz[n]]);
					max[c][xyzs] = Numerics.max(max[c][xyzs], image[c][xyzi+ngbx[n]+nix*ngby[n]+nix*niy*ngbz[n]]);
				}
			}
        }
		for (int c=0;c<nc;c++) {
			// normalize everything from previous robust estimate
			for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
				min[c][xyz] = min[c][xyz]/imrange[c];
				max[c][xyz] = max[c][xyz]/imrange[c];
			}
		}
        return;
    }// image min max
    */
    
	/**
	 *	compute the intensity priors given the atlas
	 */
	/* 
    public final void computeAtlasGainFunction() {
    	gain = new float[nobj][nsx*nsy*nsz];
    	
    	// compute the local priors
    	float shape;
    	float[] intens = new float[nc];
    	float val, diff, mindiff, maxdiff;
    	float shapesum, shapemax;
    	float img;
    	int xyzs, xyzi, xyzb;
    	int xi, yi, zi;
    	float[] XA = new float[3];
    	for (int x=1;x<nsx-1;x++) for (int y=1;y<nsy-1;y++) for (int z=1;z<nsz-1;z++) {
    		xyzs = x + nsx*y + nsx*nsy*z;
        	atlas.imageToShapeCoordinates(XA, x/scaling, y/scaling, z/scaling);
        	
        	for (int n=0;n<nobj;n++) {
        		for (int c=0;c<nc;c++) {
					if (modality[c].equals("T1map")) {
						intens[c] = -1.0f;
						for (int t=0;t<atlas.getT1map()[n].length;t+=3) {
							xi = Numerics.round(x/scaling); 
							yi = Numerics.round(y/scaling);
							zi = Numerics.round(z/scaling);
							
							xyzi = xi + nix*yi + nix*niy*zi;
							
							diff = Numerics.abs(image[c][xyzi]/imrange[c]-atlas.getT1map()[n][t])/atlas.getT1map()[n][t+1];
							for (int s=1;s<1.0f/scaling;s++) {
								for (int m=0;m<NGB;m++) {
									if (  (Numerics.abs( xi+s*ngbx[m] - x/scaling ) <= 0.5f/scaling)
										&&(Numerics.abs( yi+s*ngby[m] - y/scaling ) <= 0.5f/scaling)
										&&(Numerics.abs( zi+s*ngbz[m] - z/scaling ) <= 0.5f/scaling) ) {
									
										xyzb = xyzi + s*ngbx[m] + nix*s*ngby[m] + nix*niy*s*ngbz[m];
										diff = Numerics.min(diff, Numerics.abs(image[c][xyzb]/imrange[c]-atlas.getT1map()[n][t])/atlas.getT1map()[n][t+1]);
									}
								}
							}
							intens[c] = Numerics.max(intens[c], atlas.getT1map()[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
						}
					} else if (modality[c].equals("FlatImg")) {
						intens[c] = -1.0f;
						for (int t=0;t<atlas.getFlatImg()[n].length;t+=3) {
							xi = Numerics.round(x/scaling); 
							yi = Numerics.round(y/scaling);
							zi = Numerics.round(z/scaling);
							
							xyzi = xi + nix*yi + nix*niy*zi;
							
							diff = Numerics.abs(image[c][xyzi]/imrange[c]-atlas.getFlatImg()[n][t])/atlas.getFlatImg()[n][t+1];
							for (int s=1;s<1.0f/scaling;s++) {
								for (int m=0;m<NGB;m++) {
									if (  (Numerics.abs( xi+s*ngbx[m] - x/scaling ) <= 0.5f/scaling)
										&&(Numerics.abs( yi+s*ngby[m] - y/scaling ) <= 0.5f/scaling)
										&&(Numerics.abs( zi+s*ngbz[m] - z/scaling ) <= 0.5f/scaling) ) {
									
										xyzb = xyzi + s*ngbx[m] + nix*s*ngby[m] + nix*niy*s*ngbz[m];
										diff = Numerics.min(diff, Numerics.abs(image[c][xyzb]/imrange[c]-atlas.getFlatImg()[n][t])/atlas.getFlatImg()[n][t+1]);
									}
								}
							}
							intens[c] = Numerics.max(intens[c], atlas.getFlatImg()[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
						}
					} else {
						System.out.print("!");
					}
				}
				
				if (n==0) val = ImageInterpolation.linearInterpolation(atlas.getShape(n),1.0f,XA[0],XA[1],XA[2],nax,nay,naz);
				else val = ImageInterpolation.linearInterpolation(atlas.getShape(n),0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
				shape = (val - sigmaS)/(val + sigmaS - 2.0f*sigmaS*val);
				/*
				// use semi-product, or simple sums?? both!!
				float sum = 0.0f;
				for (int c=0;c<nc;c++) sum += intens[c];
				
				if (shape>0 && sum>0) {
					gain[n][xyzs] = shape*sum;
				} else if (shape<0 && sum<0) {
					gain[n][xyzs] = -shape*sum;
				} else {
					gain[n][xyzs] = 0;	
				}
				//gain[n][xyzs] = (shape + sum)/(nc+1);
				*/
				/*
				// full semi-product (most restrictive)
				float sum = intens[0];
				for (int c=1;c<nc;c++) {
						if (intens[c]>0 && sum>0) {
						sum = intens[c]*sum;
					} else if (intens[c]<0 && sum<0) {
						sum = -intens[c]*sum;
					} else {
						sum = 0;	
					}
				}
				if (shape>0 && sum>0) {
					gain[n][xyzs] = shape*sum;
				} else if (shape<0 && sum<0) {
					gain[n][xyzs] = -shape*sum;
				} else {
					gain[n][xyzs] = 0;	
				}
				*/
				/*
				// use gain composition
				if (nc==1) {
					float w1 = Numerics.bounded((1+shape)/(2+shape+intens[0]),0,1);
					float w2 = Numerics.bounded((1+intens[0])/(2+shape+intens[0]),0,1);
					gain[n][xyzs] = w1*intens[0] + w2*shape;
				} else if (nc==2) {
					float w1 = Numerics.bounded((1+intens[1])*(1+shape)/((1+intens[0])*(1+intens[1]) + (1+intens[0])*(1+shape) + (1+intens[1])*(1+shape)),0,1);
					float w2 = Numerics.bounded((1+intens[0])*(1+shape)/((1+intens[0])*(1+intens[1]) + (1+intens[0])*(1+shape) + (1+intens[1])*(1+shape)),0,1);
					float w3 = Numerics.bounded((1+intens[0])*(1+intens[1])/((1+intens[0])*(1+intens[1]) + (1+intens[0])*(1+shape) + (1+intens[1])*(1+shape)),0,1);
					gain[n][xyzs] = w1*intens[0] + w2*intens[1] + w3*shape;
				} else {
					gain[n][xyzs] = shape;	
				}
				
			}
		}  
		// quantify the results
		float[] positive = new float[nobj];
		float[] negative = new float[nobj];
		float[]	zero = new float[nobj];
        for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
        	for (int n=0;n<nobj;n++) {
        		if (gain[n][xyz]>0) positive[n] += gain[n][xyz];
        		else if (gain[n][xyz]<0) negative[n] += gain[n][xyz];
        		else zero[n]++;
        	}
        }
        System.out.println("positive gain: "+atlas.displayVector(positive));
        System.out.println("negative gain: "+atlas.displayVector(negative));
        System.out.println(" zero    gain: "+atlas.displayVector(zero));
    }
    */
    public final void normalizeGainFunctions() {
        for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
        	float sum = 0.0f;
        	/*
        	float pos = 0.0f;
        	float neg = 0.0f;
        	float sump = 0.0f;
        	float proba;
        	*/
        	for (int n=0;n<nobj;n++) {
        		sum += gain[n][xyz];
        		/*
        		sump += p0*(1.0f+gain[n][xyz])/(1.0f-gain[n][xyz]+2.0f*p0*gain[n][xyz]);
        		if (gain[n][xyz]>0) pos += gain[n][xyz];
        		else if (gain[n][xyz]<0) neg += gain[n][xyz];
        		*/
        	}
        	for (int n=0;n<nobj;n++) {
        		// sum over everything ? assumes only values with p>1/2 are used...
        		gain[n][xyz] = 2.0f*(1.0f+gain[n][xyz])/(nobj+sum)-1.0f;
        		/*
        		// with a different threshold, not as nice!!
        		proba = p0*(1.0f+gain[n][xyz])/(1.0f-gain[n][xyz]+2.0f*p0*gain[n][xyz]);
        		gain[n][xyz] = (proba/sump-p0)/(proba/sump+p0-2.0f*proba/sump*p0);
        		// sum over same sign ?
         		if (gain[n][xyz]>0) gain[n][xyz] = gain[n][xyz]/pos;
        		else if (gain[n][xyz]<0) gain[n][xyz] = -gain[n][xyz]/neg;
        		*/
        	}
        }    	
    	return;	
    }
    
	/**
	 *	compute the intensity priors given the atlas
	 */
    public final void computeBestGainFunction() {
    	//gain = new float[nobj][nsx*nsy*nsz];
    	bestgain = new float[napprox][nsx*nsy*nsz];
    	bestlabel = new byte[napprox][nsx*nsy*nsz];
    	
    	byte flat=-1;
    	byte t1map = -1;
    	byte t1ridges=-1;
    	for (byte n=0;n<nc;n++) {
			if (modality[n].equals("FlatImg")) flat=n;
			else if (modality[n].equals("T1map")) t1map=n;
			else if (modality[n].equals("T1ridges")) t1ridges=n;
		}

    	
    	// compute the local priors
    	float shape;
    	float[] intens = new float[nobj];
    	float[] contrast = new float[3];
    	float val, diff, mindiff, maxdiff, proba;
    	float shapesum, shapemax;
    	float img;
    	int xyzs, xyzi, xyzb;
    	int xi, yi, zi;
    	float[] XA = new float[3];
    	for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
    		xyzs = x + nsx*y + nsx*nsy*z;
        	atlas.imageToShapeCoordinates(XA, x/scaling, y/scaling, z/scaling);

			xi = Numerics.round(x/scaling); 
			yi = Numerics.round(y/scaling);
			zi = Numerics.round(z/scaling);
			
			xyzi = xi + nix*yi + nix*niy*zi;

			for (int n=0;n<nobj;n++) {
				intens[n] = -1.0f;
			}
			
			// look around the entire neighborhood, but with coupled values
			for (int xid = Numerics.floor(xi-pvsize);xid<xi+pvsize+1;xid++) if (xid>=0 && xid<nix) {
				for (int yid = Numerics.floor(yi-pvsize);yid<yi+pvsize+1;yid++) if (yid>=0 && yid<niy) {
					for (int zid = Numerics.floor(zi-pvsize);zid<zi+pvsize+1;zid++) if (zid>=0 && zid<niz) {
						if ( (xid-xi)*(xid-xi) + (yid-yi)*(yid-yi) + (zid-zi)*(zid-zi) <= pvsize*pvsize) {
							//processing
							int xyzid = xid + nix*yid + nix*niy*zid;
							
							for (int n=0;n<nobj;n++) {
								// compute the exisiting contrasts
								if (flat>-1) {
									contrast[flat] = -1.0f;
									for (int t=0;t<atlas.getMap(atlas.ISO)[n].length;t+=3) {
										diff = Numerics.abs(image[flat][xyzid]/imrange[flat]-atlas.getMap(atlas.ISO)[n][t])/atlas.getMap(atlas.ISO)[n][t+1];
										contrast[flat] = Numerics.max(contrast[flat], atlas.getMap(atlas.ISO)[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
									}
								}
								if (t1map>-1) {
									contrast[t1map] = -1.0f;
									for (int t=0;t<atlas.getMap(atlas.T1M)[n].length;t+=3) {
										diff = Numerics.abs(image[t1map][xyzid]/imrange[t1map]-atlas.getMap(atlas.T1M)[n][t])/atlas.getMap(atlas.T1M)[n][t+1];
										contrast[t1map] = Numerics.max(contrast[t1map], atlas.getMap(atlas.T1M)[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
									}
								}
								if (t1ridges>-1) {
									proba = Numerics.min(1.0f,Numerics.abs(image[t1ridges][xyzid]/imrange[t1ridges]));
									contrast[t1ridges] = Numerics.max(-1.0f, (proba-atlas.getMap(atlas.sPV)[n][1])
																		/(proba + atlas.getMap(atlas.sPV)[n][1]*(1.0f-2.0f*proba) ) );
								}
								/*
								if (nc==1) {
									intens[n] = contrast[0];
								} else if (nc==2) {
									if (flat>-1 && t1map>-1) {
										float wf = Numerics.bounded((1+contrast[t1map])/(2+contrast[flat]+contrast[t1map]),0,1);
										float wt = Numerics.bounded((1+contrast[flat])/(2+contrast[t1map]+contrast[flat]),0,1);
										intens[n] = Numerics.max(intens[n], wf*contrast[flat] + wt*contrast[t1map]);
									} else if (t1map>-1 && t1ridges>-1) {
										intens[n] = Numerics.max(intens[n], (1.0f+contrast[t1ridges])/2.0f*atlas.getMap(atlas.pT1R)[n][0]
																			+(1.0f-contrast[t1ridges])/2.0f*contrast[t1map] );
									} else if (flat>-1 && t1ridges>-1) {
										intens[n] = Numerics.max(intens[n], (1.0f+contrast[t1ridges])/2.0f*atlas.getMap(atlas.pT1R)[n][0]
																			+(1.0f-contrast[t1ridges])/2.0f*contrast[flat] );
									}
								} else if (nc==3) {
									if (flat>-1 && t1map>-1 && t1ridges>-1) {
										float wf = Numerics.bounded((1+contrast[t1map])/(2+contrast[flat]+contrast[t1map]),0,1);
										float wt = Numerics.bounded((1+contrast[flat])/(2+contrast[t1map]+contrast[flat]),0,1);
										intens[n] = Numerics.max(intens[n], (1.0f+contrast[t1ridges])/2.0f
																				*atlas.getMap(atlas.pT1R)[n][0]
																			+(1.0f-contrast[t1ridges])/2.0f
																				*(wf*contrast[flat] + wt*contrast[t1map]) );
									}
								}
								*/
								intens[n] = 1.0f;
								
								if (flat>-1) intens[n] = Numerics.min(intens[n], contrast[flat]);
								if (t1map>-1) intens[n] = Numerics.min(intens[n], contrast[t1map]);
								
								if (t1ridges>-1) intens[n] = (1.0f+contrast[t1ridges])/2.0f*atlas.getMap(atlas.sPV)[n][0]
															 +(1.0f-contrast[t1ridges])/2.0f*intens[n];
							}
						}
					}
				}
			}
			
			// combine with shape	
			for (int n=0;n<nobj;n++) {			
				if (n==0) val = ImageInterpolation.linearInterpolation(atlas.getShape(n),1.0f,XA[0],XA[1],XA[2],nax,nay,naz);
				else val = ImageInterpolation.linearInterpolation(atlas.getShape(n),0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
				shape = (val - sigmaS)/(val + sigmaS - 2.0f*sigmaS*val);
				
				// use gain composition
				float w1 = Numerics.bounded((1+shape)/(2+shape+intens[n]),0,1);
				float w2 = Numerics.bounded((1+intens[n])/(2+shape+intens[n]),0,1);
				intens[n] = w1*intens[n] + w2*shape;
			}
			// keep only the napprox best
			for (int n=0;n<napprox;n++) {
				byte nmax=0;
				for (byte m=1;m<nobj;m++) if (intens[m]>intens[nmax]) {
					nmax = m;
				}
				bestlabel[n][xyzs] = nmax;
				bestgain[n][xyzs] = intens[nmax];
				intens[nmax] = -1.0f;
			}
			/*
			for (int n=0;n<nobj;n++) gain[n][xyzs] = intens[n];
			*/
		}  
    }

    public final void normalizeBestGainFunctions() {
        for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
        	float sum = 0.0f;
        	for (int n=0;n<napprox;n++) {
        		sum += 1+bestgain[n][xyz];
        	}
        	for (int n=0;n<napprox;n++) {
        		// sum over everything ? assumes only values with p>1/2 are used...
        		bestgain[n][xyz] = 2.0f*(1.0f+bestgain[n][xyz])/sum-1.0f;
        	}
        }    	
    	return;	
    }
    
	/**
	 *	compute the intensity priors given the atlas
	 */
    public final void computeApproxPartialVolumes(byte[][] mgdmlb, float[][] mgdmlv, int nlb, float pvscale) {
    	float[] posterior = new float[nobj];
    	float[] prior = new float[nobj];
    	float sum;
    	
    	for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
    		int xyzs = x + nsx*y + nsx*nsy*z;
			
    		for (int n=0;n<nobj;n++) posterior[n] = 0.0f;
    		sum = 0.0f;
    		for (int n=0;n<nlb;n++) if (mgdmlb[n][xyzs]>-1) {
    			
    			if (n==0) {
    				if (mgdmlv[n][xyzs]>=0)	sum += mgdmlv[n][xyzs];
    				else sum += pvscale;
    				posterior[mgdmlb[n][xyzs]] = 0.5f + 0.5f*Numerics.min(sum/pvscale, 1.0f);
    			} else {
    				posterior[mgdmlb[n][xyzs]] = 0.5f - 0.5f*Numerics.min(sum/pvscale, 1.0f);
    				if (mgdmlv[n][xyzs]>=0)	sum += mgdmlv[n][xyzs];
    				else sum += pvscale;
    			}
    			
    			//posterior[mgdmlb[n][xyzs]] = (nlb-n)/(float)nlb;
    		}
    		// also compute the other ones? or just the closest other neighbor?
    		for (int n=0;n<nobj;n++) if (posterior[n]==0) {
    			posterior[n] = 0.5f - 0.5f*Numerics.min(sum/pvscale, 1.0f);
    		}
    		
    		// get the gains
    		for (int n=0;n<nobj;n++) prior[n] = 0.0f;
    		for (int n=0;n<napprox;n++) prior[bestlabel[n][xyzs]] = Numerics.max(0.0f, bestgain[n][xyzs]);
    		
    		// combine a priori gains with a posteriori segmentations
    		sum = 0.0f;
			for (int n=0;n<nobj;n++) {			
				posterior[n] = 2.0f*prior[n]*posterior[n]/Numerics.max(0.001f,(prior[n]+posterior[n]));
				sum += posterior[n];
			}
			/*
			// build a membership function??
			if (sum>0.001) {
				for (int n=0;n<nobj;n++) gain[n][xyzs]/= sum;
			}
			*/
			// keep only the best ones
			for (int n=0;n<napprox;n++) {
				byte nmax=0;
				for (byte m=1;m<nobj;m++) if (posterior[m]>posterior[nmax]) {
					nmax = m;
				}
				bestlabel[n][xyzs] = nmax;
				bestgain[n][xyzs] = posterior[nmax];
				posterior[nmax] = -1.0f;
			}
		}  
    }

    public final void computeCurvatureWeights(float sigma) {
    	curvweight = new float[nsx*nsy*nsz];

    	int t1ridges=-1;
    	for (byte n=0;n<nc;n++) {
			if (modality[n].equals("T1ridges")) t1ridges=n;
		}
        for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
        	if (t1ridges==-1) curvweight[xyz] = 1.0f;
        	else {
        		float proba = Numerics.min(1.0f,Numerics.abs(image[t1ridges][xyz]/imrange[t1ridges]));
				curvweight[xyz] = Numerics.bounded(proba/sigma, 0.1f, 1.0f);	
        	}
        }    	
    	return;	
    }
	/**
	 *	compute the intensity priors given the atlas
	 */
	 /*
    public final void computeBoundaryForceField(float factor) {
    	field = new float[3][nsx*nsy*nsz];
    	score = new float[nsx*nsy*nsz];
    	
    	// compute the local priors
    	for (int x=1;x<nsx-1;x++) for (int y=1;y<nsy-1;y++) for (int z=1;z<nsz-1;z++) {
        	int xyzs = x + nsx*y + nsx*nsy*z;
        	
        	score[xyzs] = 0.0f;
        	for (byte k = 0; k<26; k++) {
				int xyzn = neighborIndex(k, xyzs);
		
				float mindiff = 0.0f;
				float maxdiff = 0.0f;
				for (int c=0;c<nc;c++) {
					mindiff += (min[c][xyzs]-min[c][xyzn])*(min[c][xyzs]-min[c][xyzn])/(sigmaI*sigmaI);
					maxdiff += (max[c][xyzs]-max[c][xyzn])*(max[c][xyzs]-max[c][xyzn])/(sigmaI*sigmaI);
				}
				score[xyzs] = Numerics.max(score[xyzs], mindiff, maxdiff);
			}
			score[xyzs] = 1.0f/(1.0f + factor*score[xyzs]);
		}    
    	for (int x=1;x<nsx-1;x++) for (int y=1;y<nsy-1;y++) for (int z=1;z<nsz-1;z++) {
    		int xyzs = x + nsx*y + nsx*nsy*z;
        	
    		field[0][xyzs] = 0.5f*(score[xyzs+1]-score[xyzs-1]);
    		field[1][xyzs] = 0.5f*(score[xyzs+nsx]-score[xyzs-nsx]);
    		field[2][xyzs] = 0.5f*(score[xyzs+nsx*nsy]-score[xyzs-nsx*nsy]);
    	}
     }
     */
    
	/**
	 *	propagate the gain on neighboring regions
	 */
	 /*
    public final void growFastMarching(int nmax) {
    	
		if (debug) System.out.print("evolution: fast marching\n");

    	// computation variables
        processed = new byte[nsx*nsy*nsz]; // note: using a byte instead of boolean for the second pass
    	
        // init
        boolean[] boundary = ObjectGeometry.multiObjectBoundary(segmentation, nsx, nsy, nsz, 6);
    	
    	// first estimate the binary tree size
    	int size = ObjectStatistics.volume(boundary, nsx, nsy, nsz);
    	// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	if (debug) System.out.print("init ("+size+")\n");
        
		heap.reset();
		double curr, next;
		float speed = 0;
		int nbound = 0;
		byte lb=-1;
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (mask[xyz]) {
			processed[xyz] = 0;
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (boundary[xyz]) {
				// add to the heap
				curr = gain[segmentation[xyz]][xyz];
				next = -1.0f;
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					 if (mask[xyzn]) if (segmentation[xyzn]!=segmentation[xyz]) {
						if (gain[segmentation[xyzn]][xyz]>next) {
							lb = segmentation[xyzn];
							next = gain[lb][xyz];
						}
					}
				}
				if (next>curr) {
					speed =  (float)(1.0/(alpha+next-curr));
					
					heap.addValue(speed, xyz, lb);
					
					nbound++;
				}
			}
		}
		System.out.println("initial boundary points: "+nbound);
			
		// grow the labels and functions			
		float dist = 0.0f;
		int changed = 0;
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// if the location has been processed already, do we stop??
			if (processed[xyz]>nmax)  continue;
			
			// check if the proposed gain is better than previous
			if (gain[lb][xyz]<gain[segmentation[xyz]][xyz]) continue;
			
			// check topology
			if (!homeomorphicLabeling(xyz, lb)) continue; 
			
			// update the segmentation and distance function (?) at the current level
			processed[xyz]++;
			segmentation[xyz] = lb;
			score[xyz] = dist;
			changed++;
			//System.out.print(".");
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn]) {
					if (processed[xyzn]<processed[xyz] && segmentation[xyzn]!=lb) {
						curr = gain[segmentation[xyzn]][xyzn];
						
						next = gain[lb][xyzn];
						if (next>curr) {
							// add to the heap
							speed =  (float)(1.0/(alpha+next-curr));
							
							//System.out.print("+");
							heap.addValue(dist+speed, xyzn, lb);
							nbound++;
						}
					}
				}
			}
		}
		System.out.println("total boundary points: "+nbound);
		System.out.println("total changed points: "+changed);
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (processed[xyz]==0) {
			score[xyz] = dist;
		}

		
        return;
    }
    */

	/**
	 *	propagate the gain on neighboring regions
	 */ 
    public final void evolveFastMarching(int nmax, float g0) {
    	
    	if (nmax<1) return;
    	
		if (debug) System.out.print("evolution: fast marching\n");

    	// computation variables
        processed = new byte[nsx*nsy*nsz]; // note: using a byte instead of boolean for the second pass
    	
        // init
        boolean[] boundary = ObjectGeometry.multiObjectBoundary(segmentation, nsx, nsy, nsz, 6);
    	
    	// first estimate the binary tree size
    	int size = ObjectStatistics.volume(boundary, nsx, nsy, nsz);
    	// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	if (debug) System.out.print("init ("+size+")\n");
        
		heap.reset();
		double curr, next;
		float speed = 0;
		int nbound = 0;
		byte lb=-1;
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (mask[xyz]) {
			processed[xyz] = 0;
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (boundary[xyz]) {
				// add to the heap
				curr = gain[segmentation[xyz]][xyz];
				next = -1.0f;
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					 if (mask[xyzn]) if (segmentation[xyzn]!=segmentation[xyz]) {
						if (gain[segmentation[xyzn]][xyz]>next) {
							lb = segmentation[xyzn];
							next = gain[lb][xyz];
						}
					}
				}
				if (next>curr && next>g0) {
					speed =  (float)(1.0/(alpha+next-curr));
					
					heap.addValue(speed, xyz, lb);
					
					nbound++;
				}
			}
		}
		System.out.println("initial boundary points: "+nbound);
			
		// grow the labels and functions			
		float dist = 0.0f;
		int changed = 0;
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// if the location has been processed already, do we stop??
			if (processed[xyz]>nmax)  continue;
			
			// check if the proposed gain is better than previous
			if (gain[lb][xyz]<gain[segmentation[xyz]][xyz]) continue;
			
			// check topology
			if (!homeomorphicLabeling(xyz, lb)) continue; 
			
			// update the segmentation and distance function (?) at the current level
			processed[xyz]++;
			segmentation[xyz] = lb;
			score[xyz] = dist;
			changed++;
			//System.out.print(".");
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn]) {
					if (processed[xyzn]<processed[xyz] && segmentation[xyzn]!=lb) {
						curr = gain[segmentation[xyzn]][xyzn];
						
						// propagate gain from source ??
						/*
						float diff = 0.0f;
						for (int c=0;c<nc;c++) {
							diff += Numerics.max((max[c][xyzn]-max[c][xyz])*(max[c][xyzn]-max[c][xyz]),
												  (min[c][xyzn]-min[c][xyz])*(min[c][xyzn]-min[c][xyz]));
						}
						float similarity = 1.0f/( 1.0f + (diff)/(nc*sigmaI*sigmaI) );
		
						next = Numerics.max(gain[lb][xyzn], similarity*gain[lb][xyz]);
						gain[lb][xyzn] = (float)next;
						*/
						next = gain[lb][xyzn];
						if (next>curr && next>g0) {
							// add to the heap
							speed =  (float)(1.0/(alpha+next-curr));
							
							//System.out.print("+");
							heap.addValue(dist+speed, xyzn, lb);
							nbound++;
						}
					}
				}
			}
		}
		System.out.println("total boundary points: "+nbound);
		System.out.println("total changed points: "+changed);
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (processed[xyz]==0) {
			score[xyz] = dist;
		}

		
        return;
    }

	/**
	 *	propagate the gain on neighboring regions
	 */
    public final void evolveShapeMarching(int nmax, float p0) {
    	
    	if (nmax<1) return;
    	
		if (debug) System.out.print("evolution: fast marching\n");

    	// computation variables
        processed = new byte[nsx*nsy*nsz]; // note: using a byte instead of boolean for the second pass
    	
        // init
        boolean[] boundary = ObjectGeometry.multiObjectBoundary(segmentation, nsx, nsy, nsz, 6);
    	
    	// first estimate the binary tree size
    	int size = ObjectStatistics.volume(boundary, nsx, nsy, nsz);
    	// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MAXTREE);
    	if (debug) System.out.print("init ("+size+")\n");
        
		heap.reset();
		double curr, next;
		float speed = 0;
		int nbound = 0;
		byte lb=-1;
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (mask[xyz]) {
			processed[xyz] = 0;
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (boundary[xyz]) {
				// add to the heap
				curr = getAtlasPrior(segmentation[xyz], xyz);
				next = -1.0f;
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					 if (mask[xyzn]) if (segmentation[xyzn]!=segmentation[xyz]) {
						if (getAtlasPrior(segmentation[xyzn], xyz)>next) {
							lb = segmentation[xyzn];
							next = getAtlasPrior(segmentation[xyzn], xyz);
						}
					}
				}
				if (next>curr && next>p0) {
					//speed =  (float)(1.0/(alpha+next-curr));
					speed =  (float)next;
					
					heap.addValue(speed, xyz, lb);
					
					nbound++;
				}
			}
		}
		System.out.println("initial boundary points: "+nbound);
			
		// grow the labels and functions			
		float dist = 0.0f;
		int changed = 0;
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// if the location has been processed already, do we stop??
			if (processed[xyz]>nmax)  continue;
			
			// check if the proposed gain is better than previous
			if (getAtlasPrior(lb, xyz)<getAtlasPrior(segmentation[xyz], xyz)) continue;
			
			// check topology
			if (!homeomorphicLabeling(xyz, lb)) continue; 
			
			// update the segmentation and distance function (?) at the current level
			processed[xyz]++;
			segmentation[xyz] = lb;
			score[xyz] = dist;
			changed++;
			//System.out.print(".");
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn]) {
					if (processed[xyzn]<processed[xyz] && segmentation[xyzn]!=lb) {
						curr = getAtlasPrior(segmentation[xyzn], xyzn);
						next = getAtlasPrior(lb, xyzn);
						
						if (next>curr && next>p0) {
							// add to the heap
							//speed =  (float)(1.0/(alpha+next-curr));
							speed = (float) next;
							
							//System.out.print("+");
							heap.addValue(dist*speed, xyzn, lb);
							nbound++;
						}
					}
				}
			}
		}
		System.out.println("total boundary points: "+nbound);
		System.out.println("total changed points: "+changed);
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (processed[xyz]==0) {
			score[xyz] = dist;
		}

		
        return;
    }

	/**
	 *	propagate the gain on neighboring regions
	 */
    public final void evolveWarpMarching(float[][] forces) {
    	
    	if (debug) System.out.print("evolution: fast marching warping\n");

    	// computation variables
        processed = new byte[nsx*nsy*nsz]; // note: using a byte instead of boolean for the second pass
    	
        // init
        boolean[] boundary = ObjectGeometry.multiObjectBoundary(segmentation, nsx, nsy, nsz, 6);
    	
    	// first estimate the binary tree size
    	int size = ObjectStatistics.volume(boundary, nsx, nsy, nsz);
    	// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	if (debug) System.out.print("init ("+size+")\n");
        
		heap.reset();
		double curr, next;
		float speed = 0;
		int nbound = 0;
		byte lb=-1;
		float[][] warp = new float[3][nsx*nsy*nsz];
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) {
			int x,y,z;
			z = xyz/(nsx*nsy);
			y = (xyz-nsx*nsy*z)/nsx;
			x = xyz-nsx*nsy*z-nsx*y;
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x/scaling, y/scaling, z/scaling);
			
			warp[X][xyz] = ImageInterpolation.linearInterpolation(forces[X],0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
			warp[Y][xyz] = ImageInterpolation.linearInterpolation(forces[Y],0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
			warp[Z][xyz] = ImageInterpolation.linearInterpolation(forces[Z],0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
		}
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (mask[xyz]) {
			processed[xyz] = 0;
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (boundary[xyz]) {
				// add to the heap
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (mask[xyzn]) if (segmentation[xyzn]!=segmentation[xyz]) {
						if (k<3) {
							int d = k;
							if (warp[d][xyz]>0.5f/scaling) {
								lb = segmentation[xyz];
								next = warp[d][xyz];
								speed =  (float)(1.0/(alpha+next));
								warp[d][xyz] -= 1.0f/scaling;
								heap.addValue(speed, xyzn, lb);
								nbound++;
							}
						} else {
							int d = k-3;
							if (warp[d][xyz]<-0.5f/scaling) {
								lb = segmentation[xyz];
								next = -warp[d][xyz];
								speed =  (float)(1.0/(alpha+next));
								warp[d][xyz] += 1.0f/scaling;
								heap.addValue(speed, xyzn, lb);
								nbound++;
							}
						}							
					}
				}
			}
		}
		System.out.println("initial boundary points: "+nbound);
			
		// grow the labels and functions			
		float dist = 0.0f;
		int changed = 0;
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// check topology
			if (!homeomorphicLabeling(xyz, lb)) continue; 
			
			// update the segmentation and distance function (?) at the current level
			processed[xyz]++;
			segmentation[xyz] = lb;
			score[xyz] = dist;
			changed++;
			//System.out.print(".");
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn]) {
					if (processed[xyzn]<processed[xyz] && segmentation[xyzn]!=lb) {
						if (k<3) {
							int d = k;
							if (warp[d][xyz]>0.5f/scaling) {
								lb = segmentation[xyz];
								next = warp[d][xyz];
								speed =  (float)(1.0/(alpha+next));
								warp[d][xyz] -= 1.0f/scaling;
								heap.addValue(speed, xyzn, lb);
							}
						} else {
							int d = k-3;
							if (warp[d][xyz]<-0.5f/scaling) {
								lb = segmentation[xyz];
								next = -warp[d][xyz];
								speed =  (float)(1.0/(alpha+next));
								warp[d][xyz] += 1.0f/scaling;
								heap.addValue(speed, xyzn, lb);
							}
						}
					}
				}
			}
		}
		System.out.println("total boundary points: "+nbound);
		System.out.println("total changed points: "+changed);
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (processed[xyz]==0) {
			score[xyz] = dist;
		}
		
        return;
    }

	/**
	 *	propagate the gain on neighboring regions
	 */
    public final void growFastMarching(int nmax, int g0) {
    	
    	if (nmax<1) return;
    	
		if (debug) System.out.print("evolution: fast marching\n");

    	// computation variables
        processed = new byte[nsx*nsy*nsz]; // note: using a byte instead of boolean for the second pass
    	
        // init
        boolean[] boundary = ObjectGeometry.multiObjectBoundary(segmentation, nsx, nsy, nsz, 6);
    	
    	// first estimate the binary tree size
    	int size = ObjectStatistics.volume(boundary, nsx, nsy, nsz);
    	// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	if (debug) System.out.print("init ("+size+")\n");
        
		heap.reset();
		double curr, next;
		float speed = 0;
		int nbound = 0;
		byte lb=-1;
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (mask[xyz]) {
			processed[xyz] = 0;
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (boundary[xyz]) {
				// add to the heap
				curr = gain[segmentation[xyz]][xyz];
				next = -1.0f;
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					 if (mask[xyzn]) if (segmentation[xyzn]!=segmentation[xyz]) {
						if (gain[segmentation[xyzn]][xyz]>next) {
							lb = segmentation[xyzn];
							next = gain[lb][xyz];
						}
					}
				}
				if (next>curr && next>g0) {
					speed =  (float)(1.0/(alpha+next-curr));
					
					heap.addValue(speed, xyz, lb);
					
					nbound++;
				}
			}
		}
		System.out.println("initial boundary points: "+nbound);
			
		// grow the labels and functions			
		float dist = 0.0f;
		int changed = 0;
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// if the location has been processed already, do we stop??
			if (processed[xyz]>nmax)  continue;
			
			// check if the proposed gain is better than previous
			if (gain[lb][xyz]<gain[segmentation[xyz]][xyz]) continue;
			
			// check topology
			if (!homeomorphicLabeling(xyz, lb)) continue; 
			
			// update the segmentation and distance function (?) at the current level
			processed[xyz]++;
			segmentation[xyz] = lb;
			score[xyz] = dist;
			changed++;
			//System.out.print(".");
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn]) {
					if (processed[xyzn]<processed[xyz] && segmentation[xyzn]!=lb) {
						curr = gain[segmentation[xyzn]][xyzn];
						
						next = gain[lb][xyzn];
						if (next>curr && next>g0) {
							// add to the heap
							speed =  (float)(1.0/(alpha+next-curr));
							
							//System.out.print("+");
							heap.addValue(dist+speed, xyzn, lb);
							nbound++;
						}
					}
				}
			}
		}
		System.out.println("total boundary points: "+nbound);
		System.out.println("total changed points: "+changed);
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (processed[xyz]==0) {
			score[xyz] = dist;
		}

		
        return;
    }

	/**
	 *	propagate the gain on neighboring regions
	 */
    public final void growFastMarching2(int nmax, float g0) {
    	
    	if (nmax<1) return;
    	
		if (debug) System.out.print("evolution: fast marching\n");

    	// computation variables
        processed = new byte[nsx*nsy*nsz]; // note: using a byte instead of boolean for the second pass
    	
        // init
        boolean[] boundary = ObjectGeometry.multiObjectBoundary(segmentation, nsx, nsy, nsz, 6);
    	
    	// first estimate the binary tree size
    	int size = ObjectStatistics.volume(boundary, nsx, nsy, nsz);
    	// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	if (debug) System.out.print("init ("+size+")\n");
        
		heap.reset();
		double curr, next;
		float speed = 0;
		int nbound = 0;
		byte lb=-1;
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (mask[xyz]) {
			processed[xyz] = 0;
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (boundary[xyz]) {
				// add to the heap
				curr = gain[segmentation[xyz]][xyz];
				next = -1.0f;
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					 if (mask[xyzn]) if (segmentation[xyzn]!=segmentation[xyz]) {
						if (gain[segmentation[xyzn]][xyz]>next) {
							lb = segmentation[xyzn];
							next = gain[lb][xyz];
						}
					}
				}
				if (next>curr && next>g0) {
					speed =  (float)(1.0/(alpha+next-curr));
					
					heap.addValue(speed, xyz, lb);
					
					nbound++;
				}
			}
		}
		System.out.println("initial boundary points: "+nbound);
			
		// grow the labels and functions			
		float dist = 0.0f;
		int changed = 0;
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// if the location has been processed already, do we stop??
			if (processed[xyz]>nmax)  continue;
			
			// check if the proposed gain is better than previous
			if (gain[lb][xyz]<gain[segmentation[xyz]][xyz]) continue;
			
			// check topology
			if (!homeomorphicLabeling(xyz, lb)) continue; 
			
			// update the segmentation and distance function (?) at the current level
			processed[xyz]++;
			segmentation[xyz] = lb;
			score[xyz] = dist;
			changed++;
			//System.out.print(".");
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn]) {
					if (processed[xyzn]<processed[xyz] && segmentation[xyzn]!=lb) {
						curr = gain[segmentation[xyzn]][xyzn];
						
						next = gain[lb][xyzn];
						if (next>curr && next>g0) {
							// add to the heap
							speed =  (float)(1.0/(alpha+next-curr));
							
							//System.out.print("+");
							heap.addValue(dist+speed, xyzn, lb);
							nbound++;
						}
					}
				}
			}
		}
		System.out.println("total boundary points: "+nbound);
		System.out.println("total changed points: "+changed);
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (processed[xyz]==0) {
			score[xyz] = dist;
		}

		
        return;
    }

	/**
	 *	propagate the gain on neighboring regions
	 */
    public final void shrinkFastMarching(int nmax, float g0) {
    	
    	if (nmax<1) return;
    	
		if (debug) System.out.print("evolution: fast marching\n");

    	// computation variables
        processed = new byte[nsx*nsy*nsz]; // note: using a byte instead of boolean for the second pass
    	
        // init
        boolean[] boundary = ObjectGeometry.multiObjectBoundary(segmentation, nsx, nsy, nsz, 6);
    	
    	// first estimate the binary tree size
    	int size = ObjectStatistics.volume(boundary, nsx, nsy, nsz);
    	// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	if (debug) System.out.print("init ("+size+")\n");
        
		heap.reset();
		double curr, next;
		float speed = 0;
		int nbound = 0;
		byte lb=-1;
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (mask[xyz]) {
			processed[xyz] = 0;
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (boundary[xyz]) {
				// add to the heap
				curr = gain[segmentation[xyz]][xyz];
				if (curr<g0) {
					// find largest gain
					next = -1.0f;
					for (int n = 0; n<nobj; n++) {
						if (gain[n][xyz]>next) {
							next = gain[n][xyz];
						}
					}
					if (next>curr) {
						speed =  (float)(1.0/(alpha+next-curr));
						
						// find best neighbor
						next = -1.0f;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							 if (mask[xyzn]) if (segmentation[xyzn]!=segmentation[xyz]) {
								if (gain[segmentation[xyzn]][xyz]>next) {
									lb = segmentation[xyzn];
									next = gain[lb][xyz];
								}
							}
						}
						
						heap.addValue(speed, xyz, lb);
						
						nbound++;
					}
				}
			}
		}
		System.out.println("initial boundary points: "+nbound);
		
		// grow the labels and functions			
		float dist = 0.0f;
		int changed = 0;
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// if the location has been processed already, do we stop??
			if (processed[xyz]>nmax)  continue;
			
			// check if the proposed gain is better than previous
			if (gain[lb][xyz]<gain[segmentation[xyz]][xyz]) continue;
			
			// check topology
			if (!homeomorphicLabeling(xyz, lb)) continue; 
			
			// update the segmentation and distance function (?) at the current level
			processed[xyz]++;
			segmentation[xyz] = lb;
			score[xyz] = dist;
			changed++;
			//System.out.print(".");
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn]) {
					if (processed[xyzn]<processed[xyz] && segmentation[xyzn]!=lb) {
						curr = gain[segmentation[xyzn]][xyzn];
						
						/*
						// propagate gain from source ??
						float diff = 0.0f;
						for (int c=0;c<nc;c++) {
							diff += Numerics.max((max[c][xyzn]-max[c][xyz])*(max[c][xyzn]-max[c][xyz]),
												  (min[c][xyzn]-min[c][xyz])*(min[c][xyzn]-min[c][xyz]));
						}
						float similarity = 1.0f/( 1.0f + (diff)/(nc*sigmaI*sigmaI) );
		
						next = Numerics.max(gain[lb][xyzn], similarity*gain[lb][xyz]);
						gain[lb][xyzn] = (float)next;
						*/
						next = -1.0f;
						for (int n = 0; n<nobj; n++) {
							if (gain[n][xyzn]>next) {
								next = gain[n][xyzn];
							}
						}
						//next = gain[lb][xyzn];
						if (next>curr) {
							// add to the heap
							speed =  (float)(1.0/(alpha+next-curr));
							
							//System.out.print("+");
							heap.addValue(dist+speed, xyzn, lb);
							nbound++;
						}
					}
				}
			}
		}
		System.out.println("total boundary points: "+nbound);
		System.out.println("total changed points: "+changed);
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (processed[xyz]==0) {
			score[xyz] = dist;
		}

		
        return;
    }

	/**
	 *	propagate the gain on neighboring regions
	 */
    public final void shrinkFastMarching2(int nmax, float g0) {
    	
    	if (nmax<1) return;
    	
		if (debug) System.out.print("evolution: fast marching\n");

    	// computation variables
        processed = new byte[nsx*nsy*nsz]; // note: using a byte instead of boolean for the second pass
    	
        // init
        boolean[] boundary = ObjectGeometry.multiObjectBoundary(segmentation, nsx, nsy, nsz, 6);
    	
    	// first estimate the binary tree size
    	int size = ObjectStatistics.volume(boundary, nsx, nsy, nsz);
    	// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	if (debug) System.out.print("init ("+size+")\n");
        
		heap.reset();
		double curr, next;
		float speed = 0;
		int nbound = 0;
		byte lb=-1;
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (mask[xyz]) {
			processed[xyz] = 0;
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (boundary[xyz]) {
				// add to the heap
				curr = gain[segmentation[xyz]][xyz];
				if (curr<g0) {
					// find best neighbor
					next = -1.0f;
					for (int k = 0; k<6; k++) {
						int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
						 if (mask[xyzn]) if (segmentation[xyzn]!=segmentation[xyz]) {
							if (gain[segmentation[xyzn]][xyz]>next) {
								lb = segmentation[xyzn];
								next = gain[lb][xyz];
							}
						}
					}
					// replace current by next best neighbor
					if (next>curr) {
						speed =  (float)(1.0/(alpha+next-curr));
						heap.addValue(speed, xyz, lb);
						
						nbound++;
					}
				}
			}
		}
		System.out.println("initial boundary points: "+nbound);
		
		// grow the labels and functions			
		float dist = 0.0f;
		int changed = 0;
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// if the location has been processed already, do we stop??
			if (processed[xyz]>nmax)  continue;
			
			// check if the proposed gain is better than previous
			if (gain[lb][xyz]<gain[segmentation[xyz]][xyz]) continue;
			
			// check topology
			if (!homeomorphicLabeling(xyz, lb)) continue; 
			
			// update the segmentation and distance function (?) at the current level
			processed[xyz]++;
			segmentation[xyz] = lb;
			score[xyz] = dist;
			changed++;
			//System.out.print(".");
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn]) {
					if (processed[xyzn]<processed[xyz]) {
						curr = gain[segmentation[xyzn]][xyzn];
						if (curr<g0) {
							next = -1.0f;
							for (int l = 0; l<6; l++) {
								int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
								 if (mask[xyznb]) if (segmentation[xyznb]!=segmentation[xyzn]) {
									if (gain[segmentation[xyznb]][xyzn]>next) {
										lb = segmentation[xyznb];
										next = gain[lb][xyzn];
									}
								}
							}
							if (next>curr) {
								// add to the heap
								speed =  (float)(1.0/(alpha+next-curr));
							
								//System.out.print("+");
								heap.addValue(dist+speed, xyzn, lb);
								nbound++;
							}
						}
					}
				}
			}
		}
		System.out.println("total boundary points: "+nbound);
		System.out.println("total changed points: "+changed);
		
		for (int xyz = 0; xyz<nsx*nsy*nsz; xyz++) if (processed[xyz]==0) {
			score[xyz] = dist;
		}

		
        return;
    }

	private final int neighborIndex(byte d, int id) {
		switch (d) {
			case pX		: 	return id+1; 		
			case mX		:	return id-1;
			case pY		:	return id+nsx;
			case mY		:	return id-nsx;
			case pZ		:	return id+nsx*nsy;
			case mZ		:	return id-nsx*nsy;
			case pXpY	:	return id+1+nsx;
			case mXpY	:	return id-1+nsx;
			case pYpZ	:	return id+nsx+nsx*nsx;
			case mYpZ	:	return id-nsx+nsx*nsy;
			case pZpX	:	return id+nsx*nsy+1;	
			case mZpX	:	return id-nsx*nsy+1;
			case pXmY	:	return id+1-nsx;	
			case mXmY	:	return id-1-nsx;
			case pYmZ	:	return id+nsx-nsx*nsy;
			case mYmZ	:	return id-nsx-nsx*nsy;
			case pZmX	:	return id+nsx*nsy-1;
			case mZmX	:	return id-nsx*nsy-1;
			case pXpYpZ	:	return id+1+nsx+nsx*nsy;
			case mXpYpZ	:	return id-1+nsx+nsx*nsy;
			case pXmYmZ	:	return id+1-nsx-nsx*nsy; 
			case mXmYmZ	:	return id-1-nsx-nsx*nsy;
			case pXmYpZ	:	return id+1-nsx+nsx*nsy;
			case mXmYpZ	:	return id-1-nsx+nsx*nsy;
			case pXpYmZ	:	return id+1+nsx-nsx*nsy; 
			case mXpYmZ	:	return id-1+nsx-nsx*nsy;
			default		:	return id;
		}
	}
	
    /**
	 *  critical relation detection: groups objects with relations
	 */
    private final boolean homeomorphicLabeling(int xyz, byte lb) {
    	// if we don't check, just exit
    	if (!checkTopology) return true;
    	
		// is the new segmentation homeomorphic ? 
		
		// inside the original object ?
		if (segmentation[xyz]==lb) return true;
		
		boolean [][][] obj = new boolean[3][3][3];
		
		// does it change the topology of the new object ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lb) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = true;
		
		if (checkComposed) if (!ObjectStatistics.isWellComposed(obj,1,1,1)) return false;		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			
		// does it change the topology of the object it modifies ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (segmentation[xyz+i+j*nsx+l*nsx*nsy]==segmentation[xyz]) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = false;
		
		if (checkComposed) if (!ObjectStatistics.isWellComposed(obj,1,1,1)) return false;		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;

		// does it change the topology of a relation between the modified object and its neighbors ?
		int  Nconfiguration = 0;
		short[] lbs = new short[26];
		for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) {
			if ( (i*i+j*j+l*l>0) 
				&& (segmentation[xyz+i+j*nsx+l*nsx*nsy]!=lb) 
				&& (segmentation[xyz+i+j*nsx+l*nsx*nsy]!=segmentation[xyz]) ) {
				boolean found = false;
				for (int n=0;n<Nconfiguration;n++) 
					if (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lbs[n]) { found = true; break; }
				
				if (!found) {
					lbs[Nconfiguration] = segmentation[xyz+i+j*nsx+l*nsx*nsy];
					Nconfiguration++;
				}
			}
		}
		// pairs

		for (int n=0;n<Nconfiguration;n++) {
			// in relation with previous object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (segmentation[xyz+i+j*nsx+l*nsx*nsy]==segmentation[xyz])
					|| (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lbs[n]) ) {
					obj[1+i][1+j][1+l] = true;
				} else {
					obj[1+i][1+j][1+l] = false;
				}
			}
			obj[1][1][1] = false;
			if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		}
		for (int n=0;n<Nconfiguration;n++) {
			// in relation with new object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lb)
					|| (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lbs[n]) ) {
					obj[1+i][1+j][1+l] = true;
				} else {
					obj[1+i][1+j][1+l] = false;
				}
			}
			obj[1][1][1] = true;
			if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		}

		// triplets
		for (int n=0;n<Nconfiguration;n++) {
			for (int m=n+1;m<Nconfiguration;m++) {
				// in relation with previous object
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (segmentation[xyz+i+j*nsx+l*nsx*nsy]==segmentation[xyz])
						|| (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lbs[n])
						|| (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lbs[m]) ) {
						obj[1+i][1+j][1+l] = true;
					} else {
						obj[1+i][1+j][1+l] = false;
					}
				}
				obj[1][1][1] = false;
				if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			}
		}
		for (int n=0;n<Nconfiguration;n++) {
			for (int m=n+1;m<Nconfiguration;m++) {
				// in relation with new object
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lb)
						|| (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lbs[n]) 
						|| (segmentation[xyz+i+j*nsx+l*nsx*nsy]==lbs[m]) ) {
						obj[1+i][1+j][1+l] = true;
					} else {
						obj[1+i][1+j][1+l] = false;
					}
				}
				obj[1][1][1] = true;
				if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			}
		}

		// else, it works
		return true;
    }
    
	private final float getAtlasPrior(byte lb, int xyz) {
		int x,y,z;
		z = xyz/(nsx*nsy);
		y = (xyz-nsx*nsy*z)/nsx;
		x = xyz-nsx*nsy*z-nsx*y;
		float[] XA = new float[3];
		atlas.imageToShapeCoordinates(XA, x/scaling, y/scaling, z/scaling);
 		
		return ImageInterpolation.linearInterpolation(atlas.getShape(lb),0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
	}
    
	public final int[] getAtlasTopology() {
		int[] result = new int[nsx*nsy*nsz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyz = x + nsx*y + nsx*nsy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x/scaling, y/scaling, z/scaling);
			
			result[xyz] = ImageInterpolation.nearestNeighborInterpolation(atlas.getTemplate(),(byte)1,XA[0],XA[1],XA[2],nax,nay,naz);
		}
		return result;
	}
    
	private final void getAtlasWarping(float[] warp, int xyz, float[][] forces) {
		int x,y,z;
		z = xyz/(nsx*nsy);
		y = (xyz-nsx*nsy*z)/nsx;
		x = xyz-nsx*nsy*z-nsx*y;
		float[] XA = new float[3];
		atlas.imageToShapeCoordinates(XA, x/scaling, y/scaling, z/scaling);
 		
		warp[X] = ImageInterpolation.linearInterpolation(forces[X],0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
		warp[Y] = ImageInterpolation.linearInterpolation(forces[Y],0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
		warp[Z] = ImageInterpolation.linearInterpolation(forces[Z],0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
		
		return;
	}
    
    public final float[] exportAtlasTopology() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
        	float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x, y, z);
			result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(atlas.getTemplate(), (byte)1,XA[0], XA[1], XA[2], nax, nay, naz); 
		}
		return result;
    }
    /*
    public final float[] exportMin(int c) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
        	result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(min[c], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
     public final float[] exportMax(int c) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(max[c], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    */
     public final float[] exportMap(float[] map) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(map, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
     public final float[] exportMap(byte[] map) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(map, (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
     public final float[] exportMap(int[] map) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(map, (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    /*
    public final float[] exportGain(int n) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(gain[n], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
    public final float[][] exportGains() {
    	float[][] result = new float[nobj][nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) for (int n=0;n<nobj;n++) {
    		result[n][x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(gain[n], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    */
    /*
    public final byte[] computeGainSegmentation(float threshold) {
    	byte[] result = new byte[nsx*nsy*nsz];
    	
    	byte id;
    	float val,best;
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyz = x + nsx*y + nsx*nsy*z;
			id = 0;
    		best = gain[0][xyz]; 
    		for (byte n=1;n<nobj;n++) {
    			val = gain[n][xyz]; 
    			if (val>best) { best = val; id = n; }
    		}
    		if (best>threshold)	result[xyz] = id;
    		else result[xyz] = -1;
		}
		return result;
    }
    */
    /*
    public final float[] exportGainSegmentation(float threshold) {
    	float[] result = new float[nix*niy*niz];
    	
    	int id;
    	float val,best;
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		id = 0;
    		best = ImageInterpolation.linearClosestInterpolation(gain[0], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
    		for (int n=1;n<nobj;n++) {
    			val = ImageInterpolation.linearClosestInterpolation(gain[n], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
    			if (val>best) { best = val; id = n; }
    		}
    		if (best>threshold)	result[x+nix*y+nix*niy*z] = atlas.getLabels()[id];
    		else result[x+nix*y+nix*niy*z] = -1;
		}
		return result;
    }
    
    public final float[] exportGainSegmentation() {
    	float[] result = new float[nix*niy*niz];
    	
    	int id;
    	float val,best;
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		id = 0;
    		best = ImageInterpolation.linearClosestInterpolation(gain[0], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
    		for (int n=1;n<nobj;n++) {
    			val = ImageInterpolation.linearClosestInterpolation(gain[n], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
    			if (val>best) { best = val; id = n; }
    		}
    		result[x+nix*y+nix*niy*z] = atlas.getLabels()[id];
    	}
		return result;
    }
    */
    public final float[] exportBestGain(int obj) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		
    		result[x+nix*y+nix*niy*z] = 0.0f;
    		for (byte n=0;n<napprox;n++) {
				byte lb = ImageInterpolation.nearestNeighborInterpolation(bestlabel[n], (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz);
				if (lb==obj) result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(bestgain[n], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz);
			}
		}
		return result;
    }
    
    public final float[][] exportBestGains() {
    	float[][] result = new float[nobj][nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) for (int obj=0;obj<nobj;obj++) {
     		result[obj][x+nix*y+nix*niy*z] = 0.0f;
    		for (byte n=0;n<napprox;n++) {
				byte lb = ImageInterpolation.nearestNeighborInterpolation(bestlabel[n], (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz);
				if (lb==obj) result[obj][x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(bestgain[n], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz);
			}
		}
		return result;
    }
    public final float[] exportBestGainSegmentation() {
    	float[] result = new float[nix*niy*niz];
    	
    	int id;
         for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		id = ImageInterpolation.nearestNeighborInterpolation(bestlabel[0], (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
    		result[x+nix*y+nix*niy*z] = atlas.getLabels()[id];
    	}
		return result;
    }
    
    public final float[] exportProba(byte n) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		int xyz = Numerics.floor(x*scaling) + nsx*Numerics.floor(y*scaling) + nsx*nsy*Numerics.floor(z*scaling);
    		result[x+nix*y+nix*niy*z] = getAtlasPrior(n, xyz);
		}
		return result;
    }
    
    public final float[] exportSegmentation() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(segmentation, (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
 
    public final int[] rescaleSegmentation() {
    	int[] result = new int[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(segmentation, (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
 
    public final int[] labelSegmentation() {
    	int[] result = new int[nix*niy*niz];
    	
    	for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
    		if (segmentation[xyz]!=-1) {
				result[xyz] = atlas.getLabels()[segmentation[xyz]];
			} else {
				result[xyz] = 0;
			}
		}
		return result;
    }
 
    public final float[] exportScore() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(score, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
    public final float[] exportProcessed() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(processed, (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
    public final float[] exportField(int n) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(field[n], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
    public final float[] exportImage(float[] img) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(img, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
    public final float[] exportImage(byte[] img) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(img, (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    

}
