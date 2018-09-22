package de.mpg.cbs.core.intensity;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This algorithm uses a simple fast marching method to perform
 *  phase unwrapping
 *
 *	@version    Mai 2018
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class IntensityFastMarchingUnwrapping {
	
	// object types

	private	static	final	byte	EMPTY = -1;
	private	static	final	byte	OBJ = 1;
	private	static	final	byte	BG = 0;
	
	// fast marching flags
	private final static byte X = 0;
    private final static byte Y = 1;
    private final static byte Z = 2;
    	
    // numerical quantities
	private final static float	UNKNOWN=1e15f;
	
	//private static int[] xoff;
    //private static int[] yoff;
    //private static int[] zoff;

	// data and membership buffers
	private 	float[] 		phase;  			// original phase image
	private		boolean[]		mask;				// masking regions not used in computations
    //private 	float[] 		magnitude;  		// original magintude image
	private 	int[] 			wrapcount;   	    // number of wraps counted
	private static	int 		nx,ny,nz, nxyz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		BinaryHeapPair	heap;				// the heap used in fast marching
	
	// parameters
	private float[] inputImage;
	//private float[] magImage;
	private int[] maskImage = null;
	private String postprocess = "none";
	public  static final String[] postProcessTypes = {"none", "TV-approximation", "TV-residuals"};;
    private     int             nquadrant = 3;      // number of quadrants to seed the unwrapping
    private     float           tvscale = 0.0f;     // factor for optional post-processing
	
	private float[] correctImage;
	private int[]   countImage;
    //private int[]   labelImage;
    private float[] scoreImage;

	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;
	

	public final void setPhaseImage(float[] val) { inputImage = val; }
	//public final void setMagnitudeImage(float[] val) { magImage = val; }
	public final void setMaskImage(int[] val) { maskImage = val; }
	
	public final void setTVScale(float val) { tvscale = val; }
	public final void setTVPostProcessing(String val) { postprocess = val; }
	public final void setQuadrantNumber(int val) { nquadrant = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	
	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity.devel"; }
	public final String getLabel() { return "Fast Marching Unwrapping"; }
	public final String getName() { return "FastMarchingUnwrapping"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Fast marching method for unwrapping phase images, based on (Abdul-Rahman et al., 2005)"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.3"; };

	// create outputs
	public final float[] getCorrectedImage() { return correctImage; }
    public final int[] getCountImage() { return countImage; }
	//public final int[] getLabelImage() { return labelImage; }
	public final float[] getScoreImage() { return scoreImage; }
    
	public final void execute() {
	    
		boolean[][] boundaries;
		int[] label;
		float[] reliability;
		try {
		    mask = new boolean[nx*ny*nz];
		    if (maskImage==null) {
                // no mask
                for (int xyz=0;xyz<nxyz;xyz++) mask[xyz] = true;
            } else {
                for (int xyz=0;xyz<nxyz;xyz++) if (maskImage[xyz]>0) mask[xyz] = true;
                maskImage = null;
            }

            boundaries = new boolean[6][nxyz];
            for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
                int xyz = x+nx*y+nx*ny*z;
                if (x==0) boundaries[Ngb.mX][xyz] = true;
                if (y==0) boundaries[Ngb.mY][xyz] = true;
                if (z==0) boundaries[Ngb.mZ][xyz] = true;
                if (x==nx-1) boundaries[Ngb.pX][xyz] = true;
                if (y==ny-1) boundaries[Ngb.pY][xyz] = true;
                if (z==nz-1) boundaries[Ngb.pZ][xyz] = true;
            }

            wrapcount = new int[nxyz];
            label = new int[nxyz];
            reliability = new float[nxyz];
            
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeapPair(nx*ny+ny*nz+nz*nx, BinaryHeapPair.MAXTREE);
			
		} catch (OutOfMemoryError e){
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initialization\n");
		
		// get the min-max scale of the phase image
		float Pmin = 1e13f;
		float Pmax = -1e13f;
		//float Mmax = -1e13f;
		phase = inputImage;
		//magnitude = magImage;
		for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
		    if (phase[xyz]<Pmin) Pmin = phase[xyz];
		    if (phase[xyz]>Pmax) Pmax = phase[xyz];
		    //if (magnitude[xyz]>Mmax) Mmax = magnitude[xyz];
		}
		//double Pscale = (Pmax-Pmin)/(2.0*FastMath.PI);
		
		// simplify: bring everything in [0,1] internally
		for (int xyz=0;xyz<nxyz;xyz++) phase[xyz] = (phase[xyz]-Pmin)/(Pmax-Pmin);
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("initialization\n");	

        double Rmean = 0.0;
        double Rnum = 0.0;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
		    int xyz = x+nx*y+nx*ny*z;
		    if (mask[xyz]) {
		        double dphix = Numerics.modulo(phase[xyz-1]-phase[xyz], 1.0f)     - Numerics.modulo(phase[xyz]-phase[xyz+1], 1.0f);
		        double dphiy = Numerics.modulo(phase[xyz-nx]-phase[xyz], 1.0f)    - Numerics.modulo(phase[xyz]-phase[xyz+nx], 1.0f);
		        double dphiz = Numerics.modulo(phase[xyz-nx*ny]-phase[xyz], 1.0f) - Numerics.modulo(phase[xyz]-phase[xyz+nx*ny], 1.0f);
		        double dphi1 = Numerics.modulo(phase[xyz-1-nx]-phase[xyz], 1.0f)     - Numerics.modulo(phase[xyz]-phase[xyz+1+nx], 1.0f);
		        double dphi2 = Numerics.modulo(phase[xyz-nx-nx*ny]-phase[xyz], 1.0f) - Numerics.modulo(phase[xyz]-phase[xyz+nx+nx*ny], 1.0f);
		        double dphi3 = Numerics.modulo(phase[xyz-nx*ny-1]-phase[xyz], 1.0f)  - Numerics.modulo(phase[xyz]-phase[xyz+nx*ny+1], 1.0f);
		        double dphi4 = Numerics.modulo(phase[xyz+1-nx]-phase[xyz], 1.0f)     - Numerics.modulo(phase[xyz]-phase[xyz-1+nx], 1.0f);
		        double dphi5 = Numerics.modulo(phase[xyz+nx-nx*ny]-phase[xyz], 1.0f) - Numerics.modulo(phase[xyz]-phase[xyz-nx+nx*ny], 1.0f);
		        double dphi6 = Numerics.modulo(phase[xyz+nx*ny-1]-phase[xyz], 1.0f)  - Numerics.modulo(phase[xyz]-phase[xyz-nx*ny+1], 1.0f);
		        double dphi7 = Numerics.modulo(phase[xyz-1-nx-nx*ny]-phase[xyz], 1.0f) - Numerics.modulo(phase[xyz]-phase[xyz+1+nx+nx*ny], 1.0f);
		        double dphi8 = Numerics.modulo(phase[xyz+1-nx-nx*ny]-phase[xyz], 1.0f) - Numerics.modulo(phase[xyz]-phase[xyz-1+nx+nx*ny], 1.0f);
		        double dphi9 = Numerics.modulo(phase[xyz-1+nx-nx*ny]-phase[xyz], 1.0f) - Numerics.modulo(phase[xyz]-phase[xyz+1-nx+nx*ny], 1.0f);
		        double dphi0 = Numerics.modulo(phase[xyz-1-nx+nx*ny]-phase[xyz], 1.0f) - Numerics.modulo(phase[xyz]-phase[xyz+1+nx-nx*ny], 1.0f);
		        double diff = FastMath.sqrt(dphix*dphix+dphiy*dphiy+dphiz*dphiz+dphi1*dphi1+dphi2*dphi2+dphi3*dphi3+dphi4*dphi4+dphi5*dphi5+dphi6*dphi6+dphi7*dphi7+dphi8*dphi8+dphi9*dphi9+dphi0*dphi0);
		        
		        if (diff>0) reliability[xyz] = 1.0f/(float)diff;
		        else reliability[xyz] = 0.0f;
		        
		        Rmean += reliability[xyz];
		        Rnum ++;
		    }
		}
		if (Rnum>0) Rmean /= Rnum;
		if (debug) BasicInfo.displayMessage("mean reliability: "+Rmean+"\n");		
        
		
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        // add to the heap
		heap.reset();
		heap.setMaxTree();
		for (int xyz=0;xyz<nxyz;xyz++) {
		    wrapcount[xyz] = 0;
		    label[xyz] = 0;
		}
        /*
		// init: any value above a certain intensity threshold?
		boolean[] init = ObjectExtraction.objectFromImage(reliability, nx,ny,nz, 2.0f*(float)Rmean, ObjectExtraction.SUPERIOR);
		// erode to keep only nicely connected regions
		//initmag = Morphology.erodeObject(initmag, nx,ny,nz, 2,2,2);
        int[] lbls = ObjectLabeling.connected6Object3D(init, nx, ny, nz);
        int nlb = ObjectLabeling.countLabels(lbls, nx,ny,nz);
        */
        // better: just build quadrants
        int nlb = nquadrant*nquadrant*nquadrant;
        int[] seeds = new int[nlb];
        for (int x=0; x<nx/nquadrant; x++) for (int y=0; y<ny/nquadrant; y++) for (int z = 0; z<nz/nquadrant; z++) {
            for (int dx=0;dx<nquadrant;dx++) for (int dy=0;dy<nquadrant;dy++) for (int dz=0;dz<nquadrant;dz++) {
                int xyz = x+dx*Numerics.floor(nx/nquadrant)+nx*(y+dy*Numerics.floor(ny/nquadrant))+nx*ny*(z+dz*Numerics.floor(nz/nquadrant));
                int qid = dx+dy*nquadrant+dz*nquadrant*nquadrant;
                if (reliability[xyz]>reliability[seeds[qid]]) {
                    seeds[qid] = xyz;
                }
            }
        }
        for (int l=0;l<nlb;l++) {
            int xyz = seeds[l];
		    for (byte k = 0; k<6; k++) if (!boundaries[k][xyz]) {
				//int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
				// all the neighbors
				if (mask[xyzn]) {
				    heap.addValue(reliability[xyz]+reliability[xyzn],xyz,xyzn);
				}
			}
		}
		if (debug) BasicInfo.displayMessage("number of seeds: "+nlb+"\n");	
		/*
        int seed = 0;
        for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
            if (reliability[xyz]>reliability[seed]) {
                seed = xyz;
            }
        }
        for (byte k = 0; k<6; k++) if (!boundaries[k][seed]) {
            int xyzn = Ngb.neighborIndex(k, seed, nx,ny,nz);
            heap.addValue(reliability[seed]+reliability[xyzn],seed,xyzn);
        }
		/*
        for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
		    for (byte k = 0; k<6; k++) if (!boundaries[k][xyz]) {
				//int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
				// all the neighbors
				if (mask[xyzn]) {
				    heap.addValue(reliability[xyz]+reliability[xyzn],xyz,xyzn);
				}
			}
		}*/
		int nlabels=0;
		int nmerge=0;
		while ( heap.isNotEmpty() ) {
        	//System.out.print(".");
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId1();
        	int xyzn = heap.getFirstId2();
			heap.removeFirst();

			// check if the region has been processed and labels are the same			    
			//if (label[xyz]!=0 && label[xyzn]!=0)			    
			if (label[xyz]!=0 && label[xyz]==label[xyzn])			    
			    continue;

            int offset = Numerics.wrap(phase[xyz]-phase[xyzn], 1.0f);
			    
			// check: both unlabeled, one unlabeled, all labeled
			if (label[xyz]==0 && label[xyzn]==0) {
			    //System.out.print("\n"+lbmag[xyz]+":");
                // new label
			    nlabels++;
                label[xyz] = nlabels;
			    
			    label[xyzn] = nlabels;
			    phase[xyzn] += offset;
                wrapcount[xyzn] += offset;
            } else if (label[xyz]==0) {
                //System.out.print("-"+label[xyz]+"-");
                // update only the empty one
			    label[xyz] = label[xyzn];
			    wrapcount[xyz] = wrapcount[xyzn];
			    phase[xyz] -= offset;
                wrapcount[xyz] -= offset;
            } else if (label[xyzn]==0) {
                //System.out.print("="+label[xyz]+"=");
                // update only the empty one
			    label[xyzn] = label[xyz];
			    wrapcount[xyzn] = wrapcount[xyz];
			    phase[xyzn] += offset;
                wrapcount[xyzn] += offset;
            } else {
                nmerge++;
                // both have non-zero labels: merge (keep the first component, arbitrarily)
                //System.out.print(label[xyz]+" <- "+label[xyzn]+": "+phase[xyz]+", "+phase[xyzn]+" = "+offset+"\n");
                int prev = label[xyzn];
                for (int xyzl=0;xyzl<nxyz;xyzl++) if (label[xyzl]==prev) {
                    label[xyzl] = label[xyz];
                    phase[xyzl] += offset;
                    wrapcount[xyzl] += offset;
                }
            }
            // find neighbors
            for (byte k = 0; k<6; k++) if (!boundaries[k][xyz]) {
				int xyzl = Ngb.neighborIndex(k, xyz, nx,ny,nz);
				if (mask[xyzl] && label[xyzl]==0) {
				   heap.addValue(reliability[xyz]+reliability[xyzl],xyz,xyzl);
				}
			}
			for (byte k = 0; k<6; k++) if (!boundaries[k][xyzn]) {
				int xyzl = Ngb.neighborIndex(k, xyzn, nx,ny,nz);
				//if (mask[xyzl] && label[xyzl]!=label[xyzn]) {
				if (mask[xyzl] && label[xyzl]==0) {
				    heap.addValue(reliability[xyzn]+reliability[xyzl],xyzn,xyzl);
				}
			}
		}
		if (debug) BasicInfo.displayMessage("number of labels: "+nlabels+"\n");		
		if (debug) BasicInfo.displayMessage("number of merges: "+nmerge+"\n");	
		
		// for testing: add the TV step                            
		if (!postprocess.equals("none")) {
            TotalVariation1D algo = new TotalVariation1D(phase,null,nx,ny,nz, tvscale, 0.125f, 0.00001f, 500);
            algo.setScaling(1.0f);
            //algo.solveWrapped();
            //if (postprocess.equals("TV-approximation")) correctImage = algo.exportResultWrapped();
            //else if (postprocess.equals("TV-residuals")) correctImage = algo.exportResidualWrapped();
            algo.solve();
            if (postprocess.equals("TV-approximation")) correctImage = algo.exportResult();
            else if (postprocess.equals("TV-residuals")) correctImage = algo.exportResidual();
            
        } else {
            correctImage = phase;
        }
        // rescale into original values
        for (int xyz=0;xyz<nxyz;xyz++) correctImage[xyz] = Pmin + correctImage[xyz]*(Pmax-Pmin);
        
		// output
		scoreImage = reliability;
		countImage = wrapcount;
		
		return;
    }
    
}

