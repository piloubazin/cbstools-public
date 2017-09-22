package de.mpg.cbs.methods;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;

import org.apache.commons.math3.util.FastMath;

/**
 * ACE : Anatomically Consistent Enhancement Algorithm (with convenient option wrappers)
 * 
 * @author Pierre-Louis Bazin
 * 
 */
public class AnatomicallyConsistentEnhancement {
	
	private byte[] wmmask;
	private float[] wmimg;
	private float[] gmimg;
	private float[] csfimg;

	private float[] phi;
	private float[] acegm;
	private float[] acecsf;
	private byte[] skel;
	private static	int 		nx,ny,nz, nxyz;   		// images dimensions
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;
	
	// computation variables for topology
	private boolean[][][] obj = new boolean[3][3][3];
	private CriticalPointLUT lut;
	private String	lutdir = null;
	
	private static final double		SQ3 = FastMath.sqrt(3)/2.0;
	
	/**
	 * constructor
	 */
	public AnatomicallyConsistentEnhancement(byte[] mask_, float[] wm_, float[] gm_, float[] csf_, int nx_, int ny_, int nz_,
											String connectivityType_, String lutdir_) {
	 
		wmmask = mask_;
		wmimg = wm_;
		gmimg = gm_;
		csfimg = csf_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		nxyz = nx*ny*nz;
		
		lutdir = lutdir_;
		
		phi = new float[nxyz];
		skel = new byte[nxyz];
		acegm = new float[nxyz];
		acecsf = new float[nxyz];
		
		lut = new CriticalPointLUT(lutdir,"critical266LUT.raw.gz",200);
		if (!lut.loadCompressedPattern()) {
			System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
			BasicInfo.displayMessage("Problem loading the algorithm's LUT from: "+lut.getFilename()+"\n");
			return;
		} else {
			//if (debug) System.out.println("LUT loaded from: "+lut.getFilename());
		}
		return;
	}
	
	public final float[] getACEprobaGM() { return acegm; }
	public final float[] getACEprobaCSF() { return acecsf; }
	public final byte[] getACEskeleton() { return skel; }
	
	public final void computeSulcalRidges(float ridgeThreshold) {
		
		// new, simplified version of the algorithm
		float[] cortex = new float[nxyz];
		float[] speed = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (wmmask[xyz]>0) phi[xyz] = 1.0f-2.0f*Numerics.max(0.51f,wmimg[xyz]);
			else phi[xyz] = 1.0f-2.0f*Numerics.min(0.49f,wmimg[xyz]);
			cortex[xyz] = 1.0f-2.0f*Numerics.max(wmmask[xyz],(wmimg[xyz]+gmimg[xyz]));
			//speed[xyz] = 1.0f-0.9f*csfimg[xyz];
			speed[xyz] = 1.0f-0.9f*(1.0f-wmimg[xyz]-gmimg[xyz]);
		}
		// note: a processing mask could be nice
		
		// 1. define domain: gm+wm+some distance (we grow a little more for boundaries)
		float maxdist = 3.0f;
		cortex = ObjectTransforms.fastMarchingDistanceFunction(cortex, maxdist+2.0f, nx,ny,nz);
		phi = ObjectTransforms.fastMarchingDistanceFunction(phi, maxdist+2.0f, nx,ny,nz);
		
		// 2. region of interest 
		boolean[] domain = new boolean[nxyz];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			domain[xyz] = (cortex[xyz]<=maxdist);
		}
		
		// 1b. smooth wm field? speed function?
		SmoothGdm smoothgdm = new SmoothGdm(phi, nx, ny, nz, 1.0f, 1.0f, 1.0f,
												domain, 0.4f, 0.1f, "no", null);
			
		smoothgdm.evolveNarrowBand(100, 0.001f);
		phi = smoothgdm.exportLevelset();
		smoothgdm.finalize();
		smoothgdm = null;
			
		// 2. Fast marching with CSF speeds
		phi = fastMarchingWeightedFunction(phi, speed, domain, nx,ny,nz);
		
		// 3. Compute central difference gradient magnitude
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			skel[xyz] = 0;
			acegm[xyz] = gmimg[xyz];
			acecsf[xyz] = csfimg[xyz];
			if (phi[xyz]>SQ3 && cortex[xyz]<=maxdist) {
				double dx = 0.5*(phi[xyz+1] - phi[xyz-1]);
				double dy = 0.5*(phi[xyz+nx] - phi[xyz-nx]);
				double dz = 0.5*(phi[xyz+nx*ny] - phi[xyz-nx*ny]);
				double mag = FastMath.sqrt(dx*dx+dy*dy+dz*dz);
				if (mag*speed[xyz]<ridgeThreshold) {
					skel[xyz] = 1;
					acegm[xyz] = (float)(mag*speed[xyz]*gmimg[xyz]);
					acecsf[xyz] += (gmimg[xyz]-acegm[xyz]);
				}
			}
		}
		return;
	}
	
	public final void thinSulcalSkeleton(boolean removeIsolated, int minimumSize) {
		// 4. thin the skeleton
		skel = medialThinning(skel, phi, nx, ny, nz);
		boolean[] skelb = new boolean[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			skelb[xyz] = false;
			if (skel[xyz]==0) {
				acegm[xyz] = gmimg[xyz];
				acecsf[xyz] = csfimg[xyz];
			} else skelb[xyz] = true;
		}
		
		if (removeIsolated) {
			// 5. look for small structures
			int[] lbl = ObjectLabeling.connected26Object3D(skelb, nx,ny,nz);
			int nlb = ObjectLabeling.countLabels(lbl, nx,ny,nz);
			int[] size = ObjectStatistics.volumes(lbl, nlb, nx,ny,nz);

			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (lbl[xyz]>0) {
					if (size[lbl[xyz]-1]<minimumSize) {
						acegm[xyz] = gmimg[xyz];
						acecsf[xyz] = csfimg[xyz];
						skel[xyz] = 0;
					}
				}
			}
		}
		return;
	}
		
	private final float[] fastMarchingWeightedFunction(float[] levelset, float[] speed, boolean[] mask, int nx, int ny, int nz)  {
        // computation variables
        boolean[] object = new boolean[nxyz]; // note: using a byte instead of boolean for the second pass
		boolean[] processed = new boolean[nxyz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		BinaryHeap2D heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
				        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        //if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
        // initialize mask and processing domain
        float maxlvl = Numerics.max(nx/2.0f,ny/2.0f,nz/2.0f);
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x+nx*y+nx*ny*z;
        	object[xyz] = (levelset[xyz]<=0);
			if (x==0 || x==nx-1 || y==0 || y==ny-1 || z==0 || z==nz-1) mask[x+nx*y+nx*ny*z] = false;
			if (!mask[xyz]) { // inside the masked region: either fully inside or fully outside
				if (object[xyz]) levelset[xyz] = -maxlvl;
				else levelset[xyz] = maxlvl;
			}
		}
		// initialize the heap from boundaries
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
        	int xyz = x+nx*y+nx*ny*z;
        	processed[xyz] = false;
        	// search for boundaries
        	for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				if (object[xyzn]!=object[xyz]) {
					// we assume the levelset value is correct at the boundary
					
					// add to the heap with previous value
					byte lb = 0;
					if (object[xyzn]) lb = 1;
					heap.addValue(Numerics.abs(levelset[xyzn]),xyzn,lb);
                }
            }
        }
		//if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        float maxdist = 0.0f;
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// update the distance functions at the current level
			if (lb==1) levelset[xyz] = -dist;
			else levelset[xyz] = dist;
			processed[xyz]=true; // update the current level
 			
			// find new neighbors
			for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn] && !processed[xyzn]) if (object[xyzn]==object[xyz]) {
					// compute new distance based on processed neighbors for the same object
					for (byte l=0; l<6; l++) {
						nbdist[l] = -1.0f;
						nbflag[l] = false;
						int xyznb = Ngb.neighborIndex(l, xyzn, nx, ny, nz);
						// note that there is at most one value used here
						if (mask[xyznb] && processed[xyznb]) if (object[xyznb]==object[xyz]) {
							nbdist[l] = Numerics.abs(levelset[xyznb]);
							nbflag[l] = true;
						}			
					}
					float newdist = weightedMarchingDistance(nbdist, nbflag, speed[xyzn]);
					
					// add to the heap
					heap.addValue(newdist,xyzn,lb);
				}
			}			
		}
		//if (debug) BasicInfo.displayMessage("done\n");		
		
       return levelset;
     }


	/**
     * the Fast marching distance computation 
     * (!assumes a 6D array with opposite coordinates stacked one after the other)
     *  and a strictly positive speed function
     * 
     */
    private final float weightedMarchingDistance(float[] val, boolean[] flag, float speed) {
   
        float s, s2; // s = a + b +c; s2 = a*a + b*b +c*c
        float tmp;
        int count;
        s = 0;
        s2 = 0;
        count = 0;

        for (int n=0; n<6; n+=2) {
			if (flag[n] && flag[n+1]) {
				tmp = Numerics.min(val[n], val[n+1]); // Take the smaller one if both are processed
				s += tmp;
				s2 += tmp*tmp;
				count++;
			} else if (flag[n]) {
				s += val[n]; // Else, take the processed one
				s2 += val[n]*val[n];
				count++;
			} else if (flag[n+1]) {
				s += val[n+1];
				s2 += val[n+1]*val[n+1];
				count++;
			}
		}
         // count must be greater than zero since there must be at least one processed pt in the neighbors
        if (count==0) System.err.print("!");
        if (s*s-count*(s2-1.0/(speed*speed))<0) {
        	System.err.print(":");
        	tmp = 0;
        	for (int n=0;n<6;n++) if (flag[n]) tmp = Numerics.max(tmp,val[n]);
        	for (int n=0;n<6;n++) if (flag[n]) tmp = Numerics.min(tmp,val[n]);
        } else {
			tmp = (s + (float)FastMath.sqrt((double) (s*s-count*(s2-1.0/(speed*speed)))))/count;
		}
        // The larger root
        return tmp;
    }

    public final byte[] medialThinning(byte[] skeleton, float[] funct, int nx, int ny, int nz) {
		
    	lut = new CriticalPointLUT(lutdir,"critical266LUT.raw.gz",200);
		if (!lut.loadCompressedPattern()) {
			System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
			BasicInfo.displayMessage("Problem loading the algorithm's LUT from: "+lut.getFilename()+"\n");
			return null;
		} else {
			//if (debug) System.out.println("LUT loaded from: "+lut.getFilename());
		}
		
		// from the boundary, find the middle (approx)
		if (verbose) BasicInfo.displayMessage("fast marching init\n");		
		BinaryHeap3D heap = new BinaryHeap3D(nx*ny+ny*nz+nz*nx, BinaryHeap3D.MINTREE);
		
		heap.reset();
		boolean[][][] processed = new boolean[nx][ny][nz];
		boolean[][][] medial = new boolean[nx][ny][nz];
		byte[][][] skel = new byte[nx][ny][nz];
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) for (short z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (skeleton[xyz]>0) {
				skel[x][y][z] = 1;
				processed[x][y][z] = false;
				medial[x][y][z] = false;
					
				// search for boundaries
				boolean isboundary=false;
				for (int k = 0; k<6 && !isboundary; k++) {
					int ngb = x+Ngb.x[k] + nx*(y+Ngb.y[k]) + nx*ny*(z+Ngb.z[k]);
					if (skeleton[ngb]==0) {
						isboundary=true;
					}
				}
				// add to the heap
				if (isboundary) {
					heap.addValue(funct[xyz], x, y, z);
					if (debug) BasicInfo.displayMessage(".");
				}
			} else {
				skel[x][y][z] = 0;
			}
		}
		if (debug) BasicInfo.displayMessage("init\n");		
		
        // thin the object
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			float dist = heap.getFirst();
			int x = heap.getFirstX();
			int y = heap.getFirstY();
			int z = heap.getFirstZ();
			heap.removeFirst();

			// if already processed, pass
			if (processed[x][y][z])  continue;
			
			// check for simple point
			if (isSimplePoint(skel, x, y, z)) {
				if (isMedialEndpoint(skel, x, y, z)) {
					// label as medial surface endpoint
					medial[x][y][z] = true;
					if (debug) BasicInfo.displayMessage("-");
				} else {
					if (debug) BasicInfo.displayMessage("x");
					skel[x][y][z] = 0;
					
					// find new neghbors: 6C
					for (int k = 0; k<6; k++) {
						int ngb = x+Ngb.x[k] + nx*(y+Ngb.y[k]) + nx*ny*(z+Ngb.z[k]);
						if (skel[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]]==1 && !processed[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]]) {
							// add to the heap
							heap.addValue(funct[ngb], x+Ngb.x[k], y+Ngb.y[k], z+Ngb.z[k]);
						}
					}
				}
						
				// discard point
				processed[x][y][z]=true; 
			} else {
				if (debug) BasicInfo.displayMessage("!");	
			}
		}
		// label the medial endpoints
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) for (short z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			skeleton[xyz] = skel[x][y][z];
		}
		
		if (verbose) BasicInfo.displayMessage("done\n");		
		
		return skeleton;
     }

     private final boolean isMedialEndpoint(byte[][][] image, int x, int y, int z) {
		// we assume a binary image with values 0 outside and 1 inside
		for (byte d=0;d<13;d++) {
			float val = 0.0f;
			if (d==0) {		
				val=(image[x][y-1][z] 		+image[x][y+1][z]
					+image[x][y][z-1] 		+image[x][y][z+1]
					+image[x][y-1][z-1] 	+image[x][y-1][z+1]
					+image[x][y+1][z-1]		+image[x][y+1][z+1]);
			} else if (d==1) {
				val=(image[x-1][y][z]		+image[x+1][y][z]
					+image[x][y][z-1]		+image[x][y][z+1]
					+image[x-1][y][z-1]		+image[x-1][y][z+1]
					+image[x+1][y][z-1]		+image[x+1][y][z+1]);
			} else if (d==2) { 			
				val=(image[x-1][y][z]		+image[x+1][y][z]
					 +image[x][y-1][z] 		+image[x][y+1][z]
					 +image[x-1][y-1][z]	+image[x-1][y+1][z]
					 +image[x+1][y-1][z]	+image[x+1][y+1][z]);
			} else if (d==3) { 			
				val=(image[x-1][y+1][z]		+image[x+1][y-1][z]
					 +image[x][y][z-1]		+image[x-1][y+1][z-1]
					 +image[x+1][y-1][z-1] 	+image[x][y][z+1]
					 +image[x-1][y+1][z+1]	+image[x+1][y-1][z+1]);
			} else if (d==4) { 			
				val=(image[x][y+1][z-1]		+image[x][y-1][z+1]
					 +image[x-1][y][z]		+image[x-1][y+1][z-1]
					 +image[x-1][y-1][z+1]	+image[x+1][y][z]	
					 +image[x+1][y+1][z-1]	+image[x+1][y-1][z+1]);
			} else if (d==5) { 			
				val=(image[x+1][y][z-1]		+image[x-1][y][z+1]
					 +image[x][y-1][z]		+image[x+1][y-1][z-1]
					 +image[x-1][y-1][z+1]	+image[x][y+1][z]		
					 +image[x+1][y+1][z-1]	+image[x-1][y+1][z+1]);
			} else if (d==6) { 			
				val=(image[x-1][y-1][z]		+image[x+1][y+1][z]
					 +image[x][y][z-1]		+image[x-1][y-1][z-1]
					 +image[x+1][y+1][z-1]	+image[x][y][z+1]		
					 +image[x-1][y-1][z+1]	+image[x+1][y+1][z+1]);
			} else if (d==7) { 			
				val=(image[x][y-1][z-1]		+image[x][y+1][z+1]
					 +image[x-1][y][z]		+image[x-1][y-1][z-1]
					 +image[x-1][y+1][z+1]	+image[x+1][y][z]	
					 +image[x+1][y-1][z-1]	+image[x+1][y+1][z+1]);
			} else if (d==8) { 			
				val=(image[x-1][y][z-1]		+image[x+1][y][z+1]
					 +image[x][y-1][z]		+image[x-1][y-1][z-1]
					 +image[x+1][y-1][z+1]	+image[x][y+1][z]	
					 +image[x-1][y+1][z-1]	+image[x+1][y+1][z+1]);
			} else if (d==9) { 			
				val=(image[x-1][y-1][z+1]	+image[x-1][y+1][z-1]
					 +image[x+1][y-1][z-1]	+image[x+1][y+1][z-1]
					 +image[x-1][y+1][z+1]	+image[x+1][y-1][z+1]);
			} else if (d==10) { 			
				val=(image[x+1][y-1][z+1]	+image[x+1][y+1][z-1]
					 +image[x-1][y-1][z-1]	+image[x-1][y+1][z-1]
					 +image[x-1][y-1][z+1]	+image[x+1][y+1][z+1]);
			} else if (d==11) { 			
				val=(image[x-1][y+1][z+1]	+image[x+1][y+1][z-1]
					 +image[x-1][y-1][z-1]	+image[x+1][y-1][z-1]
					 +image[x-1][y-1][z+1]	+image[x+1][y+1][z+1]);
			} else if (d==12) { 			
				val=(image[x-1][y+1][z+1]	+image[x+1][y-1][z+1]
					 +image[x-1][y-1][z-1]	+image[x+1][y-1][z-1]
					 +image[x-1][y+1][z-1]	+image[x+1][y+1][z+1]);
			}
			if (val==1) return true;
		}
		return false;
	}
	
	private final boolean isSimplePoint(byte[][][] image, int x, int y, int z) {
    	// does it change the topology of the new object ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (image[x+i][y+j][z+l]==1) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = true;
		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			
		return true;
    }

}
