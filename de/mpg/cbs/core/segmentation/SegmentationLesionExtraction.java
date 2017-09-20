package de.mpg.cbs.core.segmentation;


import java.net.URL;
import java.util.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.special.Erf;


/*
 * @author Pierre-Louis Bazin
 */
public class SegmentationLesionExtraction {

	// jist containers
	private float[] probaImage;
	private int[] segImage;
	private float[] mgdmImage;
	private float[] priorImage = null;
	
	private int nx, ny, nz, nxyz, nc;
	private float rx, ry, rz;

	private String atlasParam;
	private String regionParam;
	
	private float	gmdistanceParam = 1.0f;	
	private float	csfdistanceParam = 2.0f;	
	private float	lesdistanceParam = 2.0f;	
	private float	minprobaThreshold = 0.84f;	// set as a distance to the estimated GM cluster, should be fixed
	private float	maxprobaThreshold = 0.97f;	// to remove dirty WM, otherwise use the other one
	private float	minSize = 27.0f; // to remove small lesions
	
	private float[] pvRegionImage;
	private float[] lesionsizeImage;
	private float[] lesionprobaImage;
	private float[] boundarypvImage;
	private float[] scoreImage;
	private int[] labelImage;
		
	
	// input parameters
	public final void setSegmentationImage(int[] val) { segImage = val; }
	public final void setLevelsetBoundaryImage(float[] val) { mgdmImage = val; }
	public final void setProbaImage(float[] val) { probaImage = val; }
	public final void setLocationPriorImage(float[] val) { priorImage = val; }
		
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setComponents(int c) { nc=c; }
	
	public final void setAtlasFile(String val) { atlasParam = val; }
	
	public final void setGMPartialVolumingDistance(float val) { gmdistanceParam = val; }
	public final void setCSFPartialVolumingDistance(float val) { csfdistanceParam = val; }
	public final void setLesionClusteringDistance(float val) { lesdistanceParam = val; }
	
	public final void setMinProbabilityThreshold(float val) { minprobaThreshold = val; }
	public final void setMaxProbabilityThreshold(float val) { maxprobaThreshold = val; }
	public final void setMinimumSize(float val) { minSize = val; }
	
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Segmentation.devel"; }
	public final String getLabel() {return "Extract Lesions"; }
	public final String getName() {return "ExtractLesions"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Netherlands Institute for Neuroscience | Spinoza Centre for Neuroimaging"; }
	public final String getDescription() { return "Extracts lesions from a probability image and a pre-segmentation with MGDM."; }
		
	public final String getVersion() { return "3.1.1"; }

	// output parameters
	public final float[] getRegionPrior() { return pvRegionImage; }
	public final float[] getLesionSize() { return lesionsizeImage; }
	public final float[] getLesionProba() { return lesionprobaImage; }
	public final float[] getBoundaryPartialVolume() { return boundarypvImage; }
	public final int[] getLesionLabels() { return labelImage; }
	public final float[] getLesionScore() { return scoreImage; }
	
	
	public final void execute(){
				
		// load mask and build boolean signature for each region
		BasicInfo.displayMessage("Load atlas\n");
	
		SimpleShapeAtlas2 atlas = new SimpleShapeAtlas2(atlasParam);

		int maxlb = 0;
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getLabels()[nobj]>maxlb) maxlb = atlas.getLabels()[nobj];
		}
		// define CRWM region, inter-ventricular region
		
		BitSet isRegion = new BitSet(maxlb);
		BitSet isVentL = new BitSet(maxlb);
		BitSet isVentR = new BitSet(maxlb);
		
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getNames()[nobj].equals("Cerebrum-WML")) isRegion.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("Cerebrum-WMR")) isRegion.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("CerebralWM")) isRegion.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("Cerebrum-WM")) isRegion.set(atlas.getLabels()[nobj]);
			// add the subcortex back? not very clean..
			//else if (atlas.getNames()[nobj].equals("CaudateL")) isRegion.set(atlas.getLabels()[nobj]);
			//else if (atlas.getNames()[nobj].equals("PutamenL")) isRegion.set(atlas.getLabels()[nobj]);
			//else if (atlas.getNames()[nobj].equals("ThalamusL")) isRegion.set(atlas.getLabels()[nobj]);
			//else if (atlas.getNames()[nobj].equals("GlobusPallidusL")) isRegion.set(atlas.getLabels()[nobj]);
			//else if (atlas.getNames()[nobj].equals("CaudateR")) isRegion.set(atlas.getLabels()[nobj]);
			//else if (atlas.getNames()[nobj].equals("PutamenR")) isRegion.set(atlas.getLabels()[nobj]);
			//else if (atlas.getNames()[nobj].equals("ThalamusR")) isRegion.set(atlas.getLabels()[nobj]);
			//else if (atlas.getNames()[nobj].equals("GlobusPallidusR")) isRegion.set(atlas.getLabels()[nobj]);
		}
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getNames()[nobj].equals("VentricleL")) isVentL.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("ChoroidPlexusL")) isVentL.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("Ventricle3")) isVentL.set(atlas.getLabels()[nobj]);	
		}
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getNames()[nobj].equals("VentricleR")) isVentR.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("ChoroidPlexusR")) isVentR.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("Ventricle3")) isVentR.set(atlas.getLabels()[nobj]);	
		}
		
		// signed distance functions
		System.out.println("region mask");
		float[] lvlreg = new float[nxyz];
		float[] lvlvenL = new float[nxyz];
		float[] lvlvenR = new float[nxyz];
		// build levelset boundaries
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (mgdmImage[xyz]==-1) lvlreg[xyz] = gmdistanceParam+1.0f;
			else if (isRegion.get(segImage[xyz])) lvlreg[xyz] = -mgdmImage[xyz];
			else if (priorImage!=null && priorImage[xyz]>0.5f) lvlreg[xyz] = 1.0f - 2.0f*priorImage[xyz];
			else if (priorImage!=null) lvlreg[xyz] = Numerics.max(0.001f, Numerics.min(mgdmImage[xyz], 1.0f - 2.0f*priorImage[xyz]));
			else lvlreg[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
			
			if (mgdmImage[xyz]==-1) lvlvenL[xyz] = csfdistanceParam+1.0f;
			else if (isVentL.get(segImage[xyz])) lvlvenL[xyz] = -mgdmImage[xyz];
			else lvlvenL[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
			
			if (mgdmImage[xyz]==-1) lvlvenR[xyz] = csfdistanceParam+1.0f;
			else if (isVentR.get(segImage[xyz])) lvlvenR[xyz] = -mgdmImage[xyz];
			else lvlvenR[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
		}
		lvlreg = ObjectTransforms.fastMarchingDistanceFunction(lvlreg,nx,ny,nz);
		lvlvenL =  ObjectTransforms.fastMarchingDistanceFunction(lvlvenL,nx,ny,nz);
		lvlvenR =  ObjectTransforms.fastMarchingDistanceFunction(lvlvenR,nx,ny,nz);
		
		// PV for the WM/Gm and WM/Ventricle interface (ventricles have more PV than GM, typically)
		pvRegionImage = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			pvRegionImage[xyz] = Numerics.min(Numerics.bounded(0.0f - lvlreg[xyz]/gmdistanceParam, 0.0f, 1.0f),
											  Numerics.bounded(0.0f + lvlvenL[xyz]/csfdistanceParam, 0.0f, 1.0f),
											  Numerics.bounded(0.0f + lvlvenR[xyz]/csfdistanceParam, 0.0f, 1.0f) );
		}
		
		// PV for inter ventricle region
		float[] pvVentImage = new float[nxyz];
		float maxdist = 0.0f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (-lvlvenL[xyz]>maxdist) maxdist = -lvlvenL[xyz];
			if (-lvlvenR[xyz]>maxdist) maxdist = -lvlvenR[xyz];
		}
		for (int xyz=0;xyz<nxyz;xyz++) {
			// if closer than 2*maxdist to both ventricles, inside region; else out
			if (lvlvenL[xyz]>0 && lvlvenL[xyz]<2*maxdist && lvlvenR[xyz]>0 && lvlvenR[xyz]<2*maxdist)
				pvVentImage[xyz] = Numerics.bounded(1.0f - Numerics.max(lvlvenL[xyz],lvlvenR[xyz])/(2.0f*maxdist), 0.0f, 1.0f);
			else
				pvVentImage[xyz] = 0.0f;
		}
		// merge both regions
		for (int xyz=0;xyz<nxyz;xyz++) {
			pvRegionImage[xyz] = Numerics.min(pvRegionImage[xyz], 1.0f-pvVentImage[xyz]);
		}
		
		// crop proba image to only include regions with positive WM proba
		// combine raw proba with pv estimates?
		float minpv = 0.5f;
		float maxpv = 1.0f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (pvRegionImage[xyz]<minpv) probaImage[xyz] = 0.0f;
			else probaImage[xyz] = Numerics.min((pvRegionImage[xyz]-minpv)/(maxpv-minpv), 1.0f)*probaImage[xyz];
		}
				
		// step 1: find clusters of high probability
		boolean[] clusters = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (probaImage[xyz]>=minprobaThreshold) clusters[xyz] = true;
		}
		/*
		// basic approach: just count the number of connected components
		int[] lesionlabels = ObjectLabeling.connected6Object3D(clusters, nx, ny, nz);
		int nlabels = ObjectLabeling.countLabels(lesionlabels, nx,ny,nz);
		*/
		// other idea: find all local maxima
		// pre-process with gaussian kernel for scale, reduce/avoid exact value match
		float[] smoothproba = ImageFilters.separableConvolution(probaImage, nx,ny,nz, ImageFilters.separableGaussianKernel(lesdistanceParam,lesdistanceParam,lesdistanceParam));
		boolean[] maxima = new boolean[nxyz];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (clusters[xyz]) {
				boolean ismax = false;
				// select on original values, not smoothed (likely lowered)
				if (probaImage[xyz]>=maxprobaThreshold) ismax = true;
				for (byte d=0;d<26 && ismax;d++) {
					int ngb = Ngb.neighborIndex(d,xyz, nx, ny, nz);
					if (smoothproba[ngb]>=smoothproba[xyz]) ismax = false; 
				}
				maxima[xyz] = ismax;
			}
		}
		// then grow from closest proba competitively?
		int[] lesionlabels = growMaximaRegions(maxima, probaImage, minprobaThreshold);
		int nlabels = ObjectLabeling.countLabels(lesionlabels, nx,ny,nz);
		
			
		// step 2: measure size, average proba, partial voluming, etc
		float[] lesionsize = new float[nlabels];
		float[] lesionproba = new float[nlabels];
		float[] boundarypv = new float[nlabels];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (lesionlabels[xyz]>0) {
				int lb = lesionlabels[xyz]-1;
				lesionsize[lb] ++;
				lesionproba[lb] += (probaImage[xyz]-minprobaThreshold)/(1.0f-minprobaThreshold);
				boundarypv[lb] += pvRegionImage[xyz];
			}
		}
		
		
		// step 3. keep plausible lesions, discard the rest
		// if p(inter-ventricle)>0: discard
		for (int lb=0;lb<nlabels;lb++) {
			float pwm = boundarypv[lb]/lesionsize[lb];
			float psize = 1.0f - (float)FastMath.exp(-lesionsize[lb]/minSize);
			float ples = lesionproba[lb]/lesionsize[lb];
			boundarypv[lb] = pwm;
			lesionsize[lb] = psize;
			lesionproba[lb] = ples;
		}
			
		// step 4. build interesting maps
		lesionsizeImage = new float[nxyz];
		lesionprobaImage = new float[nxyz];
		boundarypvImage = new float[nxyz];
		labelImage = new int[nxyz];
		scoreImage = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (lesionlabels[xyz]>0) {
				int lb = lesionlabels[xyz]-1;
				lesionsizeImage[xyz] = lesionsize[lb];
				lesionprobaImage[xyz] = lesionproba[lb];
				boundarypvImage[xyz] = boundarypv[lb];
				if (maxima[xyz]) labelImage[xyz] = -2;
				else labelImage[xyz] = lesionlabels[xyz];
				//finalImage[xyz] = ventriclepv[lb]*pvRegionImage[xyz]*(probaImage[xyz]-minprobaThreshold)/(1.0f-minprobaThreshold);
				//finalImage[xyz] = ventriclepv[lb]*pvRegionImage[xyz]*lesionproba[lb]/lesionsize[lb]*(probaImage[xyz]-minprobaThreshold)/(1.0f-minprobaThreshold);
				//finalImage[xyz] = (float)FastMath.sqrt(lesionproba[lb]*lesionsize[lb]*(probaImage[xyz]-minprobaThreshold)/(1.0f-minprobaThreshold));
				scoreImage[xyz] = lesionproba[lb]*lesionsize[lb];
			} else if (clusters[xyz]) {
				lesionprobaImage[xyz] = (probaImage[xyz]-minprobaThreshold)/(1.0f-minprobaThreshold);
				boundarypvImage[xyz] = pvRegionImage[xyz];
				if (maxima[xyz]) labelImage[xyz] = -2;
				else labelImage[xyz] = -1;
			}
		}		
		return;
	}
	
	private final int[] growMaximaRegions(boolean[] maxima, float[] proba, float minp) {
		
		int[] regions = new int[nxyz];
		BinaryHeapPair heap = new BinaryHeapPair(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MAXTREE);
		heap.reset();
		// initialize the heap from maxima
		int nmax=0;       
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
        	int xyz = x+nx*y+nx*ny*z;
        	if (maxima[xyz]) {
        		nmax++;
        		regions[xyz] = nmax;
				// search for boundaries
				for (byte k = 0; k<6; k++) {
					int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
					if (!maxima[xyzn] && proba[xyzn]>minp) {
						// add to the heap with previous value
						//heap.addValue(proba[xyzn]*Ngb.neighborDistanceRatio(k),xyzn,nmax);
						//heap.addValue(proba[xyzn],xyzn,nmax);
						heap.addValue(proba[xyz]*proba[xyzn],xyzn,nmax);
					}
				}
            } else {
            	regions[xyz] = -1;
            }
        }
        // grow all at once until empty or unreachable
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId1();
        	int lb = heap.getFirstId2();
			heap.removeFirst();

			if (regions[xyz]>0)  continue;
			
			// update the distance functions at the current level
			regions[xyz] = lb; // update the current level
 			
			// find new neighbors
			for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				
				// must be in outside the object or its processed neighborhood
				if (regions[xyzn]==-1 && proba[xyzn]>minp) {
					// add to the heap (6-connected go first)
					//heap.addValue(proba[xyzn]*Ngb.neighborDistanceRatio(k),xyzn,lb);
					//heap.addValue(proba[xyzn],xyzn,lb);
					heap.addValue(dist*proba[xyzn],xyzn,lb);
				}
			}
		}
		for (int xyz=0;xyz<nxyz;xyz++) if (regions[xyz]==-1) regions[xyz]=0;
		return regions;
	}
	
}
