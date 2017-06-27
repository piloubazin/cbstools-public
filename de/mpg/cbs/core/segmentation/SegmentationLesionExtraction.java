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
	
	private int nx, ny, nz, nxyz, nc;
	private float rx, ry, rz;

	private String atlasParam;
	private String regionParam;
	public static final String[] regionTypes = {"crwm", "csf", "cbwm"};
	private String backgroundParam;
	public static final String[] backgroundTypes = {"crgm", "brain", "cbgm"};
	private float	distanceParam = 1.0f;
		
	private float	probaThreshold = 0.667f;
	private float	minSize = 3.0f;
	
	private float[] pvRegionImage;
	private float[] pvVentImage;
	private float[] lesionsizeImage;
	private float[] lesionprobaImage;
	private float[] boundarypvImage;
	private float[] ventriclepvImage;
		
	
	// input parameters
	public final void setSegmentationImage(int[] val) { segImage = val; }
	public final void setLevelsetBoundaryImage(float[] val) { mgdmImage = val; }
	public final void setProbaImage(float[] val) { probaImage = val; }
		
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setComponents(int c) { nc=c; }
	
	public final void setAtlasFile(String val) { atlasParam = val; }
	
	public final void setEnhancedRegion(String val) { regionParam = val; }
	public final void setContrastBackground(String val) { backgroundParam = val; }
	public final void setPartialVolumingDistance(float val) { distanceParam = val; }
	
	public final void setProbabilityThreshold(float val) { probaThreshold = val; }
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
	public final float[] getInterVentriclePrior() { return pvVentImage; }
	public final float[] getLesionSize() { return lesionsizeImage; }
	public final float[] getLesionProba() { return lesionprobaImage; }
	public final float[] getBoundaryPartialVolume() { return boundarypvImage; }
	public final float[] getVentriclePartialVolume() { return ventriclepvImage; }
	
	
	public final void execute(){
				
		// load mask and build boolean signature for each region
		Interface.displayMessage("Load atlas\n");
	
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
		}
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getNames()[nobj].equals("VentricleL")) isVentL.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("ChoroidPlexusL")) isVentL.set(atlas.getLabels()[nobj]);
		}
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getNames()[nobj].equals("VentricleR")) isVentR.set(atlas.getLabels()[nobj]);
			else if (atlas.getNames()[nobj].equals("ChoroidPlexusR")) isVentR.set(atlas.getLabels()[nobj]);
		}
		
		// signed distance functions
		System.out.println("region mask");
		float[] lvlreg = new float[nxyz];
		float[] lvlvenL = new float[nxyz];
		float[] lvlvenR = new float[nxyz];
		// build levelset boundaries
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (mgdmImage[xyz]==-1) lvlreg[xyz] = distanceParam+1.0f;
			else if (isRegion.get(segImage[xyz])) lvlreg[xyz] = -mgdmImage[xyz];
			else lvlreg[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
			
			if (mgdmImage[xyz]==-1) lvlvenL[xyz] = distanceParam+1.0f;
			else if (isVentL.get(segImage[xyz])) lvlvenL[xyz] = -mgdmImage[xyz];
			else lvlvenL[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
			
			if (mgdmImage[xyz]==-1) lvlvenR[xyz] = distanceParam+1.0f;
			else if (isVentR.get(segImage[xyz])) lvlvenR[xyz] = -mgdmImage[xyz];
			else lvlvenR[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
		}
		lvlreg = ObjectTransforms.fastMarchingDistanceFunction(lvlreg,nx,ny,nz);
		lvlvenL =  ObjectTransforms.fastMarchingDistanceFunction(lvlvenL,nx,ny,nz);
		lvlvenR =  ObjectTransforms.fastMarchingDistanceFunction(lvlvenR,nx,ny,nz);
		
		// PV for the WM
		pvRegionImage = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			pvRegionImage[xyz] = Numerics.bounded(0.5f - 0.5f*lvlreg[xyz]/distanceParam, 0.0f, 1.0f);
		}
		
		// PV for inter ventricle region
		pvVentImage = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			pvVentImage[xyz] = Numerics.bounded(1.0f - 0.25f*lvlvenL[xyz]/distanceParam - 0.25f*lvlvenR[xyz]/distanceParam, 0.0f, 1.0f);
		}
		
		// step 1: find clusters of high probability
		boolean[] clusters = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (probaImage[xyz]>probaThreshold && isRegion.get(segImage[xyz])) clusters[xyz] = true;
		}
		int[] lesionlabels = ObjectLabeling.connected6Object3D(clusters, nx, ny, nz);
		int nlabels = ObjectLabeling.countLabels(lesionlabels, nx,ny,nz);
		
		// step 2: measure size, average proba, partial voluming, etc
		float[] lesionsize = new float[nlabels];
		float[] lesionproba = new float[nlabels];
		float[] boundarypv = new float[nlabels];
		float[] ventriclepv = new float[nlabels];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (lesionlabels[xyz]>0) {
				int lb = lesionlabels[xyz]-1;
				lesionsize[lb] ++;
				lesionproba[lb] += probaImage[xyz];
				boundarypv[lb] += pvRegionImage[xyz];
				ventriclepv[lb] += pvVentImage[xyz];
			}
		}
		
		// step 3. build interesting maps
		lesionsizeImage = new float[nxyz];
		lesionprobaImage = new float[nxyz];
		boundarypvImage = new float[nxyz];
		ventriclepvImage = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (lesionlabels[xyz]>0) {
				int lb = lesionlabels[xyz]-1;
				lesionsizeImage[xyz] = lesionsize[lb];
				lesionprobaImage[xyz] = lesionproba[lb];
				boundarypvImage[xyz] = boundarypv[lb];
				ventriclepvImage[xyz] = ventriclepv[lb];
			}	
		}

		return;
	}
	
}
