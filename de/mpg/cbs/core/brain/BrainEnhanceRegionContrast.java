package de.mpg.cbs.core.brain;


import java.net.URL;
import java.util.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class BrainEnhanceRegionContrast {

	// jist containers
	private float[] intensImage;
	private int[] segImage;
	private float[] mgdmImage;
	
	private int nx, ny, nz, nxyz, nc;
	private float rx, ry, rz;

	private String atlasParam;
	private String regionParam;
	private static final String[] regionTypes = {"crwm", "csf", "cbwm"};
	private String backgroundParam;
	private static final String[] backgroundTypes = {"crgm", "brain", "cbgm"};
	private float	distanceParam = 1.0f;
	private boolean includeBg = false;
		
	private byte[] regionImage;
	private float[] probaRegionImage;
	private float[] probaBackgroundImage;
		
	private String regionName;
	private String backgroundName;
		
	
	// input parameters
	public final void setSegmentationImage(int[] val) { segImage = val; }
	public final void setLevelsetBoundaryImage(float[] val) { mgdmImage = val; }
	public final void setIntensityImage(float[] val) { intensImage = val; }
		
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setComponents(int c) { nc=c; }
	
	public final void setAtlasFile(String val) { atlasParam = val; }
	
	public final void setEnhancedRegion(String val) { regionParam = val; }
	public final void setContrastBackground(String val) { backgroundParam = val; }
	public final void setPartialVolumingDistance(float val) { distanceParam = val; }
	public final void setIncludeBackground(boolean val) { includeBg = val; }
		
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Brain Processing.devel"; }
	public final String getLabel() {return "Enhance Region Contrast"; }
	public final String getName() {return "EnhanceRegionContrast"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Netherlands Institute for Neuroscience | Spinoza Centre for Neuroimaging"; }
	public final String getDescription() { return "Enhances the contrast between selected regions from a MGDM brain segmentation."; }
		
	public final String getVersion() { return "3.1.1"; }

	// output parameters
	public final byte[] getRegionMask() { return regionImage; }
	public final float[] getRegionProbability() { return probaRegionImage; }
	public final float[] getBackgroundProbability() { return probaBackgroundImage; }
	
	public final String getRegionName() { return regionName; }
	public final String getBackgroundName() { return backgroundName; }
	
	public final void execute(){
				
		// load mask and build boolean signature for each region
		Interface.displayMessage("Load atlas\n");
	
		SimpleShapeAtlas2 atlas = new SimpleShapeAtlas2(atlasParam);

		int maxlb = 0;
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getLabels()[nobj]>maxlb) maxlb = atlas.getLabels()[nobj];
		}
		
		System.out.println("Extracting region: "+regionParam);
		
		BitSet isRegion = new BitSet(maxlb);
		BitSet isBackground = new BitSet(maxlb);
		structureName = null; 
		insideName = null;
		backgroundName = null;
		if (regionParam.equals("crwm")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebrum-WML")) isRegion.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-WMR")) isRegion.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CerebralWM")) isRegion.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-WM")) isRegion.set(atlas.getLabels()[nobj]);
				
				regionName = "_crwm";
			}
		}
		
		if (backgroundParam.equals("crgm")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebrum-GML")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-GMR")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-GM")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CerebralGM")) isBackground.set(atlas.getLabels()[nobj]);
				
				backgroundName = "_crgm";
			}
		}
		System.out.println("(output extensions: "+regionName+", "+backgroundName+")");
		
		// 1. define regions of interest
		
		
		// 2. estimate pv factors
		
		
		// 3. build contrast maps
		
		
		
		// process the data: probabilities
		probaStructureImage = new float[nx*ny*nz];
		probaInsideImage = new float[nx*ny*nz];
		probaBackgroundImage = new float[nx*ny*nz];
		System.out.println("Computing extracted probabilities");
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			probaStructureImage[xyz] = 0.0f;
			probaInsideImage[xyz] = 0.0f;
			probaBackgroundImage[xyz] = 0.0f;
			for (int c=0;c<nc;c++) {
				int xyzc = x+nx*y+nx*ny*z+nx*ny*nz*c;
				if (labelImage[xyzc]>-1) {
					if (isStructure.get(labelImage[xyzc])) {
						probaStructureImage[xyz] = Numerics.max(probaStructureImage[xyz], functionImage[xyzc]);
					} else if (isInside.get(labelImage[xyzc])) {
						probaInsideImage[xyz] = Numerics.max(probaInsideImage[xyz], functionImage[xyzc]);
					} else if (isBackground.get(labelImage[xyzc])) {
						probaBackgroundImage[xyz] = Numerics.max(probaBackgroundImage[xyz], functionImage[xyzc]);
					}
				}
			}
			if (densityParam) {
				// modulate by a sigmoid to lower or higher values away from the boundaries
				// find the closest label for each region
				boolean done = false;
				float factorS = -1;
				float factorI = -1;
				float factorB = -1;
				float dist = 0.0f;
				for (int n=0;n<nmgdm && !done;n++) {
					if (n==0) dist = -mgdmfull.getFunctions()[n][xyz];
					else if (n==1) dist *= -1;
					else dist += mgdmfull.getFunctions()[n-1][xyz];
					
					if (factorS==-1 && isStructure.get(atlas.getLabels()[mgdmfull.getLabels()[n][xyz]])) {
						factorS = (float)(2.0/(1.0+FastMath.exp(dist/scalingParam) ) );
					}
					if (factorI==-1 && isInside.get(atlas.getLabels()[mgdmfull.getLabels()[n][xyz]])) {
						factorI = (float)(2.0/(1.0+FastMath.exp(dist/scalingParam) ) );
					}
					if (factorB==-1 && isBackground.get(atlas.getLabels()[mgdmfull.getLabels()[n][xyz]])) {
						factorB = (float)(2.0/(1.0+FastMath.exp(dist/scalingParam) ) );
					}
					if (factorS!=-1 && factorI !=-1 && factorB!=-1) done = true;
				}
				if (factorS!=-1) probaStructureImage[xyz] *= factorS;
				if (factorI!=-1) probaInsideImage[xyz] *= factorI;
				if (factorB!=-1) probaBackgroundImage[xyz] *= factorB;
			}
			if (normalizeParam) {
				if (isStructure.get(segImage[xyz]) ) {
					float sum = Numerics.max(1e-3f,probaStructureImage[xyz]+probaInsideImage[xyz],probaStructureImage[xyz]+probaBackgroundImage[xyz]);
					probaStructureImage[xyz] /= sum;
					probaInsideImage[xyz] /= sum;
					probaBackgroundImage[xyz] /= sum;
				} else
				if (isInside.get(segImage[xyz]) ) {
					float sum = Numerics.max(1e-3f,probaInsideImage[xyz]+probaStructureImage[xyz],probaInsideImage[xyz]+probaBackgroundImage[xyz]);
					probaStructureImage[xyz] /= sum;
					probaInsideImage[xyz] /= sum;
					probaBackgroundImage[xyz] /= sum;
				} else
				if (isBackground.get(segImage[xyz]) ) {
					float sum = Numerics.max(1e-3f,probaBackgroundImage[xyz]+probaInsideImage[xyz],probaBackgroundImage[xyz]+probaStructureImage[xyz]);
					probaStructureImage[xyz] /= sum;
					probaInsideImage[xyz] /= sum;
					probaBackgroundImage[xyz] /= sum;
				} else {
					probaStructureImage[xyz] = 0.0f;
					probaInsideImage[xyz] = 0.0f;
					probaBackgroundImage[xyz] = 0.0f;
				}					
			}
			// remap everything into [0,1] in case it's not there
			probaStructureImage[xyz] = Numerics.bounded(probaStructureImage[xyz], 0.0f, 1.0f);
			probaInsideImage[xyz] = Numerics.bounded(probaInsideImage[xyz], 0.0f, 1.0f);
			probaBackgroundImage[xyz] = Numerics.bounded(probaBackgroundImage[xyz], 0.0f, 1.0f);
		}
		
		// process the data: segmentation
		System.out.println("Computing extracted segmentations");
		segStructureImage = new byte[nx*ny*nz];
		segInsideImage = new byte[nx*ny*nz];
		segBackgroundImage = new byte[nx*ny*nz];
		for (int xyz=0;xyz<nxyz;xyz++) if (segImage[xyz]>-1) {
			if (isStructure.get(segImage[xyz])) {
				segStructureImage[xyz] = 1;
			} else if (isInside.get(segImage[xyz])) {
				segInsideImage[xyz] = 1;
			} else if (isBackground.get(segImage[xyz])) {
				segBackgroundImage[xyz] = 1;
			}
		}

		// process the data: levelsets 
		// TODO: (could use the full representation and/or reinitialize => do it before the probabilities)
		System.out.println("Computing extracted level sets");
		lvlStructureImage = new float[nxyz];
		lvlInsideImage = new float[nxyz];
		lvlBackgroundImage = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) if (segImage[xyz]>-1) {
			if (isStructure.get(segImage[xyz])) {
				lvlStructureImage[xyz] = -mgdmImage[xyz];
				lvlInsideImage[xyz] = mgdmImage[xyz];
				lvlBackgroundImage[xyz] = mgdmImage[xyz];
			} else if (isInside.get(segImage[xyz])) {
				lvlStructureImage[xyz] = mgdmImage[xyz];
				lvlInsideImage[xyz] = -mgdmImage[xyz];
				lvlBackgroundImage[xyz] = mgdmImage[xyz];
			} else if (isBackground.get(segImage[xyz])) {
				lvlStructureImage[xyz] = mgdmImage[xyz];
				lvlInsideImage[xyz] = mgdmImage[xyz];
				lvlBackgroundImage[xyz] = -mgdmImage[xyz];
			}
		}

		return;
	}
	
}
