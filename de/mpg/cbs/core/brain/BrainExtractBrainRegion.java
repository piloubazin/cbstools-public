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
public class BrainExtractBrainRegion {

	// jist containers
	private int[] segImage;
	private float[] mgdmImage;
	private int[] labelImage;
	private float[] functionImage;
	
	private int nx, ny, nz, nxyz, nc;
	private float rx, ry, rz;

	private String atlasParam;
	private String regionParam;
	private static final String[] regionTypes = {"left_cerebrum", "right_cerebrum", "cerebrum", "cerebellum", 
													"cerebellum_brainstem", "subcortex", "tissues(anat)", "tissues(func)",
													"brain_mask"};
	private boolean	normalizeParam = true;
	private boolean	densityParam = false;
	private float	scalingParam = 1.0f;
		
	private byte[] segStructureImage;
	private byte[] segInsideImage;
	private byte[] segBackgroundImage;
	private float[] lvlStructureImage;
	private float[] lvlInsideImage;
	private float[] lvlBackgroundImage;
	private float[] probaStructureImage;
	private float[] probaInsideImage;
	private float[] probaBackgroundImage;
		
	private String structureName;
	private String insideName;
	private String backgroundName;
		
	
	// input parameters
	public final void setSegmentationImage(int[] val) { segImage = val; }
	public final void setLevelsetBoundaryImage(float[] val) { mgdmImage = val; }
	public final void setMaximumMembershipImage(float[] val) { functionImage = val; }
	public final void setMaximumLabelImage(int[] val) { labelImage = val; }
		
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setComponents(int c) { nc=c; }
	
	public final void setAtlasFile(String val) { atlasParam = val; }
	
	public final void setExtractedRegion(String val) {
	    regionParam = val; 
	    // set up the structure names here in case we need to skip computations
	    if (regionParam.equals("left_cerebrum")) {
            structureName = "lcrgm";
            insideName = "lcrwm";
            backgroundName = "lcrbg";
		} else
		if (regionParam.equals("right_cerebrum")) {
            structureName = "rcrgm";
            insideName = "rcrwm";
            backgroundName = "rcrbg";
		} else
		if (regionParam.equals("cerebrum")) {
            structureName = "crgm";
            insideName = "crwm";
            backgroundName = "crbg";
		} else
		if (regionParam.equals("cerebellum")) {
            structureName = "cbgm";
            insideName = "cbwm";
            backgroundName = "cbbg";
		} else
		if (regionParam.equals("cerebellum_brainstem")) {
            structureName = "cbsgm";
            insideName = "cbswm";
            backgroundName = "cbsbg";
		} else
		if (regionParam.equals("subcortex")) {
            structureName = "subgm";
            insideName = "subwmbg";
            backgroundName = "subcsf";
		} else
		if (regionParam.equals("tissues(anat)")) {
            structureName = "angm";
            insideName = "anwm";
            backgroundName = "ancsf";
		} else
		if (regionParam.equals("tissues(func)")) {
            structureName = "fngm";
            insideName = "fnwm";
            backgroundName = "fncsf";
		} else
		if (regionParam.equals("brain_mask")) {
            structureName = "csf";
            insideName = "brain";
            backgroundName = "bg";
		}	    
	}
	public final void setNormalizeProbabilities(boolean val) { normalizeParam = val; }
	public final void setEstimateTissueDensities(boolean val) { densityParam = val; }
	public final void setPartialVolumingDistance(float val) { scalingParam = val; }
		
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Brain Processing"; }
	public final String getLabel() {return "Extract Brain Region"; }
	public final String getName() {return "ExtractBrainRegion"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Extract a selected region of interest from a MGDM brain segmentation."; }
		
	public final String getVersion() { return "3.1"; }

	// output parameters
	public final byte[] getInsideWMmask() { return segInsideImage; }
	public final byte[] getStructureGMmask() { return segStructureImage; }
	public final byte[] getBackgroundCSFmask() { return segBackgroundImage; }
	
	public final float[] getInsideWMlevelset() { return lvlInsideImage; }
	public final float[] getStructureGMlevelset() { return lvlStructureImage; }
	public final float[] getBackgroundCSFlevelset() { return lvlBackgroundImage; }
	
	public final float[] getInsideWMprobability() { return probaInsideImage; }
	public final float[] getStructureGMprobability() { return probaStructureImage; }
	public final float[] getBackgroundCSFprobability() { return probaBackgroundImage; }
	
	public final String getInsideName() { return insideName; }
	public final String getStructureName() { return structureName; }
	public final String getBackgroundName() { return backgroundName; }
	
	public final void execute(){
				
		// load mask and build boolean signature for each region
		BasicInfo.displayMessage("Load atlas\n");
	
		SimpleShapeAtlas2 atlas = new SimpleShapeAtlas2(atlasParam);

		int maxlb = 0;
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getLabels()[nobj]>maxlb) maxlb = atlas.getLabels()[nobj];
		}
		
		System.out.println("Extracting region: "+regionParam);
		
		BitSet isStructure = new BitSet(maxlb);
		BitSet isInside = new BitSet(maxlb);
		BitSet isBackground = new BitSet(maxlb);
		structureName = null; 
		insideName = null;
		backgroundName = null;
		if (regionParam.equals("left_cerebrum")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebrum-GML")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Cerebrum-WML")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CaudateL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusL")) isInside.set(atlas.getLabels()[nobj]);
			
				else isBackground.set(atlas.getLabels()[nobj]);
			}
            structureName = "lcrgm";
            insideName = "lcrwm";
            backgroundName = "lcrbg";
		} else
		if (regionParam.equals("right_cerebrum")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebrum-GMR")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Cerebrum-WMR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CaudateR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusR")) isInside.set(atlas.getLabels()[nobj]);
			
				else isBackground.set(atlas.getLabels()[nobj]);
			}
            structureName = "rcrgm";
            insideName = "rcrwm";
            backgroundName = "rcrbg";
		} else
		if (regionParam.equals("cerebrum")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebrum-GML")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-GMR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-GM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CerebralGM")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Cerebrum-WML")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CaudateL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusL")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-WMR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CaudateR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CerebralWM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-WM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricles")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("SubcorticalGM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("SubCorticalGM")) isInside.set(atlas.getLabels()[nobj]);
				
				else isBackground.set(atlas.getLabels()[nobj]);
			}
            structureName = "crgm";
            insideName = "crwm";
            backgroundName = "crbg";
		} else
		if (regionParam.equals("cerebellum")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebellum-GM")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isInside.set(atlas.getLabels()[nobj]);
				//else if (atlas.getNames()[nobj].equals("Ventricle4")) isInside.set(atlas.getLabels()[nobj]);
				
				else isBackground.set(atlas.getLabels()[nobj]);
			}
            structureName = "cbgm";
            insideName = "cbwm";
            backgroundName = "cbbg";
		} else
		if (regionParam.equals("cerebellum_brainstem")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebellum-GM")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Brainstem")) isInside.set(atlas.getLabels()[nobj]);
				//else if (atlas.getNames()[nobj].equals("Ventricle4")) isInside.set(atlas.getLabels()[nobj]);
				
				else isBackground.set(atlas.getLabels()[nobj]);
			}
            structureName = "cbsgm";
            insideName = "cbswm";
            backgroundName = "cbsbg";
		} else
		if (regionParam.equals("subcortex")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("VentricleL")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusL")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleR")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusR")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricle3")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricles")) isBackground.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("CaudateL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("StriatumL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusL")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("CaudateR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("StriatumR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusR")) isStructure.set(atlas.getLabels()[nobj]);
				
				else isInside.set(atlas.getLabels()[nobj]);
			}
            structureName = "subgm";
            insideName = "subwmbg";
            backgroundName = "subcsf";
		} else
		if (regionParam.equals("tissues(anat)")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebrum-GML")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-GMR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-GM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CerebralGM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebellum-GM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CaudateL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CaudateR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("SubcorticalGM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("SubCorticalGM")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Cerebrum-WML")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-WMR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CerebralWM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-WM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Brainstem")) isInside.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Sulcal-CSF")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleL")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusL")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleR")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusR")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricle3")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricle4")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricles")) isBackground.set(atlas.getLabels()[nobj]);
			}
            structureName = "angm";
            insideName = "anwm";
            backgroundName = "ancsf";
		} else
		if (regionParam.equals("tissues(func)")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Cerebrum-GML")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-GMR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-GM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CerebralGM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebellum-GM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CaudateL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusL")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("AmygdalaR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CaudateR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("HippocampusR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("PutamenR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ThalamusR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("GlobusPallidusR")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("SubcorticalGM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("SubCorticalGM")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Brainstem")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Cerebrum-WML")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-WMR")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("CerebralWM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebrum-WM")) isInside.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isInside.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Sulcal-CSF")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleL")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusL")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("VentricleR")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("ChoroidPlexusR")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricle3")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricle4")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Ventricles")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Arteries")) isBackground.set(atlas.getLabels()[nobj]);
			}
            structureName = "fngm";
            insideName = "fnwm";
            backgroundName = "fncsf";
		} else
		if (regionParam.equals("brain_mask")) {
			for (int nobj=0;nobj<atlas.getNumber();nobj++) {
				if (atlas.getNames()[nobj].equals("Sulcal-CSF")) isStructure.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Arteries")) isStructure.set(atlas.getLabels()[nobj]);
				
				else if (atlas.getNames()[nobj].equals("Background")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Sinuses")) isBackground.set(atlas.getLabels()[nobj]);
				else if (atlas.getNames()[nobj].equals("Dura")) isBackground.set(atlas.getLabels()[nobj]);
				
				else isInside.set(atlas.getLabels()[nobj]);
			}
            structureName = "csf";
            insideName = "brain";
            backgroundName = "bg";
		}
		//System.out.println("(output extensions: "+structureName+", "+insideName+", "+backgroundName+")");
		
		// simplified? only use first gdm for approximate scoring
		// not good enough at further boundaries
		MgdmRepresentation mgdmfull = null;
		int nmgdm = nc;
			
		if (densityParam) {
			System.out.println("Estimating densities: pre-processing");
			byte[] objlabel = atlas.getLabels();
			int nobj = objlabel.length;
			// too large: really slow...
			//float dist = 2*scaling+10.0f/rx;
			float dist = 2.0f*scalingParam+1.0f;
			/* needed?
			byte[] tmp = new byte[nobj-1];
			for (int n=0;n<nobj-1;n++) tmp[n] = objlabelImage[n+1];
			objlabel = tmp;
			nobj--;
			*/
			mgdmfull = new MgdmRepresentation(segImage, mgdmImage, nx,ny,nz, rx,ry,rz,
												objlabel, nobj, nmgdm, false, dist);
		}		
		
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
				System.out.println("re-estimate densities");
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
				System.out.println("normalize probabilities");
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
