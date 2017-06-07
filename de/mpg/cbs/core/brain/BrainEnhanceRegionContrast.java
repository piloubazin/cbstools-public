package de.mpg.cbs.core.brain;


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
public class BrainEnhanceRegionContrast {

	// jist containers
	private float[] intensImage;
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
		
	private byte[] regionImage;
	private byte[] backgroundImage;
	private float[] probaRegionImage;
	private float[] probaBackgroundImage;
	private float[] pvRegionImage;
	private float[] pvBackgroundImage;
		
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
	public final byte[] getBackgroundMask() { return backgroundImage; }
	public final float[] getRegionProbability() { return probaRegionImage; }
	public final float[] getBackgroundProbability() { return probaBackgroundImage; }
	public final float[] getRegionPartialVolume() { return pvRegionImage; }
	public final float[] getBackgroundPartialVolume() { return pvBackgroundImage; }
	
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
		regionName = null; 
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
		System.out.println("region mask");
		regionImage = new byte[nxyz];
		backgroundImage = new byte[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) if (segImage[xyz]>-1) {
			if (isRegion.get(segImage[xyz])) {
				regionImage[xyz] = 1;
			} 
			if (isBackground.get(segImage[xyz])) {
				backgroundImage[xyz] = 1;
			}
		}
		
		// 2. estimate pv factors
		System.out.println("region mask");
		pvRegionImage = new float[nxyz];
		pvBackgroundImage = new float[nxyz];
		// build levelset boundaries
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (regionImage[xyz]>0) pvRegionImage[xyz] = -mgdmImage[xyz];
			else pvRegionImage[xyz] = mgdmImage[xyz];
			
			if (backgroundImage[xyz]>0) pvBackgroundImage[xyz] = -mgdmImage[xyz];
			else pvBackgroundImage[xyz] = mgdmImage[xyz];
		}
		pvRegionImage = ObjectTransforms.fastMarchingDistanceFunction(pvRegionImage,nx,ny,nz);
		pvBackgroundImage =  ObjectTransforms.fastMarchingDistanceFunction(pvBackgroundImage,nx,ny,nz);
		// map back to linear PV values
		for (int xyz=0;xyz<nxyz;xyz++) {
			pvRegionImage[xyz] = Numerics.bounded(0.5f - 0.5f*pvRegionImage[xyz]/distanceParam, 0.0f, 1.0f);
			pvBackgroundImage[xyz] = Numerics.bounded(0.5f - 0.5f*pvBackgroundImage[xyz]/distanceParam, 0.0f, 1.0f);
		}
		
		// 3. build contrast maps
		int nreg = 0;
		for (int xyz=0;xyz<nxyz;xyz++) if (regionImage[xyz]>0) nreg++;
		int nbkg = 0;
		for (int xyz=0;xyz<nxyz;xyz++) if (backgroundImage[xyz]>0) nbkg++;
		
		double[] reg = new double[nreg];
		double[] bkg = new double[nbkg];
		double[] wreg = new double[nreg];
		double[] wbkg = new double[nbkg];
		int nr=0;
		int nb=0;
		for (int xyz=0;xyz<nxyz;xyz++) if (regionImage[xyz]>0) {
			reg[nr] = intensImage[xyz];
			wreg[nr] = pvRegionImage[xyz];
			nr++;
		} 
		for (int xyz=0;xyz<nxyz;xyz++) if (backgroundImage[xyz]>0) {
			bkg[nr] = intensImage[xyz];
			wbkg[nr] = pvBackgroundImage[xyz];
			nb++;
		}
		// robust statistics
		double dev = 100.0*0.5*Erf.erf(1.0/FastMath.sqrt(2.0));
		
		double regmed = ImageStatistics.weightedPercentile(reg, wreg, 50.0, nreg);
		double regstd = 0.5*(ImageStatistics.weightedPercentile(reg,wreg,50.0+dev,nreg) - ImageStatistics.weightedPercentile(reg,wreg,50.0-dev,nreg));
		Interface.displayMessage("Region parameter estimates: mean = "+regmed+", stdev = "+regstd+",\n");
		
		double bkgmed = ImageStatistics.weightedPercentile(bkg, wbkg, 50.0, nbkg);
		double bkgstd = 0.5*(ImageStatistics.weightedPercentile(bkg,wbkg,50.0+dev,nbkg) - ImageStatistics.weightedPercentile(bkg,wbkg,50.0-dev,nbkg));
		Interface.displayMessage("Background parameter estimates: mean = "+bkgmed+", stdev = "+bkgstd+",\n");
		
		// process the data: probabilities
		probaRegionImage = new float[nxyz];
		probaBackgroundImage = new float[nxyz];
		System.out.println("Probabilities");
		for (int xyz=0;xyz<nxyz;xyz++) {
			double preg = FastMath.exp(-0.5*(intensImage[xyz]-regmed)*(intensImage[xyz]-regmed)/regstd/regstd)/regstd;
			double pbkg = FastMath.exp(-0.5*(intensImage[xyz]-bkgmed)*(intensImage[xyz]-bkgmed)/bkgstd/bkgstd)/bkgstd;
		
			probaRegionImage[xyz] = (float)(preg/(preg+pbkg));
			probaBackgroundImage[xyz] = (float)(pbkg/(pbkg+preg));
		}

		return;
	}
	
}
