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
public class BrainDefineMultiRegionPriors {

	// jist containers
	private int[] segImage;
	private float[] mgdmImage;
	
	private int nx, ny, nz, nxyz, nc;
	private float rx, ry, rz;

	private String atlasParam;
	//private String regionParam;
	//public static final String[] regionTypes = {"inter-ventricular", "ventricle-horns", "internal-capsule"};
	//private String methodParam;
	//public static final String[] methodTypes = {"closest-distance", "equi-distance", "distance-size"};
	private float	distanceParam = 0.0f;
		
	private float[] interImage;
	private float[] hornsImage;
	private float[] icapsImage;
		
	//private String regionName;
		
	
	// input parameters
	public final void setSegmentationImage(int[] val) { segImage = val; }
	public final void setLevelsetBoundaryImage(float[] val) { mgdmImage = val; }
		
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setAtlasFile(String val) { atlasParam = val; }
	
	//public final void setDefinedRegion(String val) { regionParam = val; }
	//public final void setDefinitionMethod(String val) { methodParam = val; }
	public final void setDistanceOffset(float val) { distanceParam = val; }
		
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Brain Processing.devel"; }
	public final String getLabel() {return "Define Multi-Region Priors"; }
	public final String getName() {return "DefineMultiRegionPriors"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Netherlands Institute for Neuroscience | Spinoza Centre for Neuroimaging"; }
	public final String getDescription() { return "Defines location priors based on combinations of regions from a MGDM brain segmentation."; }
		
	public final String getVersion() { return "3.1.1"; }

	// output parameters
	public final float[] getInterVentricularPV() { return interImage; }
	public final float[] getVentricularHornsPV() { return hornsImage; }
	public final float[] getInternalCapsulePV() { return icapsImage; }
	
	//public final String getRegionName() { return regionName; }
	
	private static final byte X=0;
	private static final byte Y=1;
	private static final byte Z=2;
	private static final byte N=3;
	
	public final void execute(){
				
		// load mask and build boolean signature for each region
		BasicInfo.displayMessage("Load atlas\n");
	
		SimpleShapeAtlas2 atlas = new SimpleShapeAtlas2(atlasParam);

		int maxlb = 0;
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getLabels()[nobj]>maxlb) maxlb = atlas.getLabels()[nobj];
		}
		
		/// Define the different regions
		interImage = defineInterVentricularRegion(atlas,maxlb);
		hornsImage = defineVentricularHorns(atlas, maxlb);
		icapsImage = defineInternalCapsule(atlas, maxlb);
				
		return;
	}

	
	float[] defineInterVentricularRegion(SimpleShapeAtlas2 atlas, int maxlb) {
	    BitSet isVentL = new BitSet(maxlb);
        BitSet isVentR = new BitSet(maxlb);
        
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
        float[] lvlvenL = new float[nxyz];
        float[] lvlvenR = new float[nxyz];
        // build levelset boundaries
        for (int xyz=0;xyz<nxyz;xyz++) {
            if (mgdmImage[xyz]==-1) lvlvenL[xyz] = distanceParam+1.0f;
            else if (isVentL.get(segImage[xyz])) lvlvenL[xyz] = -mgdmImage[xyz];
            else lvlvenL[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
            
            if (mgdmImage[xyz]==-1) lvlvenR[xyz] = distanceParam+1.0f;
            else if (isVentR.get(segImage[xyz])) lvlvenR[xyz] = -mgdmImage[xyz];
            else lvlvenR[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
        }
        lvlvenL =  ObjectTransforms.fastMarchingDistanceFunction(lvlvenL,nx,ny,nz);
        lvlvenR =  ObjectTransforms.fastMarchingDistanceFunction(lvlvenR,nx,ny,nz);
        
        // PV for inter ventricle region
        float[] pvImage = new float[nxyz];
        float maxdist = 0.0f;
        for (int xyz=0;xyz<nxyz;xyz++) {
            if (-lvlvenL[xyz]>maxdist) maxdist = -lvlvenL[xyz];
            if (-lvlvenR[xyz]>maxdist) maxdist = -lvlvenR[xyz];
        }
        for (int xyz=0;xyz<nxyz;xyz++) {
            // if closer than 2*maxdist to both ventricles, inside region; else out
            if (lvlvenL[xyz]<2*maxdist && lvlvenR[xyz]<2*maxdist)
                pvImage[xyz] = Numerics.bounded(1.0f - Numerics.max(lvlvenL[xyz],lvlvenR[xyz])/(2.0f*maxdist), 0.0f, 1.0f);
            else
                pvImage[xyz] = 0.0f;
        }
        return pvImage;
    }
    
    float[] defineVentricularHorns(SimpleShapeAtlas2 atlas, int maxlb) {
        BitSet isVent = new BitSet(maxlb);
        BitSet isCtx = new BitSet(maxlb);
        BitSet isPut = new BitSet(maxlb);
        BitSet isCau = new BitSet(maxlb);
        BitSet isTha = new BitSet(maxlb);
        
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("VentricleL")) isVent.set(atlas.getLabels()[nobj]);
            else if (atlas.getNames()[nobj].equals("ChoroidPlexusL")) isVent.set(atlas.getLabels()[nobj]);
            else if (atlas.getNames()[nobj].equals("Ventricle3")) isVent.set(atlas.getLabels()[nobj]);	
            else if (atlas.getNames()[nobj].equals("VentricleR")) isVent.set(atlas.getLabels()[nobj]);
            else if (atlas.getNames()[nobj].equals("ChoroidPlexusR")) isVent.set(atlas.getLabels()[nobj]);
            else if (atlas.getNames()[nobj].equals("Ventricle3")) isVent.set(atlas.getLabels()[nobj]);	
        }
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("Cerebrum-GML")) isCtx.set(atlas.getLabels()[nobj]);
            else if (atlas.getNames()[nobj].equals("Cerebrum-GMR")) isCtx.set(atlas.getLabels()[nobj]);
            else if (atlas.getNames()[nobj].equals("Cerebrum-GM")) isCtx.set(atlas.getLabels()[nobj]);	
            else if (atlas.getNames()[nobj].equals("CerebralGM")) isCtx.set(atlas.getLabels()[nobj]);	
        }
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("PutamenR")) isPut.set(atlas.getLabels()[nobj]);
            if (atlas.getNames()[nobj].equals("PutamenL")) isPut.set(atlas.getLabels()[nobj]);
        }
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("CaudateR")) isCau.set(atlas.getLabels()[nobj]);
            if (atlas.getNames()[nobj].equals("CaudateL")) isCau.set(atlas.getLabels()[nobj]);
        }
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("ThalamusR")) isTha.set(atlas.getLabels()[nobj]);
            if (atlas.getNames()[nobj].equals("ThalamusL")) isTha.set(atlas.getLabels()[nobj]);
        }
        
        // signed distance functions
        System.out.println("region mask");
        float[] lvlven = new float[nxyz];
        float[] lvlctx = new float[nxyz];
        float[] lvlput = new float[nxyz];
        float[] lvlcau = new float[nxyz];
        float[] lvltha= new float[nxyz];
        // build levelset boundaries
        for (int xyz=0;xyz<nxyz;xyz++) {
            if (mgdmImage[xyz]==-1) lvlven[xyz] = distanceParam+1.0f;
            else if (isVent.get(segImage[xyz])) lvlven[xyz] = -mgdmImage[xyz];
            else lvlven[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
            
            if (mgdmImage[xyz]==-1) lvlctx[xyz] = distanceParam+1.0f;
            else if (isCtx.get(segImage[xyz])) lvlctx[xyz] = -mgdmImage[xyz];
            else lvlctx[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
            
            if (mgdmImage[xyz]==-1) lvlput[xyz] = distanceParam+1.0f;
            else if (isPut.get(segImage[xyz])) lvlput[xyz] = -mgdmImage[xyz];
            else lvlput[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);

            if (mgdmImage[xyz]==-1) lvlcau[xyz] = distanceParam+1.0f;
            else if (isCau.get(segImage[xyz])) lvlcau[xyz] = -mgdmImage[xyz];
            else lvlcau[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
            
            if (mgdmImage[xyz]==-1) lvltha[xyz] = distanceParam+1.0f;
            else if (isTha.get(segImage[xyz])) lvltha[xyz] = -mgdmImage[xyz];
            else lvltha[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
        }
        lvlven =  ObjectTransforms.fastMarchingDistanceFunction(lvlven,nx,ny,nz);
        lvlctx =  ObjectTransforms.fastMarchingDistanceFunction(lvlctx,nx,ny,nz);
        lvlput =  ObjectTransforms.fastMarchingDistanceFunction(lvlput,nx,ny,nz);
        lvlcau =  ObjectTransforms.fastMarchingDistanceFunction(lvlcau,nx,ny,nz);
        lvltha =  ObjectTransforms.fastMarchingDistanceFunction(lvltha,nx,ny,nz);
        
        // PV for inter ventricle region
        float[] pvImage = new float[nxyz];
        float maxdist = 0.0f;
        for (int xyz=0;xyz<nxyz;xyz++) {
            //if (-lvlven[xyz]>maxdist) maxdist = -lvlven[xyz];
            if (-lvlcau[xyz]>maxdist) maxdist = -lvlcau[xyz];
        }
        /*
        // ventricle vs. cortex : only targeting the corpus callosum region (and stuff lower down)
        for (int xyz=0;xyz<nxyz;xyz++) {
            // if closer than 2*maxdist to ventricles and cortex, inside region; else out
            if ( ( (lvlvenL[xyz]>0 && lvlvenL[xyz]<2*maxdist) || (lvlvenR[xyz]>0 && lvlvenR[xyz]<2*maxdist) ) && lvlctx[xyz]>0 && lvlctx[xyz]<2*maxdist)
                //pvRegionImage[xyz] = Numerics.bounded(1.0f - Numerics.min(lvlvenL[xyz],lvlvenR[xyz])/(2.0f*maxdist), 0.0f, 1.0f);
                pvRegionImage[xyz] = Numerics.max( (1.0f-lvlvenL[xyz]/(2.0f*maxdist)), (1.0f-lvlvenR[xyz]/(2.0f*maxdist)))*(1.0f-lvlctx[xyz]/(2.0f*maxdist));
            else
                pvRegionImage[xyz] = 0.0f;
        }
        // putamen vs cortex: not bad for CC, but specifically avoids the horns
        for (int xyz=0;xyz<nxyz;xyz++) {
            // if closer than 2*maxdist to ventricles and cortex, inside region; else out
            if (lvlctx[xyz]>0 && lvlput[xyz]>0 && lvltha[xyz]>0 && lvlput[xyz]>lvlctx[xyz] && lvltha[xyz]>lvlctx[xyz] && lvlcau[xyz]<2.0f*maxdist)
                //pvRegionImage[xyz] = Numerics.bounded(1.0f - Numerics.min(lvlvenL[xyz],lvlvenR[xyz])/(2.0f*maxdist), 0.0f, 1.0f);
                pvRegionImage[xyz] = (Numerics.min(lvlput[xyz],lvltha[xyz])-lvlctx[xyz])/(2.0f*maxdist);
            else
                pvRegionImage[xyz] = 0.0f;
        }*/
        float[] dcau = new float[4];
        float[] dput = new float[4];
        float[] dven = new float[4];
        for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
            int xyz = x+nx*y+nx*ny*z;
            dcau[X] = lvlcau[xyz+1]-lvlcau[xyz-1];
            dcau[Y] = lvlcau[xyz+nx]-lvlcau[xyz-nx];
            dcau[Z] = lvlcau[xyz+nx*ny]-lvlcau[xyz-nx*ny];
            dcau[N] = (float)FastMath.sqrt(dcau[X]*dcau[X]+dcau[Y]*dcau[Y]+dcau[Z]*dcau[Z]);
            
            dput[X] = lvlput[xyz+1]-lvlput[xyz-1];
            dput[Y] = lvlput[xyz+nx]-lvlput[xyz-nx];
            dput[Z] = lvlput[xyz+nx*ny]-lvlput[xyz-nx*ny];
            dput[N] = (float)FastMath.sqrt(dput[X]*dput[X]+dput[Y]*dput[Y]+dput[Z]*dput[Z]);
            
            dven[X] = lvlven[xyz+1]-lvlven[xyz-1];
            dven[Y] = lvlven[xyz+nx]-lvlven[xyz-nx];
            dven[Z] = lvlven[xyz+nx*ny]-lvlven[xyz-nx*ny];
            dven[N] = (float)FastMath.sqrt(dven[X]*dven[X]+dven[Y]*dven[Y]+dven[Z]*dven[Z]);
            
            float angle = 0.0f;
            if (dven[N]*dput[N]>0) angle = (dven[X]*dput[X]+dven[Y]*dput[Y]+dven[Z]*dput[Z])/(dven[N]*dput[N]);
            
            // if angle <0 : inside
            if (lvlctx[xyz]>0 && lvlven[xyz]>0 && lvlput[xyz]>0)
                pvImage[xyz] = 0.5f*(1.0f+angle);
            else
                pvImage[xyz] = 0.0f;
            
            // build proba, restrict to Caudate
            if (lvlcau[xyz]>0.0f) pvImage[xyz] = 0.0f;
        }
        return pvImage;
    }			
		
	float[] defineInternalCapsule(SimpleShapeAtlas2 atlas, int maxlb) {
        BitSet isCau = new BitSet(maxlb);
        BitSet isPut = new BitSet(maxlb);
        BitSet isTha = new BitSet(maxlb);
        BitSet isWM = new BitSet(maxlb);
        BitSet isGP = new BitSet(maxlb);
        
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("CaudateR")) isCau.set(atlas.getLabels()[nobj]);
            if (atlas.getNames()[nobj].equals("CaudateL")) isCau.set(atlas.getLabels()[nobj]);
        }
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("PutamenR")) isPut.set(atlas.getLabels()[nobj]);
            if (atlas.getNames()[nobj].equals("PutamenL")) isPut.set(atlas.getLabels()[nobj]);
        }
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("ThalamusR")) isTha.set(atlas.getLabels()[nobj]);
            if (atlas.getNames()[nobj].equals("ThalamusL")) isTha.set(atlas.getLabels()[nobj]);
        }
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("Cerebrum-WML")) isWM.set(atlas.getLabels()[nobj]);
            else if (atlas.getNames()[nobj].equals("Cerebrum-WMR")) isWM.set(atlas.getLabels()[nobj]);
            else if (atlas.getNames()[nobj].equals("Cerebrum-WM")) isWM.set(atlas.getLabels()[nobj]);	
            else if (atlas.getNames()[nobj].equals("CerebralWM")) isWM.set(atlas.getLabels()[nobj]);	
        }
        for (int nobj=0;nobj<atlas.getNumber();nobj++) {
            if (atlas.getNames()[nobj].equals("GlobusPallidusL")) isGP.set(atlas.getLabels()[nobj]);
            if (atlas.getNames()[nobj].equals("GlobusPallidusR")) isGP.set(atlas.getLabels()[nobj]);
        }
        
        // signed distance functions
        System.out.println("region mask");
        float[] lvlcau = new float[nxyz];
        float[] lvlput = new float[nxyz];
        float[] lvltha= new float[nxyz];
        // build levelset boundaries (exculde ventricle horns from caudate if already computed)
        for (int xyz=0;xyz<nxyz;xyz++) {
            if (mgdmImage[xyz]==-1) lvlcau[xyz] = distanceParam+1.0f;
            else if (isCau.get(segImage[xyz]) && hornsImage==null) lvlcau[xyz] = -mgdmImage[xyz];
            else if (isCau.get(segImage[xyz]) && hornsImage!=null && hornsImage[xyz]<0.5f) lvlcau[xyz] = Numerics.max(-mgdmImage[xyz],2.0f*(hornsImage[xyz]-0.5f));
            else lvlcau[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
            
            if (mgdmImage[xyz]==-1) lvlput[xyz] = distanceParam+1.0f;
            else if (isPut.get(segImage[xyz])) lvlput[xyz] = -mgdmImage[xyz];
            else lvlput[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
            
            if (mgdmImage[xyz]==-1) lvltha[xyz] = distanceParam+1.0f;
            else if (isTha.get(segImage[xyz])) lvltha[xyz] = -mgdmImage[xyz];
            else lvltha[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
        }
        lvlcau =  ObjectTransforms.fastMarchingDistanceFunction(lvlcau,nx,ny,nz);
        lvlput =  ObjectTransforms.fastMarchingDistanceFunction(lvlput,nx,ny,nz);
        lvltha =  ObjectTransforms.fastMarchingDistanceFunction(lvltha,nx,ny,nz);
        
        // PV for internal capsule: function of distance gradients
        float[] pvImage = new float[nxyz];
        float[] lvlic = new float[nxyz];
        for (int xyz=0;xyz<nxyz;xyz++) lvlic[xyz] = distanceParam+1;
            
        float[] dcau = new float[4];
        float[] dput = new float[4];
        float[] dtha = new float[4];
        for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
            int xyz = x+nx*y+nx*ny*z;
            dcau[X] = lvlcau[xyz+1]-lvlcau[xyz-1];
            dcau[Y] = lvlcau[xyz+nx]-lvlcau[xyz-nx];
            dcau[Z] = lvlcau[xyz+nx*ny]-lvlcau[xyz-nx*ny];
            dcau[N] = (float)FastMath.sqrt(dcau[X]*dcau[X]+dcau[Y]*dcau[Y]+dcau[Z]*dcau[Z]);
            
            dput[X] = lvlput[xyz+1]-lvlput[xyz-1];
            dput[Y] = lvlput[xyz+nx]-lvlput[xyz-nx];
            dput[Z] = lvlput[xyz+nx*ny]-lvlput[xyz-nx*ny];
            dput[N] = (float)FastMath.sqrt(dput[X]*dput[X]+dput[Y]*dput[Y]+dput[Z]*dput[Z]);
            
            dtha[X] = lvltha[xyz+1]-lvltha[xyz-1];
            dtha[Y] = lvltha[xyz+nx]-lvltha[xyz-nx];
            dtha[Z] = lvltha[xyz+nx*ny]-lvltha[xyz-nx*ny];
            dtha[N] = (float)FastMath.sqrt(dtha[X]*dtha[X]+dtha[Y]*dtha[Y]+dtha[Z]*dtha[Z]);
            
            // check for negative angles between two of the structures
            float angle = Numerics.min((dcau[X]*dput[X]+dcau[Y]*dput[Y]+dcau[Z]*dput[Z])/(dcau[N]*dput[N]), 
                                       (dput[X]*dtha[X]+dput[Y]*dtha[Y]+dput[Z]*dtha[Z])/(dput[N]*dtha[N]), 
                                       (dtha[X]*dcau[X]+dtha[Y]*dcau[Y]+dtha[Z]*dcau[Z])/(dtha[N]*dcau[N]));
            
            // if angle <0 : inside
            /*
            if (angle<0)
                pvRegionImage[xyz] = -angle;
            else
                pvRegionImage[xyz] = 0.0f;
            
            // build proba, restrict to WM, add Cau/Put/Tha/GP
            float pvdist = Numerics.min(mgdmImage[xyz]/distanceParam,1.0f);
            if (isCau.get(segImage[xyz])) pvRegionImage[xyz] = pvdist;
            else if (isPut.get(segImage[xyz])) pvRegionImage[xyz] = pvdist;
            else if (isTha.get(segImage[xyz])) pvRegionImage[xyz] = pvdist;
            else if (isGP.get(segImage[xyz])) pvRegionImage[xyz] = pvdist;
            else if (isWM.get(segImage[xyz])) pvRegionImage[xyz] *= pvdist;
            */
            // restrict to WM+GP?
            if (isGP.get(segImage[xyz]) || isWM.get(segImage[xyz]))
                pvImage[xyz] = 0.5f*(1.0f-angle);
            else 
                pvImage[xyz] = 0.0f;
            
            // only inside WM and GP + CAU,PUT,THA: rebuild levelset
            if (mgdmImage[xyz]==-1) lvlic[xyz] = distanceParam+1.0f;
            else if (isCau.get(segImage[xyz])) lvlic[xyz] = -mgdmImage[xyz];
            else if (isPut.get(segImage[xyz])) lvlic[xyz] = -mgdmImage[xyz];
            else if (isTha.get(segImage[xyz])) lvlic[xyz] = -mgdmImage[xyz];
            else if (isGP.get(segImage[xyz])) lvlic[xyz] = -mgdmImage[xyz];
            else if (isWM.get(segImage[xyz])) lvlic[xyz] = 0.5f-pvImage[xyz];
            else lvlic[xyz] = Numerics.max(0.001f, mgdmImage[xyz]);
        }
        lvlic =  ObjectTransforms.fastMarchingDistanceFunction(lvlic,nx,ny,nz);
        for (int xyz=0;xyz<nxyz;xyz++) {
            pvImage[xyz] = Numerics.bounded(0.5f-0.5f*lvlic[xyz]/distanceParam, 0.0f, 1.0f);
        }
        return pvImage;
    }

}
