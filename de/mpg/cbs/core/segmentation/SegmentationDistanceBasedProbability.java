package de.mpg.cbs.core.segmentation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class SegmentationDistanceBasedProbability {

	// jist containers
	private float[] probaImage;
	//private float scaleParam = 5.0f;
	private float ratioParam = 0.5f;
	private int[] segImage;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;
	
	private int nlabels = 0;
	
	// create inputs
	public final void setSegmentationImage(int[] val) { segImage = val; }
	//public final void setScale_mm(float val) { scaleParam = val; }
	public final void setDistanceRatio(float val) { ratioParam = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Segmentation"; }
	public final String getLabel() { return "Distance-based Probability"; }
	public final String getName() { return "DistanceBasedProbability"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Convert a segmentation map into probabilities based on distance to the inside of each object"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
	
	public final float[] getProbabilityImage() { return probaImage; }
	public final int getLabelNumber() { return nlabels; }
	
	public void execute() {
		
		
		BasicInfo.displayMessage("build from max labeling and max probability...\n");
		int[] objlb = ObjectLabeling.listOrderedLabels(segImage,nx,ny,nz);
		nlabels = (byte)objlb.length;
			
		BasicInfo.displayMessage("found "+nlabels+" labels\n");
		// create a distance-based probability map
		float[] boundary = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) boundary[xyz] = 0.5f;
		
		byte[] lbs = new byte[nlabels];
		for (int l=0;l<nlabels;l++) lbs[l] = (byte)objlb[l];
		BasicInfo.displayMessage("distance-based MGDM representation...\n");
		MgdmRepresentation mgdm = new MgdmRepresentation(segImage, boundary, nx,ny,nz, rx,ry,rz, lbs, nlabels, 4, false, 9.0f);
			
		probaImage = new float[nlabels*nxyz];
		float[] maxobjdist = new float[nlabels];
		for (int xyz=0;xyz<nxyz;xyz++) {
			for (int l=0;l<nlabels;l++) if (segImage[xyz] == objlb[l]) {
				if (mgdm.getFunctions()[0][xyz]>maxobjdist[l]) maxobjdist[l] = mgdm.getFunctions()[0][xyz];
			}
		}
		
		// background: use the largest distance from foreground object instead as basis ? Or just 1-sum?
		for (int xyz=0;xyz<nxyz;xyz++) {
			float probaFg = 0.0f;
			for (byte l=1;l<nlabels;l++) {
				float dist = mgdm.reconstructedLevelSetAt(xyz, l);
				if (dist<-ratioParam*maxobjdist[l]) probaImage[l*nxyz+xyz] = 1.0f;
				else if (dist<0) probaImage[l*nxyz+xyz] = 0.5f - 0.5f*dist/maxobjdist[l]/ratioParam;
				else if (dist<ratioParam*maxobjdist[l]) probaImage[l*nxyz+xyz] = 0.5f - 0.5f*dist/maxobjdist[l]/ratioParam;
				else probaImage[l*nxyz+xyz] = 0.0f;
				
				probaFg += probaImage[l*nxyz+xyz];
			}	
			probaImage[xyz] = 1.0f-Numerics.min(probaFg, 1.0f);
		}
		
	}
	
}