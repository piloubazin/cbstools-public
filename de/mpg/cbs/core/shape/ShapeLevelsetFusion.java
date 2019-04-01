package de.mpg.cbs.core.shape;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

//import Jama.*;

/*
 * @author Pierre-Louis Bazin
 */
public class ShapeLevelsetFusion {

	// data buffers
	private float[][] lvlImages;
	
	private int nsubj;
	//private boolean skelParam;
	private boolean topoParam = true;
	private     String	            lutdir = null;
	
	private float[] 			meanImage;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	public final void setNumberOfImages(int val) {
	    nsubj = val;
	    lvlImages = new float[nsubj][];
	}
	public final void setLevelsetImageAt(int n, float[] val) { lvlImages[n] = val; }
	
	//public static final void setFollowSkeleton(boolean val) { skelParam=val; }
	public final void setCorrectSkeletonTopology(boolean val) { topoParam=val; }
	public final void setTopologyLUTdirectory(String val) { lutdir = val; }

	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Shape.devel"; }
	public final String getLabel() { return "Levelset Fusion"; }
	public final String getName() { return "LevelsetFusion"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Combines a collection of levelset surfaces into an average with average volume and spherical topology"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.3"; };

	// create outputs
	public final float[] getLevelsetAverage() { return meanImage; }

	public void execute() {
		
		float[][] levelsets = lvlImages;
		
		System.out.println("compute average");
		float[] average = null;
        average = new float[nxyz];
        for (int xyz=0;xyz<nxyz;xyz++) {	
            average[xyz] = 0.0f;
            for (int n=0;n<nsubj;n++) average[xyz] += levelsets[n][xyz]/(float)nsubj;
		}
		// optional: correct topology here
		if (topoParam==true) {
            ShapeTopologyCorrection2 topocorrect = new ShapeTopologyCorrection2();
		
            topocorrect.setDimensions(nx,ny,nz);
            topocorrect.setResolutions(rx,ry,rz);
		
            topocorrect.setShapeImage(average);
            topocorrect.setShapeImageType("signed_distance_function");
		
            topocorrect.setPropagationDirection("background->object");
            topocorrect.setTopology("26/6");
            topocorrect.setTopologyLUTdirectory(lutdir);
            topocorrect.setMinimumDistance(0.00001f);
		
            topocorrect.execute();
		
            average = topocorrect.getCorrectedImage();
		}
		
		// compute mean volume to normalize the average to the same size
		double meanvol = 0.0;
		float mindist = 1e6f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			for (int n=0;n<nsubj;n++) if (levelsets[n][xyz]<0) meanvol += 1.0/nsubj;
			if (average[xyz]<mindist) mindist = average[xyz];
		}
		System.out.println("mean volume (voxels): "+meanvol);
		System.out.println("minimum avg. distance: "+mindist);
		
		// find appropriate threshold to have correct volume; should use a fast marching approach!
		BinaryHeap2D	heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap4D.MINTREE);
		double vol = 0.0;
        boolean[] label = new boolean[nxyz];
        for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
            int xyz=x+nx*y+nx*ny*z;
            if (average[xyz]<0 || average[xyz]==mindist) {
                vol++;
                label[xyz] = true;
            
                // check for neighbors outside
                for (byte k = 0; k<6 ; k++) {
                    int ngb = Ngb.neighborIndex(k, xyz, nx, ny, nz);
                    if (average[ngb]>=0) {
                        // add to the heap (multiple instances are OK)
                        heap.addValue(average[ngb],ngb,(byte)1);
                    }
                }
            } else {
                label[xyz] = false;
            }
        }
        // run until the volume exceeds the mean volume
        float threshold = 0.0f;
        while (heap.isNotEmpty() && vol<meanvol) {
            threshold = heap.getFirst();
            int xyz = heap.getFirstId();
            heap.removeFirst();	
            if (label[xyz]==false) {
                vol++;
                label[xyz] = true;
                // add neighbors
                for (byte k = 0; k<6; k++) {
                    int ngb = Ngb.neighborIndex(k, xyz, nx, ny, nz);
                    if (label[ngb]==false) {
                        heap.addValue(average[ngb],ngb,(byte)1);
                    }
                }
            }
        }
        System.out.println("Distance offset: "+threshold);
        for (int xyz=0;xyz<nxyz;xyz++) {
            average[xyz] -= threshold;
        }
        // result
        meanImage = average;
        
        return;
	}


}
