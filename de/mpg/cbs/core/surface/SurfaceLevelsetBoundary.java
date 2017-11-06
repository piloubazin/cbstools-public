package de.mpg.cbs.core.surface;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class SurfaceLevelsetBoundary {

	// jist containers
	private float[] lvlImage;
	private float levelParam = 0.0f;
	private int[] boundaryImage;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	// create inputs
	public final void setLevelSetImage(float[] val) { lvlImage = val; }
	public final void setLevelSetBoundaryValue(float val) { levelParam = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Surfaces"; }
	public final String getLabel() { return "Level Set Boundary"; }
	public final String getName() { return "LevelSetBoundary"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Finds a one-voxel thick boundary for a level set function with arbitrary values."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.2"; };
	
	// create outputs
	public final int[] getBoundaryImage() { return boundaryImage; }
	
	public void execute() {
		
		boundaryImage = new int[nx*ny*nz];
		
		// find regions on the boundary
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			// 6C neighbors => 26C boundary
			for (byte n=0;n<6;n++) {
			    int ngb = Ngb.neighborIndex(n, xyz, nx, ny, nz);
			    if ( (lvlImage[xyz]>levelParam) != (lvlImage[ngb]>levelParam) ) {
			        // boundary between xyz and ngb
			        if (Numerics.abs(lvlImage[xyz]-levelParam)<Numerics.abs(lvlImage[ngb]-levelParam)) {
			            // only label if the smallest of the two?
			            boundaryImage[xyz] = 1;
			        }
			    }
			}
		}
	}

}
