package de.mpg.cbs.core.surface;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class SurfaceProbabilityToLevelset {

	// jist containers
	private float[] lvlImage;
	private float scaleParam = 5.0f;
	private float[] probaImage;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	// create inputs
	public final void setProbabilityImage(float[] val) { probaImage = val; }
	public final void setScale_mm(float val) { scaleParam = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Surfaces"; }
	public final String getLabel() { return "Probability To Level Set"; }
	public final String getName() { return "ProbabilityToLevelSet"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Convert a probability map thresholded at 0.5 into a level set surface (incl. re-initialization)."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
	
	public final float[] getLevelSetImage() { return lvlImage; }
	
	public void execute() {
		
		float[] proba = probaImage;
		
		// transform the values
		float[] levelset = new float[nx*ny*nz];
		//float scale = scaleParam.getValue().floatValue();
		//System.out.println("transform data (scale: "+scale);
		
		// buidl levelset on boundary
		if (nz==1) {
            for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) {
                int xyz = x+nx*y;
                if (proba[xyz]>=0.5f && (proba[xyz+1]<0.5f || proba[xyz-1]<0.5f
                                           || proba[xyz+nx]<0.5f || proba[xyz-nx]<0.5f)) 
                    levelset[xyz] = 0.5f-proba[xyz];
                else if (proba[xyz]<0.5f && (proba[xyz+1]>=0.5f || proba[xyz-1]>=0.5f
                                               || proba[xyz+nx]>=0.5f || proba[xyz-nx]>=0.5f))
                    levelset[xyz] = 0.5f-proba[xyz];
                else if (proba[xyz]>=0.5f) levelset[xyz] = -1.0f;
                else levelset[xyz] = +1.0f;
            }		    
		} else {
            for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
                int xyz = x+nx*y+nx*ny*z;
                if (proba[xyz]>=0.5f && (proba[xyz+1]<0.5f || proba[xyz-1]<0.5f
                                           || proba[xyz+nx]<0.5f || proba[xyz-nx]<0.5f
                                           || proba[xyz+nx*ny]<0.5f || proba[xyz-nx*ny]<0.5f)) 
                    levelset[xyz] = 0.5f-proba[xyz];
                else if (proba[xyz]<0.5f && (proba[xyz+1]>=0.5f || proba[xyz-1]>=0.5f
                                               || proba[xyz+nx]>=0.5f || proba[xyz-nx]>=0.5f
                                               || proba[xyz+nx*ny]>=0.5f || proba[xyz-nx*ny]>=0.5f))
                    levelset[xyz] = 0.5f-proba[xyz];
                else if (proba[xyz]>=0.5f) levelset[xyz] = -1.0f;
                else levelset[xyz] = +1.0f;
            }
        }
		// no masking
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = true;
		}
		// expand boundary? yes!
		if (nz==1) {
            InflateGdm2D gdm = new InflateGdm2D(levelset, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no", null);
            gdm.evolveNarrowBand(0, 1.0f);
            lvlImage = gdm.getLevelSet();
		} else {
            InflateGdm gdm = new InflateGdm(levelset, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no", null);
            gdm.evolveNarrowBand(0, 1.0f);
            lvlImage = gdm.getLevelSet();
        }
		System.out.println("\n done");
	}


}
