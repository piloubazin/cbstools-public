package de.mpg.cbs.core.registration;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class RegistrationSurfaceDataToGroupwiseTemplate {

	// jist containers
	private float[]		sourceContrastImage;
	private int[]		sourceMaskImage = null;
	private float[]		sourceLevelsetImage;
	private float[] 	sourceMappingImage = null;
	private float[]		templateLevelsetImage;
	private float[]		templateMappingImage = null;
	
	private String  	methodParam;
	public static final String[]	mappingTypes = {"projected","raw"};
	
	private String interpParam="NN";	
	public static final String[] interpTypes = {"NN","linear"};
		
	private float[] 	mappedDataImage;
	private int[]  		mappedMaskImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private int ntx, nty, ntz, ntxyz;
	private int nsx, nsy, nsz, nst, nsxyz;
	private int nox, noy, noz, noxyz;
	private float rtx, rty, rtz;
	private float rsx, rsy, rsz;
	private float rox, roy, roz;

	// create inputs
	public final void setSourceContrastImage(float[] val) { sourceContrastImage = val; }
	public final void setSourceMaskImage(int[] val) { sourceMaskImage = val; }
	public final void setSourceLevelsetImage(float[] val) { sourceLevelsetImage = val; }
	public final void setSourceMappingImage(float[] val) { sourceMappingImage = val; }
	public final void setTemplateLevelsetImage(float[] val) { templateLevelsetImage = val; }
	public final void setTemplateMappingImage(float[] val) { templateMappingImage = val; }

	public final void setSourceDimensions(int x, int y, int z, int t) { nsx=x; nsy=y; nsz=z; nst=t; nsxyz=nsx*nsy*nsz; }
	public final void setSourceDimensions(int[] dim) { nsx=dim[0]; nsy=dim[1]; nsz=dim[2]; nst=dim[3]; nsxyz=nsx*nsy*nsz; }
	
	public final void setSourceResolutions(float x, float y, float z) { rsx=x; rsy=y; rsz=z; }
	public final void setSourceResolutions(float[] res) { rsx=res[0]; rsy=res[1]; rsz=res[2]; }

	public final void setTemplateDimensions(int x, int y, int z) { ntx=x; nty=y; ntz=z; ntxyz=ntx*nty*ntz; }
	public final void setTemplateDimensions(int[] dim) { ntx=dim[0]; nty=dim[1]; ntz=dim[2]; ntxyz=ntx*nty*ntz; }
	
	public final void setTemplateResolutions(float x, float y, float z) { rtx=x; rty=y; rtz=z; }
	public final void setTemplateResolutions(float[] res) { rtx=res[0]; rty=res[1]; rtz=res[2]; }

	public final void setOutputDimensions(int x, int y, int z) { nox=x; noy=y; noz=z; noxyz=nox*noy*noz; }
	public final void setOutputDimensions(int[] dim) { nox=dim[0]; noy=dim[1]; noz=dim[2]; noxyz=nox*noy*noz; }
	
	public final void setOutputResolutions(float x, float y, float z) { rox=x; roy=y; roz=z; }
	public final void setOutputResolutions(float[] res) { rox=res[0]; roy=res[1]; roz=res[2]; }

	public final void setSurfaceMappingMethod(String val) { methodParam = val; }
	public final void setInterpolation(String val) { interpParam = val; }
	
	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Registration"; }
	public final String getLabel() { return "Surface Data to Groupwise Template"; }
	public final String getName() { return "SurfaceDataToGroupwiseTemplate"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for Neuroimaging, Netherlands Institute for Neuroscience, Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Sample subject data along a surface into a groupwiese template surface, possibly remapped to canonical space"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.2"; };

	// create outputs		
	public final float[] getMappedData() { return mappedDataImage; }
	public final int[] getmappedDataMask() { return mappedMaskImage; }
		
	public final void execute(){

		// create template mapping if it doesn't exist
		if (templateMappingImage==null) {
			templateMappingImage = new float[3*ntxyz];
			for (int x=0; x<ntx; x++) for (int y=0; y<nty; y++) for (int z = 0; z<ntz; z++) {
				int xyz = x + ntx*y + ntx*nty*z;
				templateMappingImage[xyz+X*ntxyz] = x;
				templateMappingImage[xyz+Y*ntxyz] = y;
				templateMappingImage[xyz+Z*ntxyz] = z;
			}
			nox = ntx; noy = nty; noz = ntz; noxyz = ntxyz;
		}

		// create source mapping if it doesn't exist
		if (sourceMappingImage==null) {
			sourceMappingImage = new float[3*ntxyz];
			for (int x=0; x<ntx; x++) for (int y=0; y<nty; y++) for (int z = 0; z<ntz; z++) {
				int xyz = x + ntx*y + ntx*nty*z;
				sourceMappingImage[xyz+X*ntxyz] = x;
				sourceMappingImage[xyz+Y*ntxyz] = y;
				sourceMappingImage[xyz+Z*ntxyz] = z;
			}
		}
		
		// create source mask for data=0 if it doesn't exist
		boolean[] sourceMask = new boolean[nsxyz];
		if (sourceMaskImage==null) {
			for (int xyz=0;xyz<nsxyz;xyz++) {
				sourceMask[xyz] = (sourceContrastImage[xyz]!=0);
			}
		} else {
			for (int xyz=0;xyz<nsxyz;xyz++) {
				sourceMask[xyz] = (sourceMaskImage[xyz]!=0);
			}
		}			
		
		// define sampling space: template boundary
		boolean[] boundary = new boolean[noxyz];
		for (int x=1;x<nox-1;x++) for (int y=1;y<noy-1;y++) for (int z=1;z<noz-1;z++) {
			int xyz = x+nox*y+nox*noy*z;
			float xt = templateMappingImage[xyz+X*noxyz];
			float yt = templateMappingImage[xyz+Y*noxyz];
			float zt = templateMappingImage[xyz+Z*noxyz];
			float levelxyz = ImageInterpolation.linearInterpolation(templateLevelsetImage, 1e13f, xt,yt,zt, ntx,nty,ntz);
			// 6C neighbors => 26C boundary
			for (byte n=0;n<6;n++) {
			    int ngb = Ngb.neighborIndex(n, xyz, nox, noy, noz);
				float xn = templateMappingImage[ngb+X*noxyz];
				float yn = templateMappingImage[ngb+Y*noxyz];
				float zn = templateMappingImage[ngb+Z*noxyz];
			    float levelngb = ImageInterpolation.linearInterpolation(templateLevelsetImage, 1e13f, xn,yn,zn, ntx,nty,ntz);
			    if ( (levelxyz>0) != (levelngb>0) ) {
			        // boundary between xyz and ngb
			        if (Numerics.abs(levelxyz)<Numerics.abs(levelngb)) {
			            // only label if the smallest of the two
			            boundary[xyz] = true;
			        }
			    }
			}
		}
		
		byte LINEAR = 1;
		byte NEAREST = 2;
		byte interp = NEAREST;
		if (interpParam.equals("linear")) interp = LINEAR;

		byte RAW = 1;
		byte PROJECTED = 2;
		byte sampling = PROJECTED;
		if (methodParam.equals("raw")) sampling = RAW;
		
		float maskval = 1e13f;
		
		// sample from output space to template to source
		mappedDataImage = new float[noxyz*nst];
		mappedMaskImage = new int[noxyz];
		float[] proj = new float[3];
		for (int xyz=0;xyz<noxyz;xyz++) if (boundary[xyz]) {
			// get template coordinates
			float xt = templateMappingImage[xyz+X*noxyz];
			float yt = templateMappingImage[xyz+Y*noxyz];
			float zt = templateMappingImage[xyz+Z*noxyz];
			
			// find the corresponding point in original subject space
			float xs = ImageInterpolation.linearInterpolation(sourceMappingImage, X*ntxyz, maskval, xt,yt,zt, ntx,nty,ntz);
			float ys = ImageInterpolation.linearInterpolation(sourceMappingImage, Y*ntxyz, maskval, xt,yt,zt, ntx,nty,ntz);
			float zs = ImageInterpolation.linearInterpolation(sourceMappingImage, Z*ntxyz, maskval, xt,yt,zt, ntx,nty,ntz);
			
			// project to the closest surface point or sample raw data
			if (sampling==PROJECTED) {
				projectToLevelset(sourceLevelsetImage, new float[]{xs,ys,zs}, proj);
				xs = proj[X];	
				ys = proj[Y];	
				zs = proj[Z];	
			}
			if (interp==NEAREST) {
				for (int t=0;t<nst;t++) {
					mappedDataImage[xyz+t*noxyz] = ImageInterpolation.nearestNeighborInterpolation(sourceContrastImage, t*nsxyz, sourceMask, maskval, xs,ys,zs, nsx,nsy,nsz);
				}
			} else if (interp==LINEAR) {
				for (int t=0;t<nst;t++) {
					mappedDataImage[xyz+t*noxyz] = ImageInterpolation.linearInterpolation(sourceContrastImage, t*nsxyz, sourceMask, maskval, xs,ys,zs, nsx,nsy,nsz);
				}
			}
			// update the mask
			if (mappedDataImage[xyz]!=maskval) mappedMaskImage[xyz] = 1;
			else {
				mappedMaskImage[xyz] = 0;
				for (int t=0;t<nst;t++) mappedDataImage[xyz+t*noxyz] = 0.0f;
			}
		}
		return;
	}

	private final float projectToLevelset(float[] levelset, float[] pt0, float[] pt) {
		
		// best gradient approximation among various choices (most regular)
		double I = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z], nsx, nsy, nsz);
		double Imx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X]-1.0f, pt0[Y], pt0[Z], nsx, nsy, nsz);
		double Ipx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X]+1.0f, pt0[Y], pt0[Z], nsx, nsy, nsz);
		double Imy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y]-1.0f, pt0[Z], nsx, nsy, nsz);
		double Ipy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y]+1.0f, pt0[Z], nsx, nsy, nsz);
		double Imz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z]-1.0f, nsx, nsy, nsz);
		double Ipz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z]+1.0f, nsx, nsy, nsz);
		
		//if (Double.isNaN(I) || Double.isNaN(Imx) || Double.isNaN(Ipx) || Double.isNaN(Imy) || Double.isNaN(Ipy) || Double.isNaN(Imz) || Double.isNaN(Ipz))
		//	System.out.println("NaN: levelset ("+pt0[X]+", "+pt0[Y]+", "+pt0[Z]+")");
		
		double length = I;

		double Dx = 0.5*(Ipx-Imx);
		double Dy = 0.5*(Ipy-Imy);
		double Dz = 0.5*(Ipz-Imz);
		
		double grad = FastMath.sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
		
		//if (Double.isNaN(grad)) System.out.println("NaN: gradient ("+Dx+", "+Dy+", "+Dz+")");

		double res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z], nsx, nsy, nsz);
		
		int t=0;

		if (grad > 0.01) {
			// closed form approximation to the closest point on the N-th layer surface 
			// (accurate if close enough to approximate the surface by a plane)
			pt[X] = (float)(pt0[X] - length*Dx/grad);
			pt[Y] = (float)(pt0[Y] - length*Dy/grad);
			pt[Z] = (float)(pt0[Z] - length*Dz/grad);
			
			//if (Float.isNaN(pt[X]) || Float.isNaN(pt[Y]) || Float.isNaN(pt[Z]))
			//	System.out.println("NaN: coord ("+pt0[X]+","+pt0[Y]+", "+pt0[Z]+"| "+length+"x "+Dx+", "+Dy+", "+Dz+" / "+grad);
			
			// correction: go back halfway and recompute gradient, etc: should converge toward the surface
			// approaches using small steps are not very good (unstable)
			// approaches doing the full correction from the new solution are not very stable either
			res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nsx, nsy, nsz);
			while (res*res>0.0001 && t<100 && grad>0.01) {
				t++;
				//System.out.print("-");
				
				// correction: come back halfway, recompute (going further back does not improve the result)
				pt[X] = (float)(pt[X] + 0.5*length*Dx/grad);
				pt[Y] = (float)(pt[Y] + 0.5*length*Dy/grad);
				pt[Z] = (float)(pt[Z] + 0.5*length*Dz/grad);
				
				//if (Float.isNaN(pt[X]) || Float.isNaN(pt[Y]) || Float.isNaN(pt[Z]))
				//	System.out.println("NaN: coord 0.5 ("+pt[X]+","+pt[Y]+", "+pt[Z]+"| "+length+"x "+Dx+", "+Dy+", "+Dz+" / "+grad);
			
				// correction: project once more onto the curve (with / without minimization of distance to origin point)
				I = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nsx, nsy, nsz);
				Imx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]-1.0f, pt[Y], pt[Z], nsx, nsy, nsz);
				Ipx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]+1.0f, pt[Y], pt[Z], nsx, nsy, nsz);
				Imy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]-1.0f, pt[Z], nsx, nsy, nsz);
				Ipy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]+1.0f, pt[Z], nsx, nsy, nsz);
				Imz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]-1.0f, nsx, nsy, nsz);
				Ipz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]+1.0f, nsx, nsy, nsz);
				
				//if (Double.isNaN(I) || Double.isNaN(Imx) || Double.isNaN(Ipx) || Double.isNaN(Imy) || Double.isNaN(Ipy) || Double.isNaN(Imz) || Double.isNaN(Ipz))
				//	System.out.println("NaN: levelset ("+pt[X]+", "+pt[Y]+", "+pt[Z]+")");

				length = I;
				Dx = 0.5*(Ipx-Imx);
				Dy = 0.5*(Ipy-Imy);
				Dz = 0.5*(Ipz-Imz);
				
				grad = FastMath.sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
				//if (Double.isNaN(grad)) System.out.println("NaN: gradient ("+Dx+", "+Dy+", "+Dz+")");

				if (grad > 0.01) {
					pt[X] = (float)(pt[X] - length*Dx/grad);
					pt[Y] = (float)(pt[Y] - length*Dy/grad);
					pt[Z] = (float)(pt[Z] - length*Dz/grad);
				
					//if (Float.isNaN(pt[X]) || Float.isNaN(pt[Y]) || Float.isNaN(pt[Z]))
					//	System.out.println("NaN: coord 2 ("+pt[X]+","+pt[Y]+", "+pt[Z]+"| "+length+"x "+Dx+", "+Dy+", "+Dz+" / "+grad);
				}
				res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nsx, nsy, nsz);
			}

		} else {
			//System.out.print("o");
			pt[X] = pt0[X];
			pt[Y] = pt0[Y];
			pt[Z] = pt0[Z];
		}
		//System.out.println("x");
				
		return (float)res;
	}

}
