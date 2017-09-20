package de.mpg.cbs.core.registration;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.special.Erf;


/*
 * @author Pierre-Louis Bazin
 */
public class RegistrationTargetBasedReorientation {

	// jist containers
	private float[]	headImage;
	private float[]	extraImage;
	private byte[] targetImage;
	
	private float distanceParam = 35.0f;
	private float sizeParam = 30.0f;
	
	public final byte mX = 1; 
	public final byte pX = 2; 
	public final byte mY = 3; 
	public final byte pY = 4; 
	public final byte mZ = 5; 
	public final byte pZ = 6; 
	private byte neckdirParam = 0;
	
	private float[] reorientImage = null;
	private byte[] locatorImage = null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;
	
	private int nlabels = 0;
	
	// create inputs
	public final void setHeadProbabilityImage(float[] val) { headImage = val; }
	public final void setExtraImage(float[] val) { extraImage = val; }
	public final void setTargetImage(byte[] val) { targetImage = val; }
	public final void setDistance_mm(float val) { distanceParam = val; }
	public final void setNeckDirection(byte val) { neckdirParam = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Registration.devel"; }
	public final String getLabel() { return "Target-based Reorientation"; }
	public final String getName() { return "TargetBasedReorientation"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Netherlands Institute for Neuroscience"; }
	public final String getDescription() { return "Finds the closest location on the head to reach a target and reorient the image"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.1"; };
	
	public final float[] getReorientedImage() { return reorientImage; }
	public final byte[] getLocatorImage() { return locatorImage; }
	
	public	static	final	byte	X = 0;
	public	static	final	byte	Y = 1;
	public	static	final	byte	Z = 2;
	
	
	public void execute() {
		
		// 1. Separate bg / fg in original image
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (headImage[xyz]<0.5) bgmask[xyz] = true;
		}

		// 2. Get level set surfaces
		
		// main background piece
		bgmask = ObjectLabeling.largestObject(bgmask, nx, ny, nz, 6);
		bgmask = ObjectMorphology.erodeObject(bgmask, nx,ny,nz, 3,3,3);
		bgmask = ObjectLabeling.largestObject(bgmask, nx, ny, nz, 6);
		bgmask = ObjectMorphology.dilateObject(bgmask, nx,ny,nz, 3,3,3);
		
		// erase small protrusions
		bgmask = ObjectMorphology.dilateObject(bgmask, nx,ny,nz, 3,3,3);
		bgmask = ObjectMorphology.erodeObject(bgmask, nx,ny,nz, 3,3,3);
		for (int xyz=0;xyz<nxyz;xyz++) bgmask[xyz] = (!bgmask[xyz]);
		
		// gdm smoothing?
		float[] gdm = ObjectTransforms.fastMarchingDistanceFunction(bgmask, nx, ny, nz);
	
		boolean[] mask = ObjectMorphology.dilateObject(bgmask, nx,ny,nz, 3,3,3);
		SmoothGdm evolve = new SmoothGdm(gdm, gdm, nx, ny, nz, rx, ry, rz, mask, 0.1f, 0.4f, "no",null);
		evolve.fastMarchingReinitialization(true);
		BasicInfo.displayMessage("level set segmentation...\n");
		evolve.evolveNarrowBand(500, 0.001f);
		
		gdm = evolve.getLevelSet();	
		for (int xyz=0;xyz<nxyz;xyz++) bgmask[xyz] = (gdm[xyz]<0);
		
		// extend to avoid target placement in the neck
		if (neckdirParam>0) {
			bgmask = ObjectMorphology.erodeObject(bgmask, nx,ny,nz, 3,3,3);
			if (neckdirParam==mX) {
				for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
					boolean max = false;
					for (int x=nx-1;x>=0;x--) {
						int xyz = x+nx*y+nx*ny*z;
						if (bgmask[xyz]) max = bgmask[xyz];
						bgmask[xyz] = max;
					}
				}
			} else
			if (neckdirParam==pX) {
				for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					boolean max = false;
					for (int x=0;x<nx;x++) {
						int xyz = x+nx*y+nx*ny*z;
						if (bgmask[xyz]) max = bgmask[xyz];
						bgmask[xyz] = max;
					}
				}
			}
			if (neckdirParam==mY) {
				for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
					boolean max = false;
					for (int y=ny-1;y>=0;y--) {
						int xyz = x+nx*y+nx*ny*z;
						if (bgmask[xyz]) max = bgmask[xyz];
						bgmask[xyz] = max;
					}
				}
			} else
			if (neckdirParam==pY) {
				for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
					boolean max = false;
					for (int y=0;y<ny;y++) {
						int xyz = x+nx*y+nx*ny*z;
						if (bgmask[xyz]) max = bgmask[xyz];
						bgmask[xyz] = max;
					}
				}
			}
			if (neckdirParam==mZ) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
					boolean max = false;
					for (int z=nz-1;z>=0;z--) {
						int xyz = x+nx*y+nx*ny*z;
						if (bgmask[xyz]) max = bgmask[xyz];
						bgmask[xyz] = max;
					}
				}
			} else
			if (neckdirParam==pZ) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
					boolean max = false;
					for (int z=0;z<nz;z++) {
						int xyz = x+nx*y+nx*ny*z;
						if (bgmask[xyz]) max = bgmask[xyz];
						bgmask[xyz] = max;
					}
				}
			}
			bgmask = ObjectMorphology.dilateObject(bgmask, nx,ny,nz, 3,3,3);
		}
		float[] bgdist = ObjectTransforms.fastMarchingDistanceFunction(bgmask, nx, ny, nz);
		
		boolean[] tgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) if (targetImage[xyz]>0) tgmask[xyz] = true;
		float[] tgdist = ObjectTransforms.fastMarchingDistanceFunction(tgmask, nx, ny, nz);
		
		// 3. Find intersection(s): get the minimum target distance on boundary
		boolean[] roi = new boolean[nxyz];
		float dist = distanceParam/rx;
		float mindist = nx+ny+nz;
		int id = 0;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (Numerics.abs(bgdist[xyz])<1.0f) if (Numerics.abs(tgdist[xyz])<mindist) {
			//if (Numerics.abs(bgdist[xyz])<1.0f) if (Numerics.abs(tgdist[xyz]-dist)<mindist) {
				// find closest point to the ideal distance? or the closest in general?
				mindist = Numerics.abs(tgdist[xyz]);
				id = xyz;
			}
		}
		roi[id] = true;
		BasicInfo.displayMessage("intersection center ("+id+")\n");
		
		// find centroid
		float xc = 0.0f, yc = 0.0f, zc = 0.0f, nc = 0.0f;
		float[] p0 = new float[3];
		float[] pC = new float[3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (roi[xyz]) { xc += x; yc += y; zc += z; nc++; }
		}
		if (nc>0) {
			xc /= nc;
			yc /= nc;
			zc /= nc;
		
			// project onto head level set
			pC[X] = xc; pC[Y] = yc; pC[Z] = zc;
			projectToLevelset(bgdist, pC, p0, false);
		} else {
			p0[X] = nx/2;
			p0[Y] = ny/2;
			p0[Z] = nz/2;
		}
		// find target centroid
		float xt = 0.0f, yt = 0.0f, zt = 0.0f, nt = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (tgmask[xyz]) { xt += x; yt += y; zt += z; nt++; }
		}
		xt /= nt;
		yt /= nt;
		zt /= nt;
		
		BasicInfo.displayMessage("intersection center ("+xc+", "+yc+", "+zc+")\n");
		BasicInfo.displayMessage("projection center ("+p0[X]+", "+p0[Y]+", "+p0[Z]+")\n");
		BasicInfo.displayMessage("target center ("+xt+", "+yt+", "+zt+")\n");
				
		xc = p0[X];
		yc = p0[Y];
		zc = p0[Z];
		
		// build a circular target patch
		boolean[] pcmask = new boolean[nxyz];
		float pcdist2 = Numerics.square(0.5f*sizeParam/rx);
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (Numerics.abs(bgdist[xyz])<1.8f) if ( (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc) < pcdist2) {
				pcmask[xyz] = true;
			}
		}
		
		// move target to match distance
		float dist0 = (float)FastMath.sqrt( (xc-xt)*(xc-xt) + (yc-yt)*(yc-yt) + (zc-zt)*(zc-zt) );	
		BasicInfo.displayMessage("distance: "+dist0+"\n");
		xt = xc + dist/dist0*(xt-xc);
		yt = yc + dist/dist0*(yt-yc);
		zt = zc + dist/dist0*(zt-zc);
		
		BasicInfo.displayMessage("moved target center ("+xt+", "+yt+", "+zt+")\n");
		
		/*
		// 4. search for optimal placement (location + orientation)?
		int delta = 5;
		for (int dx=-delta;dx<=delta;dx++) for (int dy=-delta;dy<=delta;dy++) for (int dz=-delta;dz<delta;dz++) {
			// find local tangent plane, normal (given some size...)
			
		}
		*/
		
		// 5. Reorient data
		double[] rotx = new double[3];
		double[] roty = new double[3];
		double[] rotz = new double[3];
		
		if (nc>0) {
			rotx[X] = xc - xt;
			rotx[Y] = yc - yt;
			rotx[Z] = zc - zt;
		} else {
			rotx[X] = 0.0f;
			rotx[Y] = 0.0f;
			rotx[Z] = 1.0f;
		}
		double normx = FastMath.sqrt(rotx[X]*rotx[X] + rotx[Y]*rotx[Y] + rotx[Z]*rotx[Z]);
		rotx[X] /= normx;
		rotx[Y] /= normx;
		rotx[Z] /= normx;
		
		if (rotx[Y]*rotx[Y]<=rotx[Z]*rotx[Z]) {
			roty[X] = 0.0-rotx[Y]*rotx[X];
			roty[Y] = 1.0-rotx[Y]*rotx[Y];
			roty[Z] = 0.0-rotx[Y]*rotx[Z];
		} else {
			roty[X] = 0.0-rotx[Z]*rotx[X];
			roty[Y] = 0.0-rotx[Z]*rotx[Y];
			roty[Z] = 1.0-rotx[Z]*rotx[Z];
		}
		double normy = FastMath.sqrt(roty[X]*roty[X] + roty[Y]*roty[Y] + roty[Z]*roty[Z]);
		roty[X] /= normy;
		roty[Y] /= normy;
		roty[Z] /= normy;
		
		// vector product to get the rotation matrix
		rotz[X] = rotx[Y]*roty[Z]-rotx[Z]*roty[Y];
		rotz[Y] = rotx[Z]*roty[X]-rotx[X]*roty[Z];
		rotz[Z] = rotx[X]*roty[Y]-rotx[Y]*roty[X];
		double normz = FastMath.sqrt(rotz[X]*rotz[X] + rotz[Y]*rotz[Y] + rotz[Z]*rotz[Z]);
		rotz[X] /= normz;
		rotz[Y] /= normz;
		rotz[Z] /= normz;
		
		BasicInfo.displayMessage("rotation matrix \n ("+rotx[X]+", "+rotx[Y]+", "+rotx[Z]+")\n");
		BasicInfo.displayMessage(" ("+roty[X]+", "+roty[Y]+", "+roty[Z]+")\n");
		BasicInfo.displayMessage(" ("+rotz[X]+", "+rotz[Y]+", "+rotz[Z]+")\n");
		
		
		// rotate image from the center
		reorientImage = new float[nxyz];
		locatorImage = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			// transform?
			/*
			float xr = (float)(rotx[X]*(x-nx/2) + rotx[Y]*(y-ny/2) + rotx[Z]*(z-nz/2) + nx/2);
			float yr = (float)(roty[X]*(x-nx/2) + roty[Y]*(y-ny/2) + roty[Z]*(z-nz/2) + ny/2);
			float zr = (float)(rotz[X]*(x-nx/2) + rotz[Y]*(y-ny/2) + rotz[Z]*(z-nz/2) + nz/2);
			*/
			float xr = (float)(rotx[X]*(x-nx/2) + roty[X]*(y-ny/2) + rotz[X]*(z-nz/2) + nx/2);
			float yr = (float)(rotx[Y]*(x-nx/2) + roty[Y]*(y-ny/2) + rotz[Y]*(z-nz/2) + ny/2);
			float zr = (float)(rotx[Z]*(x-nx/2) + roty[Z]*(y-ny/2) + rotz[Z]*(z-nz/2) + nz/2);
		
			int xyz = x+nx*y+nx*ny*z;
			reorientImage[xyz] = ImageInterpolation.linearInterpolation(extraImage, 0.0f, xr, yr, zr, nx, ny, nz);
			//reorientImage[xyz] = bgdist[xyz];
			/* reoriented or initial */
			// get the head outline
			boolean bgval = ImageInterpolation.nearestNeighborInterpolation(bgmask, false, xr, yr, zr, nx, ny, nz);
			//float phead = ImageInterpolation.linearInterpolation(bgdist, 1.0f, xr, yr, zr, nx, ny, nz);
			if (bgval==true) locatorImage[xyz] = 1;
			// paste the target
			byte tgval = ImageInterpolation.nearestNeighborInterpolation(targetImage, (byte)0, xr, yr, zr, nx, ny, nz);
			if (tgval>0) locatorImage[xyz] = 2;
			// add the intersection roi
			boolean pcval = ImageInterpolation.nearestNeighborInterpolation(pcmask, false, xr, yr, zr, nx, ny, nz);
			if (pcval) locatorImage[xyz] = 3;
			/*
			if (bgmask[xyz]) locatorImage[xyz] = 1;
			if (targetImage[xyz]>0) locatorImage[xyz] = 2;
			if (roi[xyz]) locatorImage[xyz] = 3;
			*/
		}
			
		// indicate the source / target point (inverse transform)
		/*
		float xrc = (float)(rotx[X]*(xc-nx/2) + roty[X]*(yc-ny/2) + rotz[X]*(zc-nz/2) + nx/2);
		float yrc = (float)(rotx[Y]*(xc-nx/2) + roty[Y]*(yc-ny/2) + rotz[Y]*(zc-nz/2) + ny/2);
		float zrc = (float)(rotx[Z]*(xc-nx/2) + roty[Z]*(yc-ny/2) + rotz[Z]*(zc-nz/2) + nz/2);
		
		float xrt = (float)(rotx[X]*(xt-nx/2) + roty[X]*(yt-ny/2) + rotz[X]*(zt-nz/2) + nx/2);
		float yrt = (float)(rotx[Y]*(xt-nx/2) + roty[Y]*(yt-ny/2) + rotz[Y]*(zt-nz/2) + ny/2);
		float zrt = (float)(rotx[Z]*(xt-nx/2) + roty[Z]*(yt-ny/2) + rotz[Z]*(zt-nz/2) + nz/2);
		*/
		float xrc = (float)(rotx[X]*(xc-nx/2) + rotx[Y]*(yc-ny/2) + rotx[Z]*(zc-nz/2) + nx/2);
		float yrc = (float)(roty[X]*(xc-nx/2) + roty[Y]*(yc-ny/2) + roty[Z]*(zc-nz/2) + ny/2);
		float zrc = (float)(rotz[X]*(xc-nx/2) + rotz[Y]*(yc-ny/2) + rotz[Z]*(zc-nz/2) + nz/2);
			
		float xrt = (float)(rotx[X]*(xt-nx/2) + rotx[Y]*(yt-ny/2) + rotx[Z]*(zt-nz/2) + nx/2);
		float yrt = (float)(roty[X]*(xt-nx/2) + roty[Y]*(yt-ny/2) + roty[Z]*(zt-nz/2) + ny/2);
		float zrt = (float)(rotz[X]*(xt-nx/2) + rotz[Y]*(yt-ny/2) + rotz[Z]*(zt-nz/2) + nz/2);
			
		
		int xyzc, xyzt;
		xyzc = Numerics.round(xrc) + nx*Numerics.round(yrc) + nx*ny*Numerics.round(zrc);
		xyzt = Numerics.round(xrt) + nx*Numerics.round(yrt) + nx*ny*Numerics.round(zrt);
		
		// raw oriwentation
		//xyzc = Numerics.round(xc) + nx*Numerics.round(yc) + nx*ny*Numerics.round(zc);
		//xyzt = Numerics.round(xt) + nx*Numerics.round(yt) + nx*ny*Numerics.round(zt);
		
		
		locatorImage[xyzc] = 5;	
		
		locatorImage[xyzc-1] = 4;			
		locatorImage[xyzc+1] = 4;			
		locatorImage[xyzc-nx] = 4;			
		locatorImage[xyzc+nx] = 4;			
		locatorImage[xyzc-nx*ny] = 4;			
		locatorImage[xyzc+nx*ny] = 4;			
		
		locatorImage[xyzt] = 7;	
		
		locatorImage[xyzt-1] = 6;			
		locatorImage[xyzt+1] = 6;			
		locatorImage[xyzt-nx] = 6;			
		locatorImage[xyzt+nx] = 6;			
		locatorImage[xyzt-nx*ny] = 6;			
		locatorImage[xyzt+nx*ny] = 6;			
			
		// 6. Build region of estimated stimulation
		
		return;
	}
	
	private final float projectToLevelset(float[] levelset, float[] pt0, float[] pt, boolean iterate) {
		
		// best gradient approximation among various choices (most regular)
		double I = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z], nx, ny, nz);
		double Imx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X]-1.0f, pt0[Y], pt0[Z], nx, ny, nz);
		double Ipx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X]+1.0f, pt0[Y], pt0[Z], nx, ny, nz);
		double Imy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y]-1.0f, pt0[Z], nx, ny, nz);
		double Ipy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y]+1.0f, pt0[Z], nx, ny, nz);
		double Imz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z]-1.0f, nx, ny, nz);
		double Ipz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z]+1.0f, nx, ny, nz);
		
		//if (Double.isNaN(I) || Double.isNaN(Imx) || Double.isNaN(Ipx) || Double.isNaN(Imy) || Double.isNaN(Ipy) || Double.isNaN(Imz) || Double.isNaN(Ipz))
		//	System.out.println("NaN: levelset ("+pt0[X]+", "+pt0[Y]+", "+pt0[Z]+")");
		
		double length = I;

		double Dx = 0.5*(Ipx-Imx);
		double Dy = 0.5*(Ipy-Imy);
		double Dz = 0.5*(Ipz-Imz);
		
		double grad = FastMath.sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
		
		//if (Double.isNaN(grad)) System.out.println("NaN: gradient ("+Dx+", "+Dy+", "+Dz+")");

		double res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z], nx, ny, nz);
		
		int t=0;

		if (grad > 0.01) {
			// closed form approximation to the closest point on the N-th layer surface 
			// (accurate if close enough to approximate the surface by a plane)
			pt[X] = (float)(pt0[X] - length*Dx/grad);
			pt[Y] = (float)(pt0[Y] - length*Dy/grad);
			pt[Z] = (float)(pt0[Z] - length*Dz/grad);
			
			//if (Float.isNaN(pt[X]) || Float.isNaN(pt[Y]) || Float.isNaN(pt[Z]))
			//	System.out.println("NaN: coord ("+pt0[X]+","+pt0[Y]+", "+pt0[Z]+"| "+length+"x "+Dx+", "+Dy+", "+Dz+" / "+grad);
			
			if (iterate) {
				// correction: go back halfway and recompute gradient, etc: should converge toward the surface
				// approaches using small steps are not very good (unstable)
				// approaches doing the full correction from the new solution are not very stable either
				res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				while (res*res>0.0001 && t<100 && grad>0.01) {
					t++;
					//System.out.print(".");
					
					// correction: come back halfway, recompute (going further back does not improve the result)
					pt[X] = (float)(pt[X] + 0.5*length*Dx/grad);
					pt[Y] = (float)(pt[Y] + 0.5*length*Dy/grad);
					pt[Z] = (float)(pt[Z] + 0.5*length*Dz/grad);
					
					//if (Float.isNaN(pt[X]) || Float.isNaN(pt[Y]) || Float.isNaN(pt[Z]))
					//	System.out.println("NaN: coord 0.5 ("+pt[X]+","+pt[Y]+", "+pt[Z]+"| "+length+"x "+Dx+", "+Dy+", "+Dz+" / "+grad);
				
					// correction: project once more onto the curve (with / without minimization of distance to origin point)
					I = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
					Imx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]-1.0f, pt[Y], pt[Z], nx, ny, nz);
					Ipx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]+1.0f, pt[Y], pt[Z], nx, ny, nz);
					Imy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]-1.0f, pt[Z], nx, ny, nz);
					Ipy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]+1.0f, pt[Z], nx, ny, nz);
					Imz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]-1.0f, nx, ny, nz);
					Ipz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]+1.0f, nx, ny, nz);
					
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
					res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				}
			}
		} else {
			//System.err.print(".");
			pt[X] = pt0[X];
			pt[Y] = pt0[Y];
			pt[Z] = pt0[Z];
		}
		return (float)res;
	}

}