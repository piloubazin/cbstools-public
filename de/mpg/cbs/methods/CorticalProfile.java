package de.mpg.cbs.methods;

import java.util.*;
import java.io.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;

import org.apache.commons.math3.util.FastMath;


/**
 * Handles various cortical profile computations
 * 
 * @author Pierre-Louis Bazin
 * 
 */
public class CorticalProfile {

	private float[][] profile;
	
	private int nlayers;
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;
	private static final byte D = 3;

	// debug variables
	private static final boolean debug=false;
	private int mtrial, Mtrial;
	
	
	public CorticalProfile(int nlayers_, int nx_, int ny_, int nz_, float rx_, float ry_, float rz_) {
		nlayers = nlayers_;
		profile = new float[nlayers+1][D];
		nx = nx_;
		ny = ny_;
		nz = nz_;
		nxyz = nx*ny*nz;
		rx = rx_;
		ry = ry_;
		rz = rz_;
	}

	public final float[] getPt(int l) { return profile[l]; }
	
	public final float[][] getProfile() { return profile; }
	
	public final int getMeanTrial() { return mtrial; }
	
	public final int getMaxTrial() { return Mtrial; }
	
	
	private final float projectToLevelset(float[] levelset, float[] pt0, float[] pt) {
		
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

		} else {
			//System.err.print(".");
			pt[X] = pt0[X];
			pt[Y] = pt0[Y];
			pt[Z] = pt0[Z];
		}
		mtrial += t;
		if (t>Mtrial) Mtrial = t;
		
		return (float)res;
	}
	
	private final float projectToLevelset(float[] levelset, int offset, float[] pt0, float[] pt) {
		
		// best gradient approximation among various choices (most regular)
		double I = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt0[X], pt0[Y], pt0[Z], nx, ny, nz);
		double Imx = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt0[X]-1.0f, pt0[Y], pt0[Z], nx, ny, nz);
		double Ipx = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt0[X]+1.0f, pt0[Y], pt0[Z], nx, ny, nz);
		double Imy = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt0[X], pt0[Y]-1.0f, pt0[Z], nx, ny, nz);
		double Ipy = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt0[X], pt0[Y]+1.0f, pt0[Z], nx, ny, nz);
		double Imz = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt0[X], pt0[Y], pt0[Z]-1.0f, nx, ny, nz);
		double Ipz = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt0[X], pt0[Y], pt0[Z]+1.0f, nx, ny, nz);
		
		//if (Double.isNaN(I) || Double.isNaN(Imx) || Double.isNaN(Ipx) || Double.isNaN(Imy) || Double.isNaN(Ipy) || Double.isNaN(Imz) || Double.isNaN(Ipz))
		//	System.out.println("NaN: levelset ("+pt0[X]+", "+pt0[Y]+", "+pt0[Z]+")");
		
		double length = I;

		double Dx = 0.5*(Ipx-Imx);
		double Dy = 0.5*(Ipy-Imy);
		double Dz = 0.5*(Ipz-Imz);
		
		double grad = FastMath.sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
		
		//if (Double.isNaN(grad)) System.out.println("NaN: gradient ("+Dx+", "+Dy+", "+Dz+")");

		double res = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt0[X], pt0[Y], pt0[Z], nx, ny, nz);
		
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
			res = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
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
				I = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				Imx = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X]-1.0f, pt[Y], pt[Z], nx, ny, nz);
				Ipx = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X]+1.0f, pt[Y], pt[Z], nx, ny, nz);
				Imy = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X], pt[Y]-1.0f, pt[Z], nx, ny, nz);
				Ipy = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X], pt[Y]+1.0f, pt[Z], nx, ny, nz);
				Imz = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X], pt[Y], pt[Z]-1.0f, nx, ny, nz);
				Ipz = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X], pt[Y], pt[Z]+1.0f, nx, ny, nz);
				
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
				res = ImageInterpolation.linearInterpolation(levelset, offset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
			}

		} else {
			//System.err.print(".");
			pt[X] = pt0[X];
			pt[Y] = pt0[Y];
			pt[Z] = pt0[Z];
		}
		mtrial += t;
		if (t>Mtrial) Mtrial = t;
		
		return (float)res;
	}
	
	public final float computeTrajectory(float[][] layers, int x, int y, int z) {
		int xyz = x + nx*y + nx*ny*z;
		
		mtrial = 0;
		Mtrial = 0;
		
		// first, find the closest surface to the point
		int l1=0;
		for (int l=1;l<=nlayers;l++) {
			if (Numerics.abs(layers[l][xyz]) < Numerics.abs(layers[l1][xyz]) ) l1 = l;
		}
		// first propagation l1 -> 0
		float res = projectToLevelset(layers[l1], new float[]{x,y,z}, profile[l1]);

		for (int l=l1-1;l>=0;l--) {
			float newres = projectToLevelset(layers[l], profile[l+1], profile[l]);
			if (newres>res) res = newres;
		}
		for (int l=l1+1;l<=nlayers;l++) {
			float newres = projectToLevelset(layers[l], profile[l-1], profile[l]);
			if (newres>res) res = newres;
		}
		mtrial /= (nlayers+1);
		
		return res;
	}
	
	public final float computeTrajectory(float[] layers, int x, int y, int z) {
		int xyz = x + nx*y + nx*ny*z;
		
		mtrial = 0;
		Mtrial = 0;
		
		// first, find the closest surface to the point
		int l1=0;
		for (int l=1;l<=nlayers;l++) {
			if (Numerics.abs(layers[l*nxyz+xyz]) < Numerics.abs(layers[l1*nxyz+xyz]) ) l1 = l;
		}
		// first propagation l1 -> 0
		float res = projectToLevelset(layers, l1*nxyz, new float[]{x,y,z}, profile[l1]);

		for (int l=l1-1;l>=0;l--) {
			float newres = projectToLevelset(layers, l*nxyz, profile[l+1], profile[l]);
			if (newres>res) res = newres;
		}
		for (int l=l1+1;l<=nlayers;l++) {
			float newres = projectToLevelset(layers, l*nxyz, profile[l-1], profile[l]);
			if (newres>res) res = newres;
		}
		mtrial /= (nlayers+1);
		
		return res;
	}
	
	public final float computeTrajectory(float[][] layers, float x, float y, float z) {
		
		mtrial = 0;
		Mtrial = 0;
		
		// first, find the closest surface to the point
		int l1=0; float d1 = Numerics.abs(ImageInterpolation.linearInterpolation(layers[0], 1e9f, x,y,z, nx,ny,nz));
		for (int l=1;l<=nlayers;l++) {
			float d = Numerics.abs(ImageInterpolation.linearInterpolation(layers[l], 1e9f, x,y,z, nx,ny,nz));
			if (d < d1) { l1 = l; d1 = d; }
		}
		// first propagation l1 -> 0
		float res = projectToLevelset(layers[l1], new float[]{x,y,z}, profile[l1]);

		for (int l=l1-1;l>=0;l--) {
			float newres = projectToLevelset(layers[l], profile[l+1], profile[l]);
			if (newres>res) res = newres;
		}
		for (int l=l1+1;l<=nlayers;l++) {
			float newres = projectToLevelset(layers[l], profile[l-1], profile[l]);
			if (newres>res) res = newres;
		}
		mtrial /= (nlayers+1);
		
		return res;
	}
	
	/** total length of the profile following the individual segments */
	public final float computeLength() {
		double length = 0.0;
		for (int l=0;l<nlayers;l++) {
			length += FastMath.sqrt( (profile[l][X]-profile[l+1][X])*rx*(profile[l][X]-profile[l+1][X])*rx
								+(profile[l][Y]-profile[l+1][Y])*ry*(profile[l][Y]-profile[l+1][Y])*ry
								+(profile[l][Z]-profile[l+1][Z])*rz*(profile[l][Z]-profile[l+1][Z])*rz );
		}
		return (float)length;
	}
	
	/** average angle between profile segments (normalized in [0,1]) */
	public final float computeAvgAngle() {
		double angle = 0.0;
		for (int l=0;l<=nlayers-2;l++) {
			double prod = (profile[l+1][X]-profile[l][X])*(profile[l+2][X]-profile[l+1][X])
						 +(profile[l+1][Y]-profile[l][Y])*(profile[l+2][Y]-profile[l+1][Y])
						 +(profile[l+1][Z]-profile[l][Z])*(profile[l+2][Z]-profile[l+1][Z]);
						 
			double norm1 = FastMath.sqrt( (profile[l+1][X]-profile[l][X])*(profile[l+1][X]-profile[l][X])
									  +(profile[l+1][Y]-profile[l][Y])*(profile[l+1][Y]-profile[l][Y])
									  +(profile[l+1][Z]-profile[l][Z])*(profile[l+1][Z]-profile[l][Z]) );
			
			double norm2 = FastMath.sqrt( (profile[l+2][X]-profile[l+1][X])*(profile[l+2][X]-profile[l+1][X])
									  +(profile[l+2][Y]-profile[l+1][Y])*(profile[l+2][Y]-profile[l+1][Y])
									  +(profile[l+2][Z]-profile[l+1][Z])*(profile[l+2][Z]-profile[l+1][Z]) );
			
			angle  += FastMath.acos( Numerics.bounded(prod/(norm1*norm2), -1.0, 1.0) )/FastMath.PI;
		}
		return (float)(angle/(nlayers-2.0));
	}
	
	/** maximum angle between profile segments (normalized in [0,1]) */
	public final float computeMaxAngle() {
		double angle = 0.0;
		for (int l=0;l<=nlayers-2;l++) {
			double prod = (profile[l+1][X]-profile[l][X])*(profile[l+2][X]-profile[l+1][X])
						 +(profile[l+1][Y]-profile[l][Y])*(profile[l+2][Y]-profile[l+1][Y])
						 +(profile[l+1][Z]-profile[l][Z])*(profile[l+2][Z]-profile[l+1][Z]);
						 
			double norm1 = FastMath.sqrt( (profile[l+1][X]-profile[l][X])*(profile[l+1][X]-profile[l][X])
									  +(profile[l+1][Y]-profile[l][Y])*(profile[l+1][Y]-profile[l][Y])
									  +(profile[l+1][Z]-profile[l][Z])*(profile[l+1][Z]-profile[l][Z]) );
			
			double norm2 = FastMath.sqrt( (profile[l+2][X]-profile[l+1][X])*(profile[l+2][X]-profile[l+1][X])
									  +(profile[l+2][Y]-profile[l+1][Y])*(profile[l+2][Y]-profile[l+1][Y])
									  +(profile[l+2][Z]-profile[l+1][Z])*(profile[l+2][Z]-profile[l+1][Z]) );
			
			
			double theta = FastMath.acos( Numerics.bounded(prod/(norm1*norm2), -1.0, 1.0) )/FastMath.PI;
			if (theta>angle) angle = theta;
		}
		return (float)(angle);		
	}
	
	/** average curvature (finite difference, excludes endpoints) */
	public final float computeCurvature() {
		double curv = 0.0;
		for (int l=1;l<=nlayers-1;l++) {
			double dx = (profile[l+1][X]-profile[l-1][X])/2.0;
			double dy = (profile[l+1][Y]-profile[l-1][Y])/2.0;
			double dz = (profile[l+1][Z]-profile[l-1][Z])/2.0;
			
			double d2x = (profile[l+1][X]-2.0*profile[l][X]+profile[l-1][X]);
			double d2y = (profile[l+1][Y]-2.0*profile[l][Y]+profile[l-1][Y]);
			double d2z = (profile[l+1][Z]-2.0*profile[l][Z]+profile[l-1][Z]);
			
			double num = FastMath.sqrt( Numerics.square(d2z*dy-d2y*dz) + Numerics.square(d2x*dz-d2z*dx) + Numerics.square(d2y*dx-d2x*dy) );
			double den = FastMath.sqrt( Numerics.cube(dx*dx + dy*dy + dz*dz) );
			
			curv  += num/den;
		}
		return (float)(curv/(nlayers-2.0));
	}
	
	/** average torsion (finite difference, excludes endpoints) */
	public final float computeTorsion() {
		double tors = 0.0;
		for (int l=2;l<=nlayers-2;l++) {
			double dx = (profile[l+1][X]-profile[l-1][X])/2.0;
			double dy = (profile[l+1][Y]-profile[l-1][Y])/2.0;
			double dz = (profile[l+1][Z]-profile[l-1][Z])/2.0;
			
			double d2x = (profile[l+1][X]-2.0*profile[l][X]+profile[l-1][X]);
			double d2y = (profile[l+1][Y]-2.0*profile[l][Y]+profile[l-1][Y]);
			double d2z = (profile[l+1][Z]-2.0*profile[l][Z]+profile[l-1][Z]);
			
			double d3x = (profile[l+2][X]-2.0*profile[l+1][X]+2.0*profile[l-1][X]-profile[l-2][X])/2.0;
			double d3y = (profile[l+2][Y]-2.0*profile[l+1][Y]+2.0*profile[l-1][Y]-profile[l-2][Y])/2.0;
			double d3z = (profile[l+2][Z]-2.0*profile[l+1][Z]+2.0*profile[l-1][Z]-profile[l-2][Z])/2.0;
			
			double num = d3z*(dx*d2y-dy*d2x) + d3x*(dy*d2z-dz*d2y) + d3y*(dz*d2x-dx*d2z);
			double den = (dx*dx + dy*dy + dz*dz)*(d2x*d2x + d2y*d2y + d2z*d2z);
			
			tors  += num/den;
		}
		return (float)(tors/(nlayers-2.0));
	}
	
	/** average angle between profile segments (normalized in [0,1]) */
	public final float computeAngleWith(double vx, double vy, double vz) {
		double angle = 0.0;
		double nv = FastMath.sqrt( vx*vx + vy*vy + vz*vz );
		for (int l=1;l<=nlayers-1;l++) {
			double dx = (profile[l+1][X]-profile[l-1][X])/2.0;
			double dy = (profile[l+1][Y]-profile[l-1][Y])/2.0;
			double dz = (profile[l+1][Z]-profile[l-1][Z])/2.0;
			
			double num = vx*dx + vy*dy + vz*dz;
			double nd = FastMath.sqrt( dx*dx + dy*dy + dz*dz );
						 
			angle  += FastMath.acos( Numerics.bounded(num/(nv*nd), -1.0, 1.0) )/FastMath.PI;
		}
		return (float)(angle/(nlayers-2.0));
	}
	
	/** get normal vector at closest location to a given point */
	public final float[] computeTangentAt(double x, double y, double z) {
		double mindist = 1e9;
		int lbest = -1;
		for (int l=0;l<=nlayers;l++) {
			double dist = (profile[l][X]-x)*(profile[l][X]-x)
						 +(profile[l][Y]-y)*(profile[l][Y]-y)
						 +(profile[l][Z]-z)*(profile[l][Z]-z);
			
			 if (dist<mindist) {
				 mindist = dist;
				 lbest = l;
			 }
		}
		float[] dn = new float[3];
		if (lbest==0) {
			double dx = (profile[lbest+1][X]-profile[lbest][X]);
			double dy = (profile[lbest+1][Y]-profile[lbest][Y]);
			double dz = (profile[lbest+1][Z]-profile[lbest][Z]);
			double nd = FastMath.sqrt( dx*dx + dy*dy + dz*dz );
			
			if (nd>0.0001) {
				dn[X] = (float)(dx/nd);
				dn[Y] = (float)(dx/nd);
				dn[Z] = (float)(dx/nd);
			}
		} else
		if (lbest==nlayers) {
			double dx = (profile[lbest][X]-profile[lbest-1][X]);
			double dy = (profile[lbest][Y]-profile[lbest-1][Y]);
			double dz = (profile[lbest][Z]-profile[lbest-1][Z]);
			double nd = FastMath.sqrt( dx*dx + dy*dy + dz*dz );
			
			if (nd>0.0001) {
				dn[X] = (float)(dx/nd);
				dn[Y] = (float)(dx/nd);
				dn[Z] = (float)(dx/nd);
			}
		} else
		if (lbest>-1) {
			double dx = (profile[lbest+1][X]-profile[lbest-1][X])/2.0;
			double dy = (profile[lbest+1][Y]-profile[lbest-1][Y])/2.0;
			double dz = (profile[lbest+1][Z]-profile[lbest-1][Z])/2.0;
			double nd = FastMath.sqrt( dx*dx + dy*dy + dz*dz );
			
			if (nd>0.0001) {
				dn[X] = (float)(dx/nd);
				dn[Y] = (float)(dx/nd);
				dn[Z] = (float)(dx/nd);
			}
		}
		return dn;
	}

	/** check for numerical problems */
	public final boolean checkCoordinates() {
		boolean pb = false;
		for (int l=0;l<=nlayers;l++) {
			if (Double.isNaN(profile[l][X]) || Double.isNaN(profile[l][Y]) || Double.isNaN(profile[l][Z])) {
				System.out.println("NaN: "+l);
				pb = true;
			}
		}
		return pb;
	}
	
	
}
