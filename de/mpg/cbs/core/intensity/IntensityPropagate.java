package de.mpg.cbs.core.intensity;

import java.net.URL;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis Bazin
 */
 
public class IntensityPropagate {

	// parameters
	private		static final String[]	normtypes = {"max","mean","min"};
	
	// variables
	private float[] inputImage = null;
	private byte[] maskImage = null;
	private float[] resultImage;
	private	String	 normParam = "max";
	private float distParam = 5;
	
	private int nx, ny, nz, nc, nxyz;
	private float rx, ry, rz;

	// set inputs
	public final void setInputImage(float[] val) { inputImage = val; }
	public final void setMaskImage(byte[] val) { maskImage = val; }
	public final void setPropogationDistance(float val) { distParam = val; }
	
	// set generic inputs	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nc=1; nxyz=nx*ny*nz; }
	public final void setDimensions(int x, int y, int z, int c) { nx=x; ny=y; nz=z; nc=c; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; if (dim.length>3) nc=dim[3]; else nc=1; nxyz=nx*ny*nz; }
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	//set JIST definitions
	//to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity"; }
	public final String getLabel() { return "Intensity Propogation"; }
	public final String getName() { return "IntensityPropogation"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Propagates the values inside the mask (or non-zero) into the neighboring voxels"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
	
	//set outputs
	public float[] getResultImage() { return resultImage;}
	
	public void execute() {
	
		float[][][] image3 = null;
		float[][][][] image4 = null;
	
		if (nc==1) {
			image3 = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				image3[x][y][z] = inputImage[x+nx*y+nx*ny*z];
			}
		} else {
			image4 = new float[nx][ny][nz][nc];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int c=0;c<nc;c++) {
				image4[x][y][z][c] = inputImage[x+nx*y+nx*ny*z+nx*ny*nz*c];
			}
		}
		// use a mask by default
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			mask[x][y][z] = true;
		}
		// input mask or mask zero values
		if (maskImage!=null) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				mask[x][y][z] = (maskImage[x+nx*y+nx*ny*z]!=0);
			}
		} else {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (nc==1 && image3[x][y][z]==0) mask[x][y][z] = false;
				if (nc>1) {
					mask[x][y][z] = false;
					for (int c=0;c<nc;c++) if (image4[x][y][z][c]!=0) mask[x][y][z] = true;
				}
			}
		}
	
		// main algorithm
		int nd = Numerics.ceil(distParam/Numerics.min(rx,ry,rz));
	
		byte MIN = 1, MAX = 2, MEAN = 3;
		byte merge = MAX;
		if (normParam.equals("mean")) merge = MEAN;
		if (normParam.equals("min")) merge = MIN;
	
		float[][][] result3 = null;
		float[][][][] result4 = null;
		if (nc==1) result3 = new float[nx][ny][nz];
		else result4 = new float[nx][ny][nz][nc];
		byte[][][] count = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (nc==1) result3[x][y][z] = image3[x][y][z];
			else for (int c=0;c<nc;c++) result4[x][y][z][c] = image4[x][y][z][c];
		}
		for (int n=0;n<nd;n++) {
			System.out.println("step "+(n+1));
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				count[x][y][z] = 0;
			}
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
		
				for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) if (!mask[x+dx][y+dy][z+dz]) {
					//System.out.print(".");
					if (nc==1) {
						if (merge==MIN) {
							if (count[x+dx][y+dy][z+dz]==0) result3[x+dx][y+dy][z+dz] = result3[x][y][z];
							else result3[x+dx][y+dy][z+dz] = Numerics.min(result3[x+dx][y+dy][z+dz], result3[x][y][z]);
						} else if (merge==MAX) {
							if (count[x+dx][y+dy][z+dz]==0) result3[x+dx][y+dy][z+dz] = result3[x][y][z];
							else result3[x+dx][y+dy][z+dz] = Numerics.max(result3[x+dx][y+dy][z+dz], result3[x][y][z]);
						} else if (merge==MEAN) {
							result3[x+dx][y+dy][z+dz] += result3[x][y][z];
						}
					} else {
						for (int c=0;c<nc;c++) {
							if (merge==MIN) {
								if (count[x+dx][y+dy][z+dz]==0) result4[x+dx][y+dy][z+dz][c] = result4[x][y][z][c];
								else result4[x+dx][y+dy][z+dz][c] = Numerics.min(result4[x+dx][y+dy][z+dz][c], result4[x][y][z][c]);
							} else if (merge==MAX) {
								if (count[x+dx][y+dy][z+dz]==0) result4[x+dx][y+dy][z+dz][c] = result4[x][y][z][c];
								else result4[x+dx][y+dy][z+dz][c] = Numerics.max(result4[x+dx][y+dy][z+dz][c], result4[x][y][z][c]);
							} else if (merge==MEAN) {
								result4[x+dx][y+dy][z+dz][c] += result4[x][y][z][c];
							}
						}
					}					
					count[x+dx][y+dy][z+dz]++;
				}
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (count[x][y][z]>0) {
				if (merge==MEAN) {
					if (nc==1) result3[x][y][z] /= (float)count[x][y][z];
					else for (int c=0;c<nc;c++) result4[x][y][z][c] /= (float)count[x][y][z];
				}
				mask[x][y][z] = true;
			}
		}

		resultImage = new float[nx*ny*nz*nc];
		if (nc==1) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				resultImage[x+nx*y+nx*ny*z] = result3[x][y][z];
			}
		} else {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int c=0;c<nc;c++) {
				resultImage[x+nx*y+nx*ny*z+nx*ny*nz*c] = result4[x][y][z][c];
			}
		}
	}
}
