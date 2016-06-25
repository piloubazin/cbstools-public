package de.mpg.cbs.core.cortex;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import java.util.LinkedList;
import java.util.List;

/*
 * @author Pierre-Louis Bazin
 */
public class CortexSurfaceMeshMapping {

	// jist containers
	private float[] 	intensityImage;
	private String 		intensityName = "contrast";
	private float[] 	inflatedSurfacePoints;
	private int[] 		inflatedSurfaceTriangles;
	private float[] 	origSurfacePoints;
	private int[] 		origSurfaceTriangles;
	private String  	mappingOption = "closest_point";
	private static final String[]		mappingTypes = {"closest_point","linear_interp","highest_value"};
	
	private float[] 	mappedInfSurfacePoints;
	private int[] 		mappedInfSurfaceTriangles;
	private float[] 	mappedInfSurfaceValues;
	
	private float[] 	mappedOrgSurfacePoints;
	private int[] 		mappedOrgSurfaceTriangles;
	private float[] 	mappedOrgSurfaceValues;
	
	private int nx, ny, nz, nt, nxyz;
	private float rx, ry, rz;

	private String surfaceConvention = "mipav";
	private static final String[] conventionTypes = {"mipav","voxel"};

	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private static final byte NONE = 10;
	private static final byte CLOSEST = 11;
	private static final byte LINEAR = 12;
	private static final byte HIGHEST = 13;
	
	// create inputs
	public final void setIntensityImage(float[] val) { intensityImage = val; }
	public final void setIntensityLabelName(String val) { intensityName = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nt=1; nxyz=nx*ny*nz; }
	public final void setDimensions(int x, int y, int z, int t) { nx=x; ny=y; nz=z; nt=t; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; if (dim.length>3) nt=dim[3]; else nt=1; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setOriginalSurfacePoints(float[] val) { origSurfacePoints = val; }
	public final void setOriginalSurfaceTriangles(int[] val) { origSurfaceTriangles = val; }
	public final void setInflatedSurfacePoints(float[] val) { inflatedSurfacePoints = val; }
	public final void setInflatedSurfaceTriangles(int[] val) { inflatedSurfaceTriangles = val; }
	
	public final void setSurfaceConvention(String val) { surfaceConvention = val; }
	
	public final void setMappingMethod(String val) { mappingOption = val; }
	
	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Cortex Processing"; }
	public final String getLabel() { return "Surface Mesh Mapping"; }
	public final String getName() { return "SurfaceMeshMapping"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Map some intensity image on a surface or pair of original+inflated surfaces (use 'Image Label' as suffix for the mapped data)."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
			
	// create outputs
	public final float[] 	getMappedOriginalSurfacePoints() { return mappedOrgSurfacePoints; }
	public final int[] 		getMappedOriginalSurfaceTriangles() { return mappedOrgSurfaceTriangles; }
	public final float[] 	getMappedOriginalSurfaceValues() { return mappedOrgSurfaceValues; }
	
	public final float[] 	getMappedInflatedSurfacePoints() { return mappedInfSurfacePoints; }
	public final int[] 		getMappedInflatedSurfaceTriangles() { return mappedInfSurfaceTriangles; }
	public final float[] 	getMappedInflatedSurfaceValues() { return mappedInfSurfaceValues; }
	
	public void execute() {
		
		System.out.println("Image dimensions: "+nx+" x "+ny+" x "+nz);
		System.out.println("Image resolutions: "+rx+" x "+ry+" x "+rz);
		
		float[][][] intensity3d = null;
		float[][][][] intensity4d = null;
		if (nt==1) {
			intensity3d = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				intensity3d[x][y][z] = intensityImage[x+nx*y+nx*ny*z];
			}
		} else {
			intensity4d = new float[nx][ny][nz][nt];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++) {
				intensity4d[x][y][z][t] = intensityImage[x+nx*y+nx*ny*z+nx*ny*nz*t];
			}
		}
		intensityImage = null;
		
		int npt = origSurfacePoints.length/3;
		int ntr = origSurfaceTriangles.length/3;
		
		// main algorithm
		if (surfaceConvention.equals("mipav")) {
			for (int p=0; p<npt; p++) {
				origSurfacePoints[3*p+X] = origSurfacePoints[3*p+X]/rx;
				origSurfacePoints[3*p+Y] = (ny-1)-origSurfacePoints[3*p+Y]/ry;
				origSurfacePoints[3*p+Z] = (nz-1)-origSurfacePoints[3*p+Z]/rz;
			}
		}

		byte mapStyle = NONE;
		if (mappingOption.equals("closest_point")) mapStyle = CLOSEST;
		else if (mappingOption.equals("linear_interp")) mapStyle = LINEAR;
		else if (mappingOption.equals("highest_value")) mapStyle = HIGHEST;
		
		float[][] data = new float[npt][nt];
		for(int i=0; i<npt; i++){
			float px = origSurfacePoints[3*i+X];
			float py = origSurfacePoints[3*i+Y];
			float pz = origSurfacePoints[3*i+Z];
			
			// mapping from mesh to voxel space: multiple options!
			
			if (mapStyle==CLOSEST) {
				int x = Numerics.round(px);
				int y = Numerics.round(py);
				int z = Numerics.round(pz);
				if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
					//data[i] = new double[]{factor*Numerics.bounded(intensity[x][y][z],2*Imin-Imax,2*Imax-Imin)};
					if (nt>1) {
						for (int t=0;t<nt;t++) data[i][t] = intensity4d[x][y][z][t];
					} else {
						data[i][0] = intensity3d[x][y][z];
					}
				} else {
					if (nt>1) {
						for (int t=0;t<nt;t++) data[i][t] = 0.0f;
					} else {
						data[i][0] = 0.0f;
					}
				}
			} else if (mapStyle==LINEAR) {
				int x = Numerics.round(px);
				int y = Numerics.round(py);
				int z = Numerics.round(pz);
				if (x>0 && x<nx-2 && y>0 && y<ny-2 && z>0 && z<nz-2) {
					float dx = px - x;
					float dy = py - y;
					float dz = pz - z;
					if (nt>1) {
						for (int t=0;t<nt;t++) 						
							data[i][t] = (1-dx)*(1-dy)*(1-dz)*intensity4d[x][y][z][t]
											+dx*(1-dy)*(1-dz)*intensity4d[x+1][y][z][t]
											+(1-dx)*dy*(1-dz)*intensity4d[x][y+1][z][t]
											+(1-dx)*(1-dy)*dz*intensity4d[x][y][z+1][t]
											+dx*dy*(1-dz)*intensity4d[x+1][y+1][z][t]
											+(1-dx)*dy*dz*intensity4d[x][y+1][z+1][t]
											+dx*(1-dy)*dz*intensity4d[x+1][y][z+1][t]
											+dx*dy*dz*intensity4d[x+1][y+1][z+1][t];
					} else {
						data[i][0] = (1-dx)*(1-dy)*(1-dz)*intensity3d[x][y][z]
										+dx*(1-dy)*(1-dz)*intensity3d[x+1][y][z]
										+(1-dx)*dy*(1-dz)*intensity3d[x][y+1][z]
										+(1-dx)*(1-dy)*dz*intensity3d[x][y][z+1]
										+dx*dy*(1-dz)*intensity3d[x+1][y+1][z]
										+(1-dx)*dy*dz*intensity3d[x][y+1][z+1]
										+dx*(1-dy)*dz*intensity3d[x+1][y][z+1]
										+dx*dy*dz*intensity3d[x+1][y+1][z+1];
					}
				} else {
					if (nt>1) {
						for (int t=0;t<nt;t++) data[i][t] = 0.0f;
					} else {
						data[i][0] = 0.0f;
					}
				}				
			} else if (mapStyle==HIGHEST) {
				int x = Numerics.round(px);
				int y = Numerics.round(py);
				int z = Numerics.round(pz);
				if (x>0 && x<nx-2 && y>0 && y<ny-2 && z>0 && z<nz-2) {
					float dx = px - x;
					float dy = py - y;
					float dz = pz - z;
					
					if (nt>1) {
						for (int t=0;t<nt;t++)
							data[i][t] = Numerics.max(intensity4d[x][y][z][t],
														intensity4d[x+1][y][z][t],
														intensity4d[x][y+1][z][t],
														intensity4d[x][y][z+1][t],
														intensity4d[x+1][y+1][z][t],
														intensity4d[x][y+1][z+1][t],
														intensity4d[x+1][y][z+1][t],
														intensity4d[x+1][y+1][z+1][t]);
					} else {
						data[i][0] = Numerics.max(intensity3d[x][y][z],
													intensity3d[x+1][y][z],
													intensity3d[x][y+1][z],
													intensity3d[x][y][z+1],
													intensity3d[x+1][y+1][z],
													intensity3d[x][y+1][z+1],
													intensity3d[x+1][y][z+1],
													intensity3d[x+1][y+1][z+1]);
					}
				} else {
					if (nt>1) {
						for (int t=0;t<nt;t++) data[i][t] = 0.0f;
					} else {
						data[i][0] = 0.0f;
					}
				}				
				
			}
		}
		
		// ouptput: just point to source
		mappedOrgSurfacePoints = origSurfacePoints;
		mappedOrgSurfaceTriangles = origSurfaceTriangles;
		
		mappedInfSurfacePoints = inflatedSurfacePoints;
		mappedInfSurfaceTriangles = inflatedSurfaceTriangles;
		
		mappedOrgSurfaceValues = new float[npt*nt];
		for (int t=0;t<nt;t++) {
			for(int i=0; i<npt; i++) {
				mappedOrgSurfaceValues[i+npt*t] = data[i][t];
			}
		}
		if (mappedInfSurfacePoints!=null) {
			mappedInfSurfaceValues = new float[npt*nt];
			for (int t=0;t<nt;t++) {
				for(int i=0; i<npt; i++) {
					mappedInfSurfaceValues[i+npt*t] = data[i][t];
				}
			}
		}
	}


}
