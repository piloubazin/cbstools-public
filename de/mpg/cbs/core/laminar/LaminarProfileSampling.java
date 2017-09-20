package de.mpg.cbs.core.laminar;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class LaminarProfileSampling {

	private float[] layersImage;
	private float[] intensityImage;
	private byte[] maskImage=null;
	private String interpParam="linear";	
	
	private int nx, ny, nz, nt, nxyz;
	private float rx, ry, rz;

	private float[] mappedImage;
	private byte[] mappedmaskImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	// create inputs
	public final void setProfileSurfaceImage(float[] val) { layersImage = val; }
	public final void setIntensityImage(float[] val) { intensityImage = val; }
	public final void setCortexMask(byte[] val) { maskImage = val; }
	public final void setInterpolation(String val) { interpParam = val; }
	
	public final void setDimensions(int x, int y, int z, int t) { nx=x; ny=y; nz=z; nt=t; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nt=dim[3]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Laminar Analysis"; }
	public final String getLabel() { return "Profile Sampling"; }
	public final String getName() { return "ProfileSampling"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin","Juliane Dinse"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Sample some intensity image along a cortical profile across layer surfaces."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
			
	// create outputs
	public final float[] getProfileMappedIntensityImage() { return mappedImage; }
	public final byte[] getProfile4Dmask() { return mappedmaskImage; }
	
	public void execute(){
		
		int nlayers = nt-1;
		
		float[][] layers = new float[nlayers+1][nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			layers[l][xyz] = layersImage[xyz+nxyz*l];
		}
		layersImage = null;
		
		float[] intensity = intensityImage;
		
		// create a mask for all the regions outside of the area where layer 1 is > 0 and layer 2 is < 0
		boolean[] ctxmask = new boolean[nxyz];
		if (maskImage!=null) {
			for (int xyz=0;xyz<nxyz;xyz++) {
				ctxmask[xyz] = (maskImage[xyz]>0);
			}
		} else {
			for (int xyz=0;xyz<nxyz;xyz++) {
				ctxmask[xyz] = (layers[0][xyz]>=0.0 && layers[nlayers][xyz]<=0.0);
			}
		}
				
		// main algorithm
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		byte LINEAR = 1;
		byte NEAREST = 2;
		byte interp = LINEAR;
		if (interpParam.equals("nearest")) interp = NEAREST;
		
		float maskval = 1e13f;
		float[][][][] mapping = new float[nx][ny][nz][nlayers+1];
		byte[][][][] mappingmask = new byte[nx][ny][nz][nlayers+1];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				profile.computeTrajectory(layers, x, y, z);
				
				for (int l=0;l<=nlayers;l++) {
					// interpolate the contrast
					if (interp==NEAREST) {
						mapping[x][y][z][l] = ImageInterpolation.nearestNeighborInterpolation(intensity, ctxmask, maskval, 
																					profile.getPt(l)[X], profile.getPt(l)[Y], profile.getPt(l)[Z], 
																					nx, ny, nz);
					} else {
						mapping[x][y][z][l] = ImageInterpolation.linearInterpolation(intensity, ctxmask, maskval, 
																					profile.getPt(l)[X], profile.getPt(l)[Y], profile.getPt(l)[Z], 
																					nx, ny, nz);
					}
					if (mapping[x][y][z][l]==maskval) {
						mappingmask[x][y][z][l] = (byte)0;
						mapping[x][y][z][l] = 0.0f;
					} else {
						mappingmask[x][y][z][l] = (byte)1;
					}
				}
			}
		}
		
		// output
		mappedImage = new float[nx*ny*nz*(nlayers+1)];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			mappedImage[xyz+nxyz*l] = mapping[x][y][z][l];
		}
		
		mappedmaskImage = new byte[nx*ny*nz*(nlayers+1)];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			mappedmaskImage[xyz+nxyz*l] = mappingmask[x][y][z][l];
		}
	}


}
