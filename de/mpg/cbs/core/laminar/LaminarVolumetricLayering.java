package de.mpg.cbs.core.laminar;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class LaminarVolumetricLayering {

	private float[] gwImage;
	private float[] cgImage;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private int 	layersParam = 10;
	private int 	iterationParamNarrowBand = 500;
	private float 	minimumParamNarrowBand = 0.0005f;
	private String 	algoParam = "volume-preserving";
	private String 	dirParam = "outward";
	private String 	topologyParam = "no";
	private int		kernelParam = 3;
	private float 	ratioKernelParam = 1.0f;
	private	boolean	presmoothParam = false;
	
	private static final String[] algoTypes = {"distance-preserving", "volume-preserving"};
	private static final String[] dirTypes = {"outward", "inward"};
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	private String	lutdir = null;
	
	private float[] layeringImage;
	private byte[] labelsImage;
	private float[] surfImage;
	private float[] midsurfImage;
	private int surfImageLength;
	private int midsurfImageLength;
	
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;
	
	private static final boolean debugOutput=false;
	
	// create inputs
	public final void setInnerDistanceImage(float[] val) { gwImage = val; }
	public final void setOuterDistanceImage(float[] val) { cgImage = val; }

	public final void setNumberOfLayers(int val) { layersParam = val; }		
	public final void setMaxNarrowBandIterations(int val) { iterationParamNarrowBand = val; } 
	public final void setMinNarrowBandChange(float val) { minimumParamNarrowBand = val; } 
	public final void setLayeringMethod(String val) { algoParam = val; }
	public final void setLayeringDirection(String val) { dirParam = val; }
	public final void setCurvatureApproximationScale(int val) { kernelParam = val; }
	public final void setRatioSmoothingKernelSize(float val) { ratioKernelParam = val; }
	public final void setPresmoothCorticalSurfaces(boolean val) { presmoothParam = val; }

	public final void setTopology(String val) { topologyParam = val; }
	public final void setTopologyLUTdirectory(String val) { lutdir = val; }

			
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public  final String getPackage() { return "CBS Tools"; }
	public  final String getCategory() { return "Laminar Analysis"; }
	public  final String getLabel() { return "Volumetric Layering"; }
	public  final String getName() { return "VolumetricLayering"; }
	public  final String getVersion() { return "3.1"; }

	public  final String[] getAlgorithmAuthors() {return new String[] {"Miriam Waehnert","Pierre-Louis Bazin","Juliane Dinse"}; }
	public  final String getCitation() { return "Waehnert MD, Dinse J, Weiss M, Streicher MN, Waehnert P, Geyer S, Turner R, Bazin PL, "
												+"Anatomically motivated modeling of cortical laminae, "
												+"Neuroimage, 2013."; }
	public  final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public  final String getDescription() { return "Builds a continuous layering of the cortex following distance-preserving or volume-preserving models of cortical folding."; }
	public  final String getLongDescription() { return getDescription(); }

	// create outputs
	public final float[] getContinuousDepthMeasurement() { return layeringImage; }
	public final byte[] getDiscreteSampledLayers() { return labelsImage; }
	public final float[] getLayerBoundarySurfaces() { return surfImage; }
	public final float[] getLayerCenteredSurfaces() { return midsurfImage; }
	public final int getLayerBoundarySurfacesLength() { return surfImageLength; }
	public final int getLayerCenteredSurfacesLength() { return midsurfImageLength; }
	
		
	public void execute(){
		
		// import the image data into 1D arrays
		float[] inner = gwImage;
		float[] outer = cgImage;
		
		// segmentation masks : directly from data
		// for WM only 
		byte[] init = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			//if (bytebuffer[x][y][z]==2.0f ) init[xyz] = 1;
			if (inner[xyz]<0.0f ) init[xyz] = 1;
			else init[xyz] = 0;
		}

		// for WM and GM combined
		byte[] target = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			//if (bytebuffer[x][y][z]> 0.0f) target[xyz] = 1;
			if (outer[xyz]<0.0f) target[xyz] = 1;
			else target[xyz] = 0;
		}
		//bytebuffer = null;
	
		// for the cortex only
		boolean[] cortex = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			cortex[xyz] = (target[xyz]>0 && init[xyz]==0);
		}	

		BasicInfo.displayMessage("initialize levelset functions and computation area...\n");
		
		// create a mask for all the regions outside of the computation area
		// (add an area of 10 mm around the cortex)
		// also create masks of around boundaries for curvature computations
		boolean[] mask = new boolean[nxyz];
		boolean[] innermask = new boolean[nxyz];
		boolean[] outermask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			mask[xyz] = ( (target[xyz]>0 && init[xyz]==0) || (Numerics.abs(inner[xyz])<=10.0/rx || Numerics.abs(outer[xyz])<=10.0/rx) );
			innermask[xyz] = (Numerics.abs(inner[xyz])<=1.733f);
			outermask[xyz] = (Numerics.abs(outer[xyz])<=1.733f);
		}	
		
		
		//// Step 1: build smooth distance functions from the input (optional) ////
		if (presmoothParam) {
			BasicInfo.displayMessage("smoothing of initial boundaries\n");
			
			// smooth the levelsets by default (remove stitching artefact at the boundary)
			SmoothGdm gdmin = new SmoothGdm(inner, nx, ny, nz, rx, ry, rz,
												mask, 0.9f, 0.1f, topologyParam, lutdir);
			
			gdmin.evolveNarrowBand(iterationParamNarrowBand, minimumParamNarrowBand);
			inner = gdmin.exportLevelset();
			gdmin.finalize();
			gdmin = null;
			
			SmoothGdm gdmout = new SmoothGdm(outer, nx, ny, nz, rx, ry, rz,
												mask, 0.9f, 0.1f, topologyParam, lutdir);
			
			gdmout.evolveNarrowBand(iterationParamNarrowBand, minimumParamNarrowBand);
			outer = gdmout.exportLevelset();
			gdmout.finalize();
			gdmout = null;
		}
		
		//// Step 2: precompute curvatures with quadric approximation (for volume-preserving only) ////
		float[][] curvin = null, curvout = null;
		// compute inside the whole computation area (not just the cortex)
		if (algoParam.equals("volume-preserving")) {
			BasicInfo.displayMessage("inner curvature ("+kernelParam+")...\n");
			curvin = ImageGeometry.quadricCurvatureEstimates(inner, innermask, mask, kernelParam, nx,ny,nz);
			//curvin = ImageGeometry.principalCurvatureDirections(inner, mask, nx,ny,nz);
			
			BasicInfo.displayMessage("outer curvature ("+kernelParam+")...\n");
			curvout = ImageGeometry.quadricCurvatureEstimates(outer, outermask, mask, kernelParam, nx,ny,nz);
			//curvout = ImageGeometry.principalCurvatureDirections(outer, mask, nx,ny,nz);
		}
		
		//// Step 3: run the algorithm distance-based layering ////
		BasicInfo.displayMessage("distance-preserving evolution\n");
		
		int Nlayers = layersParam;
		float[][] layers = new float[Nlayers+1][];
		layers[0] = inner;
		layers[Nlayers] = outer;
		String modelType = algoParam;
		
		VolumetricLayeringGdm gdm = new VolumetricLayeringGdm(inner, outer, "distance-preserving", dirParam, 
																0.5f, null, null, 1.0f,
																nx, ny, nz, rx, ry, rz,
																mask, 0.9f, 0.1f, topologyParam, lutdir);
		
		if (dirParam.equals("outward")) {
			for (int t=1;t<Nlayers;t++) {
				BasicInfo.displayMessage(t+"-th layer estimation...\n");
				gdm.setFraction((float)t/(float)Nlayers);
				gdm.evolveNarrowBand(iterationParamNarrowBand, minimumParamNarrowBand);
				layers[t] = gdm.exportLevelset();
			}
		} else if (dirParam.equals("inward")) {
			for (int t=Nlayers-1;t>0;t--) {
				BasicInfo.displayMessage(t+"-th layer estimation...\n");
				gdm.setFraction((float)t/(float)Nlayers);
				gdm.evolveNarrowBand(iterationParamNarrowBand, minimumParamNarrowBand);
				layers[t] = gdm.exportLevelset();
			}
		}
		
		//// Step 4: compute profiles to get location of first, last points and length ////
		float[][] volumein = null, volumeout = null;
		if (algoParam.equals("volume-preserving")) {
		
			BasicInfo.displayMessage("volume-preserving modeling\n");
		
			// store the result in corresponding variables?
			volumein = new float[3][nxyz];
			volumeout = new float[3][nxyz];
			
			CorticalProfile profile = new CorticalProfile(Nlayers, nx, ny, nz, rx, ry, rz);
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				
				if (cortex[xyz]) {			
					// 1. Build the profile
					profile.computeTrajectory(layers, x, y, z);
					
					volumein[0][xyz] = profile.computeLength();
					volumeout[0][xyz] = volumein[0][xyz];
					
					// 2. get inner curvature values; no directions needed
					
					volumein[1][xyz] = ImageInterpolation.linearInterpolation(curvin[0], mask, 0.0f, 
																				profile.getPt(0)[X], profile.getPt(0)[Y], profile.getPt(0)[Z],
																				nx,ny,nz);
					volumein[2][xyz] = ImageInterpolation.linearInterpolation(curvin[4], mask, 0.0f, 
																				profile.getPt(0)[X], profile.getPt(0)[Y], profile.getPt(0)[Z],
																				nx,ny,nz);
						
					// 3. get outer curvature values
					
					volumeout[1][xyz] = ImageInterpolation.linearInterpolation(curvout[0], mask, 0.0f, 
																				 profile.getPt(Nlayers)[X], profile.getPt(Nlayers)[Y], profile.getPt(Nlayers)[Z],
																				 nx,ny,nz);
					volumeout[2][xyz] = ImageInterpolation.linearInterpolation(curvout[4], mask, 0.0f, 
																				 profile.getPt(Nlayers)[X], profile.getPt(Nlayers)[Y], profile.getPt(Nlayers)[Z],
																				 nx,ny,nz);	
				}
			}
		} else
		if (algoParam.equals("volume-preserving2")) {
		
			BasicInfo.displayMessage("volume-preserving modeling (approx area model)\n");
		
			// store the result in corresponding variables?
			volumein = new float[1][nxyz];
			volumeout = new float[1][nxyz];
			
			CorticalProfile profile = new CorticalProfile(Nlayers, nx, ny, nz, rx, ry, rz);
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				
				if (cortex[xyz]) {			
					// 1. Build the profile
					profile.computeTrajectory(layers, x, y, z);
					float xi = profile.getPt(0)[X];
					float yi = profile.getPt(0)[Y];
					float zi = profile.getPt(0)[Z];
					float xo = profile.getPt(Nlayers)[X];
					float yo = profile.getPt(Nlayers)[Y];
					float zo = profile.getPt(Nlayers)[Z];
					
					//if (profile.checkCoordinates())
					//	System.out.println("->NaN (1) <"+x+", "+y+", "+z+">");
					
					// 2. Get points at +/-1 at (x,y,z)
					float[] dir = profile.computeTangentAt(x, y, z);
					int u = Numerics.argminmag(dir[X], dir[Y], dir[Z]);
					int v = Numerics.argsecmag(dir[X], dir[Y], dir[Z]);
					double[] du = new double[3];
					if (u==X) {
						du[X] = 1.0-dir[X]*dir[X];
						du[Y] =    -dir[X]*dir[Y];
						du[Z] =    -dir[X]*dir[Z];
					} else if (u==Y) {
						du[X] =    -dir[Y]*dir[X];
						du[Y] = 1.0-dir[Y]*dir[Y];
						du[Z] =    -dir[Y]*dir[Z];
					} else if (u==Z) {
						du[X] =    -dir[Z]*dir[X];
						du[Y] =    -dir[Z]*dir[Y];
						du[Z] = 1.0-dir[Z]*dir[Z];
					}
					double ndu = FastMath.sqrt(du[X]*du[X]+du[Y]*du[Y]+du[Z]*du[Z]);
					if (ndu>0.001) {
						du[X] /= ndu; du[Y] /= ndu; du[Z] /= ndu;
					}
					
					double[] dv = new double[3];
					if (v==X) {
						dv[X] = 1.0-dir[X]*dir[X];
						dv[Y] =    -dir[X]*dir[Y];
						dv[Z] =    -dir[X]*dir[Z];
					} else if (v==Y) {
						dv[X] =    -dir[Y]*dir[X];
						dv[Y] = 1.0-dir[Y]*dir[Y];
						dv[Z] =    -dir[Y]*dir[Z];
					} else if (v==Z) {
						dv[X] =    -dir[Z]*dir[X];
						dv[Y] =    -dir[Z]*dir[Y];
						dv[Z] = 1.0-dir[Z]*dir[Z];
					}
					double ndv = FastMath.sqrt(dv[X]*dv[X]+dv[Y]*dv[Y]+dv[Z]*dv[Z]);
					if (ndv>0.001) {
						dv[X] /= ndv; dv[Y] /= ndv; dv[Z] /= ndv;
					}
					
					// 3. reconstruct each profile
					profile.computeTrajectory(layers, (float)(x+du[X]), (float)(y+du[Y]), (float)(z+du[Z]));
					float xiu1 = profile.getPt(0)[X];
					float yiu1 = profile.getPt(0)[Y];
					float ziu1 = profile.getPt(0)[Z];
					float xou1 = profile.getPt(Nlayers)[X];
					float you1 = profile.getPt(Nlayers)[Y];
					float zou1 = profile.getPt(Nlayers)[Z];
					
					//if (profile.checkCoordinates())
					//	System.out.println("->NaN (2) <"+x+", "+y+", "+z+"> + <"+du[X]+", "+du[Y]+", "+du[Z]+">");
				
					profile.computeTrajectory(layers, (float)(x-du[X]), (float)(y-du[Y]), (float)(z-du[Z]));
					float xiu2 = profile.getPt(0)[X];
					float yiu2 = profile.getPt(0)[Y];
					float ziu2 = profile.getPt(0)[Z];
					float xou2 = profile.getPt(Nlayers)[X];
					float you2 = profile.getPt(Nlayers)[Y];
					float zou2 = profile.getPt(Nlayers)[Z];
					
					//if (profile.checkCoordinates())
					//	System.out.println("->NaN (3) <"+x+", "+y+", "+z+"> - <"+du[X]+", "+du[Y]+", "+du[Z]+">");
					
					profile.computeTrajectory(layers, (float)(x+dv[X]), (float)(y+dv[Y]), (float)(z+dv[Z]));
					float xiv1 = profile.getPt(0)[X];
					float yiv1 = profile.getPt(0)[Y];
					float ziv1 = profile.getPt(0)[Z];
					float xov1 = profile.getPt(Nlayers)[X];
					float yov1 = profile.getPt(Nlayers)[Y];
					float zov1 = profile.getPt(Nlayers)[Z];
					
					//if (profile.checkCoordinates())
					//	System.out.println("->NaN (4) <"+x+", "+y+", "+z+"> + <"+dv[X]+", "+dv[Y]+", "+dv[Z]+">");
					
					profile.computeTrajectory(layers, (float)(x-dv[X]), (float)(y-dv[Y]), (float)(z-dv[Z]));
					float xiv2 = profile.getPt(0)[X];
					float yiv2 = profile.getPt(0)[Y];
					float ziv2 = profile.getPt(0)[Z];
					float xov2 = profile.getPt(Nlayers)[X];
					float yov2 = profile.getPt(Nlayers)[Y];
					float zov2 = profile.getPt(Nlayers)[Z];
					
					//if (profile.checkCoordinates())
					//	System.out.println("->NaN (5) <"+x+", "+y+", "+z+"> - <"+dv[X]+", "+dv[Y]+", "+dv[Z]+">");
					
					// 4. get areas
					volumein[0][xyz] = (float)(FastMath.sqrt( (xiu1-xi)*(xiu1-xi)+(yiu1-yi)*(yiu1-yi)+(ziu1-zi)*(ziu1-zi) )
												*FastMath.sqrt( (xiv1-xi)*(xiv1-xi)+(yiv1-yi)*(yiv1-yi)+(ziv1-zi)*(ziv1-zi) )
												+FastMath.sqrt( (xiu2-xi)*(xiu2-xi)+(yiu2-yi)*(yiu2-yi)+(ziu2-zi)*(ziu2-zi) )
												*FastMath.sqrt( (xiv2-xi)*(xiv2-xi)+(yiv2-yi)*(yiv2-yi)+(ziv2-zi)*(ziv2-zi) )
												+FastMath.sqrt( (xiu1-xi)*(xiu1-xi)+(yiu1-yi)*(yiu1-yi)+(ziu1-zi)*(ziu1-zi) )
												*FastMath.sqrt( (xiv2-xi)*(xiv2-xi)+(yiv2-yi)*(yiv2-yi)+(ziv2-zi)*(ziv2-zi) )
												+FastMath.sqrt( (xiu2-xi)*(xiu2-xi)+(yiu2-yi)*(yiu2-yi)+(ziu2-zi)*(ziu2-zi) )
												*FastMath.sqrt( (xiv1-xi)*(xiv1-xi)+(yiv1-yi)*(yiv1-yi)+(ziv1-zi)*(ziv1-zi) ) );
									  
					volumeout[0][xyz] = (float)(FastMath.sqrt( (xou1-xo)*(xou1-xo)+(you1-yo)*(you1-yo)+(zou1-zo)*(zou1-zo) )
												*FastMath.sqrt( (xov1-xo)*(xov1-xo)+(yov1-yo)*(yov1-yo)+(zov1-zo)*(zov1-zo) )
												+FastMath.sqrt( (xou2-xo)*(xou2-xo)+(you2-yo)*(you2-yo)+(zou2-zo)*(zou2-zo) )
												*FastMath.sqrt( (xov2-xo)*(xov2-xo)+(yov2-yo)*(yov2-yo)+(zov2-zo)*(zov2-zo) )
												+FastMath.sqrt( (xou1-xo)*(xou1-xo)+(you1-yo)*(you1-yo)+(zou1-zo)*(zou1-zo) )
												*FastMath.sqrt( (xov2-xo)*(xov2-xo)+(yov2-yo)*(yov2-yo)+(zov2-zo)*(zov2-zo) )
												+FastMath.sqrt( (xou2-xo)*(xou2-xo)+(you2-yo)*(you2-yo)+(zou2-zo)*(zou2-zo) )
												*FastMath.sqrt( (xov1-xo)*(xov1-xo)+(yov1-yo)*(yov1-yo)+(zov1-zo)*(zov1-zo) ) );
									  
				}
			}
		}


		//// Step 5: recompute the layers with the volume-preserving model ////
		BasicInfo.displayMessage("volume-preserving evolution\n");
		
		if (algoParam.equals("volume-preserving")) {
			gdm = new VolumetricLayeringGdm(inner, outer, "volume-preserving", dirParam, 0.5f, 
										volumein, volumeout, 1.0f,
										nx, ny, nz, rx, ry, rz,
										mask, 0.9f, 0.1f, topologyParam, lutdir);
			
			if (dirParam.equals("outward")) {
				for (int t=1;t<Nlayers;t++) {
					BasicInfo.displayMessage(t+"-th layer estimation...\n");
					
					gdm.setFraction((float)t/(float)Nlayers);
					// compute the ratio also outside the cortex to ensure good boundary behavior
					gdm.computeVolumetricRatio(ratioKernelParam, mask);	
					
					gdm.evolveNarrowBand(iterationParamNarrowBand, minimumParamNarrowBand);
					layers[t] = gdm.exportLevelset();
				}
			} else if (dirParam.equals("inward")) {
				for (int t=Nlayers-1;t>0;t--) {
					BasicInfo.displayMessage(t+"-th layer estimation...\n");
					
					gdm.setFraction((float)t/(float)Nlayers);
					// compute the ratio also outside the cortex to ensure good boundary behavior
					gdm.computeVolumetricRatio(ratioKernelParam, mask);	
					
					gdm.evolveNarrowBand(iterationParamNarrowBand, minimumParamNarrowBand);
					layers[t] = gdm.exportLevelset();
				}
			}
		} else
		if (algoParam.equals("volume-preserving2")) {
			gdm = new VolumetricLayeringGdm(inner, outer, "volume-preserving2", dirParam, 0.5f, 
										volumein, volumeout, 1.0f,
										nx, ny, nz, rx, ry, rz,
										mask, 0.9f, 0.1f, topologyParam, lutdir);
			
			if (dirParam.equals("outward")) {
				for (int t=1;t<Nlayers;t++) {
					BasicInfo.displayMessage(t+"-th layer estimation...\n");
					
					gdm.setFraction((float)t/(float)Nlayers);
					// compute the ratio also outside the cortex to ensure good boundary behavior
					gdm.computeVolumetricRatio2(ratioKernelParam, mask);	
					
					gdm.evolveNarrowBand(iterationParamNarrowBand, minimumParamNarrowBand);
					layers[t] = gdm.exportLevelset();
				}
			} else if (dirParam.equals("inward")) {
				for (int t=Nlayers-1;t>0;t--) {
					BasicInfo.displayMessage(t+"-th layer estimation...\n");
					
					gdm.setFraction((float)t/(float)Nlayers);
					// compute the ratio also outside the cortex to ensure good boundary behavior
					gdm.computeVolumetricRatio2(ratioKernelParam, mask);	
					
					gdm.evolveNarrowBand(iterationParamNarrowBand, minimumParamNarrowBand);
					layers[t] = gdm.exportLevelset();
				}
			}
		}
		// to be added: iterate Steps 4 and 5 until convergence (to be defined)
		
		
		//// Step 6: build a partial volume and continuous layering model ////
		BasicInfo.displayMessage("generate outputs\n");
		
		// layer value interpolation
		float[] layering = new float[nx*ny*nz];
		byte[] labels = new byte[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (layers[0][xyz]>=0 && layers[Nlayers][xyz]<=0) {
				// find the closest
				int best=0;
				for (int t=1;t<=Nlayers;t++) {
					if (Numerics.abs(layers[t][xyz])<Numerics.abs(layers[best][xyz])) {
						best = t;
					}
				}
				// find the next
				int second=-1;
				boolean inside = true;
				if (layers[best][xyz]>=0) {
					// outside : the second best is next
					if (best==Nlayers) inside = false;
					else second = best+1;
				} else {
					// inside : the second best is previous
					if (best==0) inside = false;
					else second = best-1;
				}
				if (inside) {
				
					// interpolate
					layering[xyz] = ( ((float)best/(float)Nlayers)*Numerics.abs(layers[second][xyz])
										 + ((float)second/(float)Nlayers)*Numerics.abs(layers[best][xyz]) )
										 / Numerics.max(1e-6f, Numerics.abs(layers[best][xyz])+Numerics.abs(layers[second][xyz]));
					
					labels[xyz] = (byte)Numerics.max(best,second);
				
				} else {
					labels[xyz] = 0;
				}
			}
		}
		
		// output
		BasicInfo.displayMessage("...layering\n");
		
		layeringImage = layering;
		
		BasicInfo.displayMessage("...labels\n");
		
		labelsImage = labels;
				
		BasicInfo.displayMessage("...surfaces\n");
		
		float[] layersurf = new float[nx*ny*nz*(Nlayers+1)];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int n=0;n<=Nlayers;n++) {
			int xyz = x+nx*y+nx*ny*z;
			layersurf[xyz+nx*ny*nz*n] = layers[n][xyz];
		}
		surfImage = layersurf;
		surfImageLength = Nlayers+1;
		
		// useful or confusing?? useful for traditional profile sampling
		float[] midlayersurf = new float[nx*ny*nz*Nlayers];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int n=0;n<Nlayers;n++) {
			int xyz = x+nx*y+nx*ny*z;
			midlayersurf[xyz+nx*ny*nz*n] = 0.5f*(layers[n][xyz]+layers[n+1][xyz]);
		}
		midsurfImage = midlayersurf;
		midsurfImageLength = Nlayers;
		
		BasicInfo.displayMessage("done\n");
		
	}


}
