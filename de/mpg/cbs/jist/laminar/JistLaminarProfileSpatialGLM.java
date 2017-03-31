package de.mpg.cbs.jist.laminar;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;
import Jama.Matrix;
import java.util.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarProfileSpatialGLM extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume layersImage;
	private ParamVolume intensityImage;
	private ParamVolume maskImage;
	private	ParamFloat	paramExtent;
	private	ParamFloat	paramStdev;
		
	private ParamVolume mappedImage;
	private ParamVolume residualImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private static final float HASQ3 = (float)(FastMath.sqrt(3.0)/2.0);

	
	protected void createInputParameters(ParamCollection inputParams) {
		
		//imageParams=new ParamCollection("Images");
		inputParams.add(layersImage = new ParamVolume("Profile Surface Image",null,-1,-1,-1,-1));
		inputParams.add(intensityImage = new ParamVolume("Intensity Image",null,-1,-1,-1,-1));
		inputParams.add(maskImage = new ParamVolume("Cortex Mask (opt)",null,-1,-1,-1,-1));
		maskImage.setMandatory(false);
		inputParams.add(paramStdev = new ParamFloat("Spatial weighting (mm)", 0.0f, 100.0f, 1.0f));
		paramStdev.setDescription("standard deviation of distance to local profile Gaussian prior");
		inputParams.add(paramExtent = new ParamFloat("Spatial extent (mm)", 0.0f, 100.0f, 3.0f));
		paramExtent.setDescription("maximum distance to sample away from local profile (ideally, 3x Spatial Weighting)");
		
		//inputParams.add(imageParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis.devel");
		inputParams.setLabel("Local Spatial GLM Profile");
		inputParams.setName("LocalSpatialGLMProfile");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Sample some intensity image along a cortical profile across layer surfaces assuming a signal mixture (spatial GLM).");
		
		info.setVersion("3.0.9");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(mappedImage = new ParamVolume("Profile-mapped Intensity Image",null,-1,-1,-1,-1));
		
		outputParams.add(residualImage = new ParamVolume("GLM residuals",null,-1,-1,-1,-1));
		
		outputParams.setName("layers GLM images");
		outputParams.setLabel("layers GLM images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays : TO DO
		
		ImageDataFloat	layersImg = new ImageDataFloat(layersImage.getImageData());
		ImageDataFloat	intensImg = new ImageDataFloat(intensityImage.getImageData());
		
		int nx = layersImg.getRows();
		int ny = layersImg.getCols();
		int nz = layersImg.getSlices();
		int nlayers = layersImg.getComponents()-1;
		int nxyz = nx*ny*nz;
		float rx = layersImg.getHeader().getDimResolutions()[0];
		float ry = layersImg.getHeader().getDimResolutions()[1];
		float rz = layersImg.getHeader().getDimResolutions()[2];
		
		float[][] layers = new float[nlayers+1][nxyz];
		float[][][][] buffer4 = layersImg.toArray4d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			layers[l][xyz] = buffer4[x][y][z][l];
		}
		buffer4 = null;
		layersImg = null;
		
		float[] intensity = new float[nxyz];
		float[][][] buffer3 = intensImg.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			intensity[xyz] = buffer3[x][y][z];
		}
		buffer3 = null;
		intensImg = null;
		
		// create a mask for all the regions outside of the area where layer 1 is > 0 and layer 2 is < 0
		boolean[] ctxmask = new boolean[nxyz];
		if (maskImage.getImageData()!=null) {
			ImageDataUByte	maskImg = new ImageDataUByte(maskImage.getImageData());
			byte[][][] bufferbyte = maskImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				ctxmask[xyz] = (layers[0][xyz]>=0.0 && layers[nlayers][xyz]<=0.0 && bufferbyte[x][y][z]>0);
			}
			bufferbyte = null;
			maskImg = null;
		} else {
			for (int xyz=0;xyz<nxyz;xyz++) {
				ctxmask[xyz] = (layers[0][xyz]>=0.0 && layers[nlayers][xyz]<=0.0);
			}
		}
				
		// main algorithm
		
		// 1. define partial voume for each layer, each voxel
		float[][] pvol = new float[nlayers+1][nxyz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				for (int l=0;l<=nlayers;l++) {
					pvol[l][xyz] = partialVolumeFromSurface(x, y, z, layers[l], nx, ny, nz);
				}
			}
		}
		
		// 2. build and invert the GLM for each profile / voxel?
		float delta = paramExtent.getValue().floatValue()/Numerics.min(rx,ry,rz);
		float stdev = paramStdev.getValue().floatValue()/Numerics.min(rx,ry,rz);
		
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		float maskval = 1e13f;
		float[][][][] mapping = new float[nx][ny][nz][nlayers];
		float[][][] residual = new float[nx][ny][nz];
		BitSet sampled = new BitSet(nx*ny*nz);
		float[] profiledist = new float[nx*ny*nz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				findFastMarchingProfileNeighborhood(sampled, profiledist, x,y,z, delta, layers, profile, ctxmask, nx, ny, nz, nlayers);

				int nsample = sampled.size();
				if (nsample>=nlayers) {
					double[][] glm = new double[nlayers][nsample];
					double[][] data = new double[nsample][1];
					int idx = 0;
					for (int n=0;n<nsample;n++) {
						// get the next non-zero value
						idx = sampled.nextSetBit(idx);
						// build a weighting function based on distance to the original location
						double weight = FastMath.exp(-0.5*Numerics.square(profiledist[idx]/stdev));
						
						for (int l=0;l<nlayers;l++) {
							glm[l][n] = weight*(pvol[l+1][idx]-pvol[l][idx]);
						}
						data[n][0] = weight*intensity[idx];
					}
					// invert the linear model
					Matrix mtx = new Matrix(glm);
					Matrix smp = new Matrix(data);
					Matrix val = mtx.solve(smp);
					for (int l=0;l<nlayers;l++) {
						mapping[x][y][z][l] = (float)val.get(l,0);
					}
					Matrix res = mtx.times(val).minus(smp);
					residual[x][y][z] = (float)res.normInf();
				}
			}
		}
		
		// output
		String imgname = intensityImage.getImageData().getName();
		
		ImageDataFloat mapData = new ImageDataFloat(mapping);		
		mapData.setHeader(layersImage.getImageData().getHeader());
		mapData.setName(imgname+"_glmprofiles");
		mappedImage.setValue(mapData);
		mapData = null;
		mapping = null;
		
		ImageDataFloat resData = new ImageDataFloat(residual);		
		resData.setHeader(layersImage.getImageData().getHeader());
		resData.setName(imgname+"_glmresidual");
		residualImage.setValue(resData);
		resData = null;
		residual = null;
	}

	// brute force
	void findProfileNeighborhood(BitSet sampled, int x, int y, int z, int delta, float[][] layers, CorticalProfile profile, boolean[] ctxmask, int nx, int ny, int nz, int nlayers) {
		sampled.clear();
		
		// compute profiles around a certain size
		for (int dx=-delta;dx<=delta;dx++) for (int dy=-delta;dy<=delta;dy++) for (int dz=-delta;dz<=delta;dz++) {
			profile.computeTrajectory(layers, x+dx, y+dy, z+dz);
		
			// add the underlying voxels to the estimation matrix
			for (int l=0;l<=nlayers;l++) {
				int xyzs = Numerics.round(profile.getPt(l)[X]) + nx*Numerics.round(profile.getPt(l)[Y]) + nx*ny*Numerics.round(profile.getPt(l)[Z]);
				if (ctxmask[xyzs]) sampled.set(xyzs);
			}
		}
	}
	
	// should be faster and more precise
	void findFastMarchingProfileNeighborhood(BitSet sampled, float[] dist, int x, int y, int z, float delta, float[][] layers, CorticalProfile profile, boolean[] ctxmask, int nx, int ny, int nz, int nlayers) {
		sampled.clear();
		
		// compute profile centered on point of interest
		profile.computeTrajectory(layers, x, y, z);
		
		// add the underlying voxels to the estimation matrix
		for (int l=0;l<=nlayers;l++) {
			int xyzs = Numerics.round(profile.getPt(l)[X]) + nx*Numerics.round(profile.getPt(l)[Y]) + nx*ny*Numerics.round(profile.getPt(l)[Z]);
			if (ctxmask[xyzs]) sampled.set(xyzs);
		}
		
		// find enighbors through fast marching
		BinaryHeap2D heap = new BinaryHeap2D(Numerics.ceil(2*delta)*sampled.size(), BinaryHeap2D.MINTREE);
		int idx = 0;
		for (int n=0;n<sampled.size();n++) {
			// get the next non-zero value
			idx = sampled.nextSetBit(idx);
			heap.addValue(0.0f,idx,(byte)1);
		}
		// reset the sampled values: the original ones will be re-added from the heap
		sampled.clear();
		// main loop
		while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float curdist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	heap.removeFirst();
        	
        	// for points selected multiple times
        	if (sampled.get(xyz)) continue;
        	
        	// add to sampled values
        	sampled.set(xyz);
        	dist[xyz] = curdist;
        	
        	for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) {
				int ngb = xyz+dx+nx*dy+nx*ny*dz;
				if (!sampled.get(ngb) && ctxmask[ngb]) {
					// add to the tree if close enough from first profile
					float newdist = curdist+(float)FastMath.sqrt(dx*dx+dy*dy+dz*dz);
					if (newdist<=delta) heap.addValue(newdist, ngb, (byte)1);
				}
			}
		}
	}

	float partialVolumeFromSurface(int x0, int y0, int z0, float[] levelset, int nx, int ny, int nz) {
		int xyz0 = x0+nx*y0+nx*ny*z0;
		double phi0 = levelset[xyz0];
		
		if (phi0<-HASQ3) return 1.0f;
		else if (phi0>HASQ3) return 0.0f;
		
		// use direction, distance to create a continuous plane model V(X-X0)+phi(X0) = d(X,P)
		double vx = 0.5*(levelset[xyz0+1]-levelset[xyz0-1]);
		double vy = 0.5*(levelset[xyz0+nx]-levelset[xyz0-nx]);
		double vz = 0.5*(levelset[xyz0+nx*ny]-levelset[xyz0-nx*ny]);
		double norm = FastMath.sqrt(vx*vx+vy*vy+vz*vz);
		
		// integrate the voxel cube over that model
		
		// first shot: numerically (other options: Heaviside integration, case-by-case geometric solution?)
		double volume = 0.0;
		for (double x=-0.45; x<=0.45; x+=0.1) for (double y=-0.45; y<=0.45; y+=0.1) for (double z=-0.45; z<=0.45; z+=0.1) {
			if (vx*x+vy*y+vz*z+norm*phi0<=0) volume+=0.001;
		}
		
		// second option: with explicit definition of the intersection (~marching cubes!)
		// requires intersecting all 12 edges from the voxel cube, then defining either case-by-case formulas or 
		// deriving the formula for an integral of the Heaviside function
		
		return (float)volume;
	}

	float fastApproxPartialVolumeFromSurface(int x0, int y0, int z0, float[] levelset, int nx, int ny, int nz) {
		return Numerics.bounded(0.5f-levelset[x0+nx*y0+nx*ny*z0],0.0f,1.0f);
	}	
}
