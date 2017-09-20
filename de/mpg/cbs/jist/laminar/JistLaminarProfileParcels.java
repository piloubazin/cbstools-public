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
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
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
import org.apache.commons.math3.random.*;
import java.util.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarProfileParcels extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume layersImage;
	private ParamVolume dataImage;
	private ParamVolume maskImage;
	private	ParamFloat	paramExtent;
	private	ParamFloat	paramStdev;
		
	private ParamVolume mappedImage;
	private ParamVolume distanceImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private static final float HASQ3 = (float)(FastMath.sqrt(3.0)/2.0);

	
	protected void createInputParameters(ParamCollection inputParams) {
		
		//imageParams=new ParamCollection("Images");
		inputParams.add(layersImage = new ParamVolume("Profile Surface Image",null,-1,-1,-1,-1));
		inputParams.add(maskImage = new ParamVolume("Cortex Mask (opt)",null,-1,-1,-1,-1));
		maskImage.setMandatory(false);
		inputParams.add(paramStdev = new ParamFloat("Spatial weighting (mm)", 0.0f, 100.0f, 1.0f));
		paramStdev.setDescription("standard deviation of distance to local profile Gaussian prior");
		inputParams.add(paramExtent = new ParamFloat("Spatial extent (mm)", 0.0f, 100.0f, 3.0f));
		paramExtent.setDescription("maximum distance to sample away from local profile (ideally, 3x Spatial Weighting)");
		
		//inputParams.add(imageParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis.devel");
		inputParams.setLabel("Define Cortical Parcels");
		inputParams.setName("DefineCorticalParcels");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Create random small cortical parcels of a given size.");
		
		info.setVersion("3.0.9");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(mappedImage = new ParamVolume("Parcel Image",null,-1,-1,-1,-1));
		outputParams.add(distanceImage = new ParamVolume("Distance Image",null,-1,-1,-1,-1));
		
		outputParams.setName("parcel images");
		outputParams.setLabel("parcel images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays : TO DO
		System.out.println("Import data");
		
		String name = Interface.getName(layersImage);
		ImageHeader header = Interface.getHeader(layersImage);
		
		int[] dims = Interface.getDimensions(layersImage);
		float[] res = Interface.getResolutions(layersImage);
		int nlayers = Interface.getComponents(layersImage)-1;
		int nxyz = dims[X]*dims[Y]*dims[Z];
		int nx = dims[X]; int ny = dims[Y]; int nz = dims[Z]; 
		float rx = res[X]; float ry = res[Y]; float rz = res[Z]; 
		
		// at least two layer surfaces (gwb and cgb)
		float[] layers = Interface.getFloatImage4D(layersImage);
		
		// create a mask for all the regions outside of the area where layer 1 is > 0 and layer N is < 0
		boolean[] ctxmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			ctxmask[xyz] = (layers[xyz]>=0.0 && layers[xyz+nlayers*nxyz]<=0.0);
		}
		
		if (Interface.isValid(maskImage)) {
			byte[] mask = Interface.getUByteImage3D(maskImage);
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]==0) ctxmask[xyz] = false;
		}
		
		// find if we risk reaching an image boundary
		boolean checkBoundaries = false;
		int bdist = 3;
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				if (x<bdist || x>=nx-bdist || y<bdist || y>=bdist || z<bdist || z>=bdist) checkBoundaries = true;
			}
		}
		if (checkBoundaries) System.out.println("systematically check for image boundaries");
		
		// main algorithm
		float delta = paramExtent.getValue().floatValue()/Numerics.min(rx,ry,rz);
		float thickness = 10.0f/Numerics.min(rx,ry,rz);
		
		System.out.println("Define parcels (size: "+delta+")\n");
		/*
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		int[] parcels = new int[nx*ny*nz];
		int nparcel = 0;
		BitSet sampled = new BitSet(nx*ny*nz);
		
		// make an array of available samples
		BitSet available = new BitSet(nx*ny*nz);
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) available.set(xyz);
		}
		
		BinaryHeap2D heap = new BinaryHeap2D(Numerics.ceil(2*delta*thickness), BinaryHeap2D.MINTREE);
		
		// for random sampling
		Well19937c random = new Well19937c();
		while (available.size()>0) {
			System.out.print(".");
			int repeat = random.nextInt(available.size());
			int idx=0;
			for (int n=0;n<repeat;n++) {
				idx = available.nextSetBit(idx);
			}
			int zs = Numerics.floor(idx/(nx*ny));
			int ys = Numerics.floor((idx-nx*ny*zs)/nx);
			int xs = idx-nx*ys-nx*ny*zs;
			
			// new point
			nparcel++;
			growParcel(parcels, nparcel, xs, ys, zs, delta, layers, available, profile, sampled, heap, nx,ny,nz,nlayers);
		}*/
		
		int[] parcels = new int[nx*ny*nz];
		float[] distance = new float[nx*ny*nz];
		//growAllParcels(parcels, distance, layers, ctxmask, delta, nx, ny, nz, nlayers, rx, ry, rz);
				
		// 1. estimate the number of needed parcels from voxel coverage of middle layer(s) in cubes of +/- delta in size
		int zero = Numerics.ceil(delta);
		int step = 2*zero;
		System.out.println("spacing: "+step);
		System.out.println("max labels: "+nxyz/(step*step*step));
		
		BitSet sampled = new BitSet(nx*ny*nz);
		for (int x=zero; x<nx-zero; x+=step) for (int y=zero; y<ny-zero; y+=step) for (int z=zero; z<nz-zero; z+=step) {
			int xyz = x + nx*y + nx*ny*z;
			// find closest point within window
			if (ctxmask[xyz]) sampled.set(xyz);
			else {
				int dist = 1;
				while (dist<zero) {
					for (byte dir=0;dir<26;dir++) {
						int ngb = xyz + dist*Ngb.x[dir] + dist*Ngb.y[dir]*nx + dist*Ngb.z[dir]*nx*ny;
						if (ctxmask[ngb]) {
							sampled.set(ngb);
							dir = 26;
							dist = zero+1;
						}
					}
					dist++;
				}
			}
		}
		int nparcels = sampled.cardinality();
		System.out.println("Number of parcels: "+nparcels);
		
		// 2. initalize the parcels and their corresponding profiles into a single heap
		//int[] parcels = new int[nx*ny*nz];
		//float[] distance = new float[nx*ny*nz];
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		int nctx = 0;
		for (int xyz=0;xyz<nx*ny*nz;xyz++) if (ctxmask[xyz]) nctx++;
		System.out.println("cortex mask size: "+nctx+" voxels");
		System.out.println("layers: "+nlayers);
		System.out.flush();
		
		BinaryHeapPair heap = new BinaryHeapPair(nctx, BinaryHeapPair.MINTREE);
		
		int idx=0;
		for (int n=0;n<sampled.cardinality();n++) {
			// get the next non-zero value
			idx = sampled.nextSetBit(idx+1);
		
			/*
			int zs = Numerics.floor(idx/(nx*ny));
			int ys = Numerics.floor((idx-nx*ny*zs)/nx);
			int xs = idx-nx*ys-nx*ny*zs;
			
			// compute profile centered on point of interest
			profile.computeTrajectory(layers, xs, ys, zs);
		
			// add the underlying voxels to the labelled region
			for (int l=0;l<=nlayers;l++) {
				int xyzs = Numerics.round(profile.getPt(l)[X]) + nx*Numerics.round(profile.getPt(l)[Y]) + nx*ny*Numerics.round(profile.getPt(l)[Z]);
				if (ctxmask[xyzs]) heap.addValue(0.0f,xyzs,(n+1));
			}
			*/
			heap.addValue(0.0f,idx,(n+1));
		}
		
		// grow all the heap points until fully sampled
		sampled.clear();
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			float curdist = heap.getFirst();
			int xyz = heap.getFirstId1();
			int lb = heap.getFirstId2();
			heap.removeFirst();
			
			// for points selected multiple times
			if (!sampled.get(xyz)) {
				//System.out.print(".");
				// add to sampled values
				parcels[xyz] = lb;
				distance[xyz] = curdist;
				sampled.set(xyz);
				/*
				for (byte d=0;d<26;d++) {
					int ngb = Ngb.neighborIndex(d,xyz,nx,ny,nz);
					int zn = Numerics.floor(ngb/(nx*ny));
					int yn = Numerics.floor((ngb-nx*ny*zn)/nx);
					int xn = ngb-nx*yn-nx*ny*zn;
					if (xn>=0 && xn<nx && yn>=0 && yn<ny && zn>=0 && zn<nz) {
						if (!sampled.get(ngb) && ctxmask[ngb]) {
							// add to the tree no matter what (guarantees to fill the brain with labels)
							float newdist = curdist+Ngb.neighborDistance(d);
							//if (newdist<delta) 
							heap.addValue(newdist, ngb, lb);
						}
					}
				}*/
				// add full profile to sampled values
				int zs = Numerics.floor(xyz/(nx*ny));
				int ys = Numerics.floor((xyz-nx*ny*zs)/nx);
				int xs = xyz-nx*ys-nx*ny*zs;
			
				// compute profile centered on point of interest
				profile.computeTrajectory(layers, xs, ys, zs);
		
				// add the underlying voxels to the labelled region
				for (int l=0;l<=nlayers;l++) {
					int xyzs = Numerics.round(profile.getPt(l)[X]) + nx*Numerics.round(profile.getPt(l)[Y]) + nx*ny*Numerics.round(profile.getPt(l)[Z]);
					if (checkBoundaries) {
						int zn = Numerics.floor(xyzs/(nx*ny));
						int yn = Numerics.floor((xyzs-nx*ny*zn)/nx);
						int xn = xyzs-nx*yn-nx*ny*zn;
						if (xn>=0 && xn<nx && yn>=0 && yn<ny && zn>=0 && zn<nz) {
							if (ctxmask[xyzs] && !sampled.get(xyzs)) {
								parcels[xyzs] = lb;
								distance[xyzs] = curdist;
								sampled.set(xyzs);
							}
						}
					} else {
						if (ctxmask[xyzs] && !sampled.get(xyzs)) {
							parcels[xyzs] = lb;
							distance[xyzs] = curdist;
							sampled.set(xyzs);
						}
					}
				}
				
				// add neighbors, prefer ones further from boundaries
				for (byte d=0;d<26;d++) {
					int ngb = Ngb.neighborIndex(d,xyz,nx,ny,nz);
					if (checkBoundaries) {
						int zn = Numerics.floor(ngb/(nx*ny));
						int yn = Numerics.floor((ngb-nx*ny*zn)/nx);
						int xn = ngb-nx*yn-nx*ny*zn;
						if (xn>=0 && xn<nx && yn>=0 && yn<ny && zn>=0 && zn<nz) {
							if (!sampled.get(ngb) && ctxmask[ngb]) {
								// add to the tree no matter what (guarantees to fill the brain with labels)
								float newdist = curdist+Ngb.neighborDistance(d)+Numerics.abs(layers[ngb+nlayers*nxyz]-layers[ngb])/Numerics.max(0.1f,layers[ngb+nlayers*nxyz]-layers[ngb]);
								//float newdist = curdist+Ngb.neighborDistance(d);
								//if (newdist<delta) 
								heap.addValue(newdist, ngb, lb);
							}
						}
					} else {
						if (!sampled.get(ngb) && ctxmask[ngb]) {
							// add to the tree no matter what (guarantees to fill the brain with labels)
							float newdist = curdist+Ngb.neighborDistance(d)+Numerics.abs(layers[ngb+nlayers*nxyz]-layers[ngb])/Numerics.max(0.1f,layers[ngb+nlayers*nxyz]-layers[ngb]);
							//float newdist = curdist+Ngb.neighborDistance(d);
							//if (newdist<delta) 
							heap.addValue(newdist, ngb, lb);
						}
					}
				}

			}
		}

		System.out.println("\n Output");
		
		// output
		Interface.setIntegerImage3D(parcels, dims, mappedImage, name+"_parcels", header);
		Interface.setFloatImage3D(distance, dims, distanceImage, name+"_dist", header);
		
	}

	/*
	void growParcel(int[] parcels, int parcelid, int x, int y, int z, float delta, float[][] layers, BitSet available, CorticalProfile profile, BitSet sampled, BinaryHeap2D heap, int nx, int ny, int nz, int nlayers) {
		// compute profile centered on point of interest
		profile.computeTrajectory(layers, x, y, z);
		
		// add the underlying voxels to the estimation matrix
		sampled.clear();
		for (int l=0;l<=nlayers;l++) {
			int xyzs = Numerics.round(profile.getPt(l)[X]) + nx*Numerics.round(profile.getPt(l)[Y]) + nx*ny*Numerics.round(profile.getPt(l)[Z]);
			if (available.get(xyzs)) sampled.set(xyzs);
		}
		
		// find enighbors through fast marching
		heap.reset();
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
			if (!sampled.get(xyz)) {
				// add to sampled values
				sampled.set(xyz);
				
				for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) {
					int ngb = xyz+dx+nx*dy+nx*ny*dz;
					if (!sampled.get(ngb) && available.get(ngb)) {
						// add to the tree if close enough from first profile
						float newdist = curdist+(float)FastMath.sqrt(dx*dx+dy*dy+dz*dz);
						if (newdist<=delta) heap.addValue(newdist, ngb, (byte)1);
					}
				}
			}
		}
		idx = 0;
		for (int n=0;n<sampled.size();n++) {
			// get the next non-zero value
			idx = sampled.nextSetBit(idx);
			parcels[idx] = parcelid;
			available.set(idx, false);
		}
	}

	// better approach: grow all parcels at once?
	private final void growAllParcels(int[] parcels, float[] distance, float[] layers, boolean[] ctxmask, float delta, int nx, int ny, int nz, int nlayers, float rx, float ry, float rz) {
		// 1. estimate the number of needed parcels from voxel coverage of middle layer(s) in cubes of +/- delta in size
		int zero = Numerics.ceil(delta);
		int step = 2*zero;
		BitSet sampled = new BitSet(nx*ny*nz);
		for (int x=zero; x<nx-zero; x+=step) for (int y=zero; y<ny-zero; y+=step) for (int z=zero; z<nz-zero; z+=step) {
			int xyz = x + nx*y + nx*ny*z;
			// find closest point within window
			if (ctxmask[xyz]) sampled.set(xyz);
			else {
				int dist = 1;
				while (dist<zero-1) {
					for (byte dir=0;dir<26;dir++) {
						int ngb = xyz + dist*Ngb.neighborIndex(dir, xyz, nx, ny, nz);
						if (ctxmask[ngb]) {
							sampled.set(ngb);
							dir = 26;
							dist = Numerics.ceil(delta);
						}
					}
					dist++;
				}
			}
		}
		int nparcels = sampled.size();
		System.out.println("Number of parcels: "+nparcels);
		
		// 2. initalize the parcels and their corresponding profiles into a single heap
		//int[] parcels = new int[nx*ny*nz];
		//float[] distance = new float[nx*ny*nz];
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		int nctx = 0;
		for (int xyz=0;xyz<nx*ny*nz;xyz++) if (ctxmask[xyz]) nctx++;
		BinaryHeapPair heap = new BinaryHeapPair(nctx, BinaryHeapPair.MINTREE);
		
		int idx=0;
		for (int n=0;n<sampled.size();n++) {
			// get the next non-zero value
			idx = sampled.nextSetBit(idx);
		
			int zs = Numerics.floor(idx/(nx*ny));
			int ys = Numerics.floor((idx-nx*ny*zs)/nx);
			int xs = idx-nx*ys-nx*ny*zs;
			
			// compute profile centered on point of interest
			profile.computeTrajectory(layers, xs, ys, zs);
		
			// add the underlying voxels to the labelled region
			for (int l=0;l<=nlayers;l++) {
				int xyzs = Numerics.round(profile.getPt(l)[X]) + nx*Numerics.round(profile.getPt(l)[Y]) + nx*ny*Numerics.round(profile.getPt(l)[Z]);
				heap.addValue(0.0f,xyzs,(n+1));
			}
		
		}
		
		// grow all the heap points until fully sampled
		sampled.clear();
		while (heap.isNotEmpty()) {
			// extract point with minimum distance
			float curdist = heap.getFirst();
			int xyz = heap.getFirstId1();
			int lb = heap.getFirstId2();
			heap.removeFirst();
			
			// for points selected multiple times
			if (!sampled.get(xyz)) {
				// add to sampled values
				parcels[xyz] = lb;
				distance[xyz] = curdist;
				sampled.set(xyz);
				
				for (byte d=0;d<26;d++) {
					int ngb = Ngb.neighborIndex(d,xyz,nx,ny,nz);
					if (!sampled.get(ngb)) {
						// add to the tree no matter what (guarantees to fill the brain with labels)
						float newdist = curdist+Ngb.neighborDistance(d);
						heap.addValue(newdist, ngb, lb);
					}
				}
			}
		}
		return;
	}*/
	
}
