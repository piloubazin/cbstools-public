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
	private ParamVolume intensityImage;
	private ParamVolume maskImage;
	private	ParamFloat	paramExtent;
	private	ParamFloat	paramStdev;
		
	private ParamVolume mappedImage;
	
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
		inputParams.setCategory("Laminar Analysis");
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
		
		outputParams.setName("parcel images");
		outputParams.setLabel("parcel images");
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
		
		float delta = paramExtent.getValue().floatValue()/Numerics.min(rx,ry,rz);
		
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
		
		// for random sampling
		Well19937c random = new Well19937c();
		while (available.size()>0) {
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
			growParcel(sampled, parcels, nparcel, xs, ys, zs, delta, layers, profile, ctxmask, nx,ny,nz,nlayers);
		
			// remove the sampled points from available list
			idx=0;
			for (int n=0;n<sampled.size();n++) {
				idx = sampled.nextSetBit(idx);
				available.set(idx, false);
			}
		}
		
		// output
		int[][][] map = new int[nx][ny][nz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			map[x][y][z] = parcels[xyz];
		}
		String imgname = intensityImage.getImageData().getName();
		
		ImageDataInt mapData = new ImageDataInt(map);		
		mapData.setHeader(layersImage.getImageData().getHeader());
		mapData.setName(imgname+"_parcels");
		mappedImage.setValue(mapData);
		mapData = null;
		map = null;
	}

	void growParcel(BitSet sampled, int[] parcels, int parcelid, int x, int y, int z, float delta, float[][] layers, CorticalProfile profile, boolean[] ctxmask, int nx, int ny, int nz, int nlayers) {
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
			
			for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) {
				int ngb = xyz+dx+nx*dy+nx*ny*dz;
				if (!sampled.get(ngb) && ctxmask[ngb] && parcels[xyz]==0) {
					// add to the tree if close enough from first profile
					float newdist = curdist+(float)FastMath.sqrt(dx*dx+dy*dy+dz*dz);
					if (newdist<=delta) heap.addValue(newdist, ngb, (byte)1);
				}
			}
		}
		idx = 0;
		for (int n=0;n<sampled.size();n++) {
			// get the next non-zero value
			idx = sampled.nextSetBit(idx);
			parcels[idx] = parcelid;
		}
	}

}
