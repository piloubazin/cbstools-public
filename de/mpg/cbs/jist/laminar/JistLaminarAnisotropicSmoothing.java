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
public class JistLaminarAnisotropicSmoothing extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume layersImage;
	private ParamVolume intensityImage;
	private ParamVolume maskImage;
	private	ParamFloat	paramLaminarStdev;
	private ParamFloat 	paramColumnarStdev;
	private	ParamFloat	paramLaminarExtent;
	private	ParamInteger	paramVoxelCount;
	
	private ParamOption	paramSmoothingType;
	private final String[] smoothingTypes = {"laminar","columnar","both"};
		
	private ParamOption	paramSizeType;
	private final String[] sizeTypes = {"spatial_extent","voxel_number"};
		
	private ParamVolume smoothedImage;
	
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
		
		inputParams.add(paramSmoothingType = new ParamOption("Smoothing direction",smoothingTypes));
		inputParams.add(paramLaminarStdev = new ParamFloat("Laminar weighting (mm)", 0.0f, 100.0f, 1.0f));
		paramLaminarStdev.setDescription("standard deviation of distance Gaussian prior along the depth lamination");
		inputParams.add(paramColumnarStdev = new ParamFloat("Columnar weighting (mm)", 0.0f, 100.0f, 1.0f));
		paramColumnarStdev.setDescription("standard deviation of distance Gaussian prior along the columnar profiles");
		
		inputParams.add(paramSizeType = new ParamOption("Region size",sizeTypes));
		inputParams.add(paramLaminarExtent = new ParamFloat("Spatial extent (mm)", 0.0f, 100.0f, 3.0f));
		paramLaminarExtent.setDescription("maximum distance to sample if smoothing with a fixed size (ideally, 3x spatial weighting)");
		inputParams.add(paramVoxelCount = new ParamInteger("Number of voxels", 0, 1000, 27));
		paramVoxelCount.setDescription("number of voxels to include if smoothing with constant number of samples");
		
		//inputParams.add(imageParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis.devel");
		inputParams.setLabel("Anisotropic Smoothing");
		inputParams.setName("AnisotropicSmoothing");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Smooth cortical data with geodesic region growing in laminar and/or columnar directions.");
		
		info.setVersion("3.0.9");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(smoothedImage = new ParamVolume("Smoothed Image",null,-1,-1,-1,-1));
		
		outputParams.setName("smoothing images");
		outputParams.setLabel("smoothing images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays : TO DO
		
		ImageDataFloat	layersImg = new ImageDataFloat(layersImage.getImageData());
		
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
		
		ImageDataFloat	intensImg = new ImageDataFloat(intensityImage.getImageData());
		int nt = intensImg.getComponents();
		float[] intensity = new float[nxyz*nt];
		if (nt>1) {
			buffer4 = intensImg.toArray4d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++) {
				int xyz = x+nx*y+nx*ny*z;
				intensity[xyz+nx*ny*nz*t] = buffer4[x][y][z][t];
			}
			buffer4 = null;
		} else {
			float[][][] buffer3 = intensImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				intensity[xyz] = buffer3[x][y][z];
			}
			buffer3 = null;
		}
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
		
		float columnar = paramColumnarStdev.getValue().floatValue()/Numerics.min(rx,ry,rz);
		float laminar = paramLaminarStdev.getValue().floatValue()/Numerics.min(rx,ry,rz);
		
		float extent = paramLaminarExtent.getValue().floatValue()/Numerics.min(rx,ry,rz);
		int voxels = paramVoxelCount.getValue().intValue();
		if (paramSizeType.equals("voxel_number")) extent = -1.0f;
		else voxels = -1;
		
		
		float[] distances = new float[nx*ny*nz];
		BitSet sampled = new BitSet(nx*ny*nz);
		float[] smoothed = new float[nx*ny*nz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				growAnisotropicRegion(sampled, distances, x, y, z, laminar, columnar, extent, voxels, layers, ctxmask, nx,ny,nz,nlayers);
		
				double[] num = new double[nt];
				double den = 0.0;
				int idx=0;
				for (int n=0;n<sampled.size();n++) {
					idx = sampled.nextSetBit(idx);
					double weight = FastMath.exp(-0.5*(distances[idx]*distances[idx]));
					for (int t = 0; t<nt; t++) num[t] += weight*intensity[xyz+nx*ny*nz*t];
					den += weight;
				}
				for (int t = 0; t<nt; t++) smoothed[xyz] = (float)(num[t]/den);
			}
		}
		
		// output
		float[][][] result = new float[nx][ny][nz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			result[x][y][z] = smoothed[xyz];
		}
		String imgname = intensityImage.getImageData().getName();
		
		ImageDataFloat resData = new ImageDataFloat(result);		
		resData.setHeader(layersImage.getImageData().getHeader());
		resData.setName(imgname+"_smooth");
		smoothedImage.setValue(resData);
		resData = null;
		result = null;
	}

	void growAnisotropicRegion(BitSet sampled, float[] distances, int x, int y, int z, float laminar, float columnar, float maxdist, int maxcount, float[][] layers, boolean[] ctxmask, int nx, int ny, int nz, int nlayers) {
		boolean useExtent = (maxdist>0);
		boolean useSize = (maxcount>0);
		
		maxdist = columnar*laminar*maxdist;
		
		int xyz = x+nx*y+nx*ny*z;
		
		// find neighbors through fast marching
		BinaryHeap3D heap = new BinaryHeap3D(100, BinaryHeap3D.MINTREE);
		
		// get normal vector from the combination of closest surfaces
		int minl=-1;
		for (int l=1;l<=nlayers && minl==-1;l++) {
			if (layers[l][xyz]>0) minl = l-1;
		}
		double weight = Numerics.bounded(0.5*(1.0-(layers[minl+1][xyz]+layers[minl][xyz])/(layers[minl+1][xyz]-layers[minl][xyz])),0.0,1.0);
		
		double Dx = weight*(layers[minl+1][xyz+1]-layers[minl+1][xyz-1]) + (1.0-weight)*(layers[minl][xyz+1]-layers[minl][xyz-1]);
		double Dy = weight*(layers[minl+1][xyz+nx]-layers[minl+1][xyz-nx]) + (1.0-weight)*(layers[minl][xyz+nx]-layers[minl][xyz-nx]);
		double Dz = weight*(layers[minl+1][xyz+nx*ny]-layers[minl+1][xyz-nx*ny]) + (1.0-weight)*(layers[minl][xyz+nx*ny]-layers[minl][xyz-nx*ny]);
		
		double grad = FastMath.sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
		
		for (int n=0;n<26;n++) {
			// compute distance
			byte[] dir = Ngb.directionNeighbor(n);
			double dist = (dir[X]*Dx+dir[Y]*Dy+dir[Z]*Dz)/columnar
						+FastMath.sqrt( (dir[X]*grad-dir[X]*Dx)*(dir[X]*grad-dir[X]*Dx)
										+(dir[Y]*grad-dir[Y]*Dy)*(dir[Y]*grad-dir[Y]*Dy)
										+(dir[Z]*grad-dir[Z]*Dz)*(dir[Z]*grad-dir[Z]*Dz) )/laminar;
							
			heap.addValue((float)(dist/grad), x+dir[X],y+dir[Y],z+dir[Z]);
		}
		sampled.clear();
		sampled.set(xyz);
		
		// main loop
		while (heap.isNotEmpty() && (!useSize || sampled.size()<maxcount)) {
			// extract point with minimum distance
			float curdist = heap.getFirst();
			int xi = heap.getFirstX();
			int yi = heap.getFirstY();
			int zi = heap.getFirstZ();
			heap.removeFirst();
			
			int xyzi = xi+nx*yi+nx*ny*zi;
			// for points selected multiple times
			if (sampled.get(xyzi)) continue;
			
			// add to sampled values
			sampled.set(xyzi);
			distances[xyzi] = curdist;
			
			if (!useSize || sampled.size()<maxcount) {
				minl=-1;
				for (int l=1;l<=nlayers && minl==-1;l++) {
					if (layers[l][xyzi]>0) minl = l-1;
				}
				weight = Numerics.bounded(0.5*(1.0-(layers[minl+1][xyzi]+layers[minl][xyzi])/(layers[minl+1][xyzi]-layers[minl][xyzi])),0.0,1.0);
				
				Dx = weight*(layers[minl+1][xyzi+1]-layers[minl+1][xyzi-1]) + (1.0-weight)*(layers[minl][xyzi+1]-layers[minl][xyzi-1]);
				Dy = weight*(layers[minl+1][xyzi+nx]-layers[minl+1][xyzi-nx]) + (1.0-weight)*(layers[minl][xyzi+nx]-layers[minl][xyzi-nx]);
				Dz = weight*(layers[minl+1][xyzi+nx*ny]-layers[minl+1][xyzi-nx*ny]) + (1.0-weight)*(layers[minl][xyzi+nx*ny]-layers[minl][xyzi-nx*ny]);
		
				grad = FastMath.sqrt(Dx*Dx+Dy*Dy+Dz*Dz);

				for (int n=0;n<26;n++) {
					// compute distance
					byte[] dir = Ngb.directionNeighbor(n);				
					int ngb = xyzi+dir[X]+nx*dir[Y]+nx*ny*dir[Z];
					if (!sampled.get(ngb) && ctxmask[ngb]) {
						// compute new distance
						double newdist = (dir[X]*Dx+dir[Y]*Dy+dir[Z]*Dz)/columnar
										+FastMath.sqrt( (dir[X]*grad-dir[X]*Dx)*(dir[X]*grad-dir[X]*Dx)
														+(dir[Y]*grad-dir[Y]*Dy)*(dir[Y]*grad-dir[Y]*Dy)
														+(dir[Z]*grad-dir[Z]*Dz)*(dir[Z]*grad-dir[Z]*Dz) )/laminar;
						
						if (!useExtent || curdist+newdist/grad<maxdist) heap.addValue((float)(curdist+newdist/grad),x+dir[X],y+dir[Y],z+dir[Z]);
					}
				}
			}
		}
	}

}
