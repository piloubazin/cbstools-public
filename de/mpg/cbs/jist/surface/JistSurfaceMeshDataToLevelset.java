package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;
import org.apache.commons.math3.util.FastMath;

import edu.jhu.ece.iacl.algorithms.graphics.intersector.SurfaceIntersector;
import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.ResampleLevelSet;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.SurfaceToMask;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.topology.TopologyCorrection;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.*;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.cs.cisst.algorithms.geometry.surface.*;
import edu.jhu.cs.cisst.algorithms.segmentation.gac.DistanceField3D;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;


public class JistSurfaceMeshDataToLevelset  extends ProcessingAlgorithm{
	ParamVolume originalImage;
	ParamSurface originalSurface;
	ParamVolume dataImage;
	ParamVolume maskImage;
	ParamBoolean mipavTransform;
	ParamFloat dataExtension;
	ParamInteger dataIndex;
	
	private static final float	UNKNOWN = -1.0f;
       
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(originalSurface=new ParamSurface("Surface Mesh"));
		inputParams.add(originalImage=new ParamVolume("Reference Volume"));
		originalImage.setDescription("Reference volume to use for level set representation dimensions.");
		//inputParams.add(mipavTransform=new ParamBoolean("Align to MIPAV image space", true));
		inputParams.add(dataExtension=new ParamFloat("Data Extension Distance", 0.0f, 100.0f, 5.0f));
		inputParams.add(dataIndex=new ParamInteger("Data Index", 0, 10, 0));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setName("MeshDataToVolume");
		inputParams.setLabel("Mesh Data to Volume");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Embeds data from an input surface mesh (e.g. curvature, mapped thickness or contrast) into volumetric space.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(dataImage=new ParamVolume("Surface Data Volume"));
		outputParams.add(maskImage=new ParamVolume("Data Mask Volume"));
	}
	
	protected void execute(CalculationMonitor monitor) {
		
		ImageDataFloat	originalImgData = new ImageDataFloat(originalImage.getImageData());
		int nx = originalImgData.getRows();
		int ny = originalImgData.getCols();
		int nz = originalImgData.getSlices();
		float rx = originalImgData.getHeader().getDimResolutions()[0];
		float ry = originalImgData.getHeader().getDimResolutions()[1];
		float rz = originalImgData.getHeader().getDimResolutions()[2];
		//float ox = originalImgData.getHeader().getOrigin()[0];
		//float oy = originalImgData.getHeader().getOrigin()[1];
		//float oz = originalImgData.getHeader().getOrigin()[2];
		
		EmbeddedSurface surf = originalSurface.getSurface();
		
		/* avoid messy transformations
		// maps the surface to the voxel space
		if (mipavTransform.getValue().booleanValue()) {
			surf.scaleVertices(new float[]{-1.0f, 1.0f, -1.0f});
			Point3f pt0 = new Point3f((nx-1)*rx/2.0f, (ny-1)*ry/2.0f, (nz-1)*rz/2.0f);
			surf.translate(pt0);
		}
		*/
		EmbeddedSurface rsurf = surf.clone();
		
		// scale the surface data into voxel space
		rsurf.scaleVertices(new float[]{1.0f/rx, -1.0f/ry, -1.0f/rz});
		Point3f pt1 = new Point3f(0, (ny-1), (nz-1));
		rsurf.translate(pt1);
		
		
		float[][][] data = new float[nx][ny][nz];
		float[][][] dist = new float[nx][ny][nz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) {
			dist[x][y][z] = UNKNOWN;
		}
		byte[][][] mask = new byte[nx][ny][nz];
		
		int idx = dataIndex.getValue().intValue();
		for(int i=0; i<rsurf.getVertexData().length; i++){
			Point3f p = rsurf.getVertex(i);
			
			// mapping from mesh to voxel space : how to interpolate??
			/*
			int x = Numerics.round(p.x/rx);
			int y = Numerics.round((ny-1)-p.y/ry);
			int z = Numerics.round((nz-1)-p.z/rz);
			*/
			// map the 8-closest neighbors rather
			int x = Numerics.floor(p.x);
			int y = Numerics.floor(p.y);
			int z = Numerics.floor(p.z);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				float dX = (float)FastMath.sqrt( (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z) );
				if (dist[x][y][z]==UNKNOWN) {
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				} else if (dX<dist[x][y][z]) {
					// by default use the closest surface point to a given voxel
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				}
			}
			x = Numerics.ceil(p.x);
			y = Numerics.floor(p.y);
			z = Numerics.floor(p.z);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				float dX = (float)FastMath.sqrt( (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z) );
				if (dist[x][y][z]==UNKNOWN) {
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				} else if (dX<dist[x][y][z]) {
					// by default use the closest surface point to a given voxel
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				}
			}
			x = Numerics.floor(p.x);
			y = Numerics.ceil(p.y);
			z = Numerics.floor(p.z);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				float dX = (float)FastMath.sqrt( (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z) );
				if (dist[x][y][z]==UNKNOWN) {
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				} else if (dX<dist[x][y][z]) {
					// by default use the closest surface point to a given voxel
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				}
			}
			x = Numerics.floor(p.x);
			y = Numerics.floor(p.y);
			z = Numerics.ceil(p.z);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				float dX = (float)FastMath.sqrt( (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z) );
				if (dist[x][y][z]==UNKNOWN) {
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				} else if (dX<dist[x][y][z]) {
					// by default use the closest surface point to a given voxel
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				}
			}
			x = Numerics.ceil(p.x);
			y = Numerics.ceil(p.y);
			z = Numerics.floor(p.z);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				float dX = (float)FastMath.sqrt( (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z) );
				if (dist[x][y][z]==UNKNOWN) {
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				} else if (dX<dist[x][y][z]) {
					// by default use the closest surface point to a given voxel
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				}
			}
			x = Numerics.floor(p.x);
			y = Numerics.ceil(p.y);
			z = Numerics.ceil(p.z);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				float dX = (float)FastMath.sqrt( (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z) );
				if (dist[x][y][z]==UNKNOWN) {
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				} else if (dX<dist[x][y][z]) {
					// by default use the closest surface point to a given voxel
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				}
			}
			x = Numerics.ceil(p.x);
			y = Numerics.floor(p.y);
			z = Numerics.ceil(p.z);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				float dX = (float)FastMath.sqrt( (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z) );
				if (dist[x][y][z]==UNKNOWN) {
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				} else if (dX<dist[x][y][z]) {
					// by default use the closest surface point to a given voxel
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				}
			}
			x = Numerics.ceil(p.x);
			y = Numerics.ceil(p.y);
			z = Numerics.ceil(p.z);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				float dX = (float)FastMath.sqrt( (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z) );
				if (dist[x][y][z]==UNKNOWN) {
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				} else if (dX<dist[x][y][z]) {
					// by default use the closest surface point to a given voxel
					data[x][y][z] = (float)(rsurf.getVertexData(i)[idx]);
					dist[x][y][z] = dX;
					mask[x][y][z] = 1;
				}
			}
		}
		
		BasicInfo.displayMessage("propagate surface labels\n");		

		// expand labels out for a certain distance
		if (dataExtension.getValue().floatValue()>0)
			fastMarchingPropagation(data, dist, dataExtension.getValue().floatValue(), nx, ny, nz);
		
		BasicInfo.displayMessage("generate outputs\n");		

		ImageDataFloat dataImg = new ImageDataFloat(data);		
		dataImg.setHeader(originalImage.getImageData().getHeader());
		dataImg.setName(rsurf.getName()+"_data");
		dataImage.setValue(dataImg);
		dataImg = null;
		data = null;
		
		ImageDataUByte maskImg = new ImageDataUByte(mask);		
		maskImg.setHeader(originalImage.getImageData().getHeader());
		maskImg.setName(rsurf.getName()+"_mask");
		maskImage.setValue(maskImg);
		maskImg = null;
		mask = null;
		
	}
	public final void fastMarchingPropagation(float[][][] data, float[][][] distance, float maxDist, int nx, int ny, int nz) {
		 
		// computation variables
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		float curdist,newdist;	
		boolean done, isprocessed;
		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
		BasicInfo.displayMessage("fast marching\n");		

		int[] xoff = new int[]{1, -1, 0, 0, 0, 0};
		int[] yoff = new int[]{0, 0, 1, -1, 0, 0};
		int[] zoff = new int[]{0, 0, 0, 0, 1, -1};

		BinaryHeapPair heap = new BinaryHeapPair(nx+ny+nz, BinaryHeapPair.MINTREE);
		
		BasicInfo.displayMessage("init\n");		
        // initialize the heap from boundaries
        boolean[][][] computed = new boolean[nx][ny][nz];
        for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) if (distance[x][y][z]>UNKNOWN) {
        	// just put in the tree?
        	heap.addValue(distance[x][y][z],x+nx*y+nx*ny*z,x+nx*y+nx*ny*z);
        	// reset distance for computations? no, it brings inaccuracies
        	//distance[x][y][z] = UNKNOWN;
        	/*
        	// search for boundaries
        	for (int n = 0; n<6; n++) {
				int xn = x+xoff[n];
				int yn = y+yoff[n];
				int zn = z+zoff[n];
				if (distance[xn][yn][zn]==UNKNOWN) {
					// add to the heap with previous value
					heap.addValue(1.0f,xn,yn,zn,(int)data[x][y][z]);
				}
            }
            */
        }
		BasicInfo.displayMessage("main loop\n");		

        // grow the labels and functions
        while (heap.isNotEmpty()) {
        	//BasicInfo.displayMessage(".");
        	
        	// extract point with minimum distance
        	curdist = heap.getFirst();
        	int xyz = heap.getFirstId1();
        	int xyz0 = heap.getFirstId2();
        	heap.removeFirst();
        	
        	int z = Numerics.floor(xyz/(nx*ny));
        	int y = Numerics.floor((xyz-nx*ny*z)/nx);
        	int x = (xyz-nx*ny*z-nx*y);
        	
        	int z0 = Numerics.floor(xyz0/(nx*ny));
        	int y0 = Numerics.floor((xyz0-nx*ny*z0)/nx);
        	int x0 = (xyz0-nx*ny*z0-nx*y0);

			// if more than nlb labels have been found already, this is done
			if (!computed[x][y][z]) {
			
				// update the distance functions at the current level
				
				data[x][y][z] = data[x0][y0][z0];
				distance[x][y][z] = curdist;
				computed[x][y][z] = true;
				
				// find new neighbors
				for (int n = 0; n<6; n++) {
					int xn = x+xoff[n];
					int yn = y+yoff[n];
					int zn = z+zoff[n];
					if (xn>1 && xn<nx-2 && yn>1 && yn<ny-2 && zn>1 && zn<nz-2) {
						// must be in outside the object or its processed neighborhood
						if (!computed[xn][yn][zn]) {
							// compute new distance based on processed neighbors for the same object
							for (int l=0; l<6; l++) {
								nbdist[l] = UNKNOWN;
								nbflag[l] = false;
								int xb = xn + xoff[l];
								int yb = yn + yoff[l];
								int zb = zn + zoff[l];
								// note that there is at most one value used here
								if (distance[xb][yb][zb]!=UNKNOWN) { // distances may be set even if not yet computed (on boundaries)
									nbdist[l] = distance[xb][yb][zb];
									nbflag[l] = true;
								}			
							}
							newdist = minimumMarchingDistance(nbdist, nbflag);
						
							if (newdist<=maxDist) {
								// add to the heap
								heap.addValue(newdist,xn+nx*yn+nx*ny*zn,xyz0);
							}
						}
					}
				}
			}
		}

       return;
    }

    
	/**
     * the Fast marching distance computation 
     * (!assumes a 6D array with opposite coordinates stacked one after the other)
     * 
     */
    public final float minimumMarchingDistance(float[] val, boolean[] flag) {
		double s, s2, tmp; 
		int count;
		double dist;
    
        // s = a + b +c; s2 = a*a + b*b +c*c
        s = 0;
        s2 = 0;
        count = 0;

        for (int n=0; n<6; n+=2) {
			if (flag[n] && flag[n+1]) {
				tmp = Numerics.min(val[n], val[n+1]); // Take the smaller one if both are processed
				s += tmp;
				s2 += tmp*tmp;
				count++;
			} else if (flag[n]) {
				s += val[n]; // Else, take the processed one
				s2 += val[n]*val[n];
				count++;
			} else if (flag[n+1]) {
				s += val[n+1];
				s2 += val[n+1]*val[n+1];
				count++;
			}
		}
         // count must be greater than zero since there must be at least one processed pt in the neighbors
        
        tmp = (s+FastMath.sqrt( (s*s-count*(s2-1.0f))))/count;

        // The larger root
        return (float)tmp;
    }

}
