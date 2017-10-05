package de.mpg.cbs.jist.shape;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamNumberCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
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

import Jama.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistShapeLevelsetPCA extends ProcessingAlgorithm {

	// jist containers
	private ParamVolumeCollection lvlImages;
	
	private ParamDouble ratioParam;
	private ParamDouble distanceParam;
	private ParamBoolean cropParam;
	private ParamString nameParam;
	private ParamDouble scaleParam;
	
	private ParamVolume 			meanImage;
	private ParamVolume 			pcaImage;
	private ParamVolume 			segImage;
	private	ParamNumberCollection	eigenValues;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(lvlImages = new ParamVolumeCollection("Levelset Image collection"));
		
		inputParams.add(ratioParam = new ParamDouble("Ratio of first to last kept eigenvalue", 0.0f, 1.0f, 0.1f));
		inputParams.add(distanceParam = new ParamDouble("Distance to levelset (mm)", 0.0f, 100.0f, 5.0f));
		inputParams.add(cropParam = new ParamBoolean("Output cropped atlas", false));
		inputParams.add(nameParam = new ParamString("Atlas name"));
		inputParams.add(scaleParam = new ParamDouble("Eigenmode visualization scaling (x stdev)", 0.0f, 10.0f, 1.0f));
			
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Shape.devel");
		inputParams.setLabel("Levelset PCA");
		inputParams.setName("LevelsetPCA");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Create an active shape model from segmentations via level set PCA.");
		
		info.setVersion("1.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(meanImage = new ParamVolume("Average Image",VoxelType.FLOAT));
		outputParams.add(pcaImage = new ParamVolume("PCA Component Image",VoxelType.FLOAT, -1, -1, -1, -1));
		outputParams.add(segImage = new ParamVolume("Segment Image",VoxelType.UBYTE, -1, -1, -1, -1));
		outputParams.add(eigenValues = new ParamNumberCollection("Eigenvalues"));
		
		outputParams.setName("pca images");
		outputParams.setLabel("pca images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		System.out.println("get all images");
		ImageDataFloat	lvlImg = new ImageDataFloat(lvlImages.getImageDataList().get(0));
		
		int nx = lvlImg.getRows();
		int ny = lvlImg.getCols();
		int nz = lvlImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = lvlImg.getHeader().getDimResolutions()[0];
		float ry = lvlImg.getHeader().getDimResolutions()[1];
		float rz = lvlImg.getHeader().getDimResolutions()[2];
		
		int nsubj = lvlImages.getImageDataList().size();
		
		float[][] levelsets = new float[nsubj][nx*ny*nz];
		for (int n=0;n<nsubj;n++) {
			System.out.println("image "+(n+1));
			lvlImg = new ImageDataFloat(lvlImages.getImageDataList().get(n));
			float[][][] tmp = lvlImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				levelsets[n][xyz] = tmp[x][y][z];
			}
		}
		
		// crop the data to a reasonable size
		System.out.println("compute cropping parameters");
		boolean crop = cropParam.getValue().booleanValue();
		float maxdist = distanceParam.getValue().floatValue()/Numerics.max(rx,ry,rz);
		int x0 = nx, xN = 0;
		int y0 = ny, yN = 0;
		int z0 = nz, zN = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			boolean keep=false;
			for (int n=0;n<nsubj;n++) if (levelsets[n][xyz] < maxdist) keep = true;
			if (keep) {
				x0 = Numerics.min(x0, x);
				xN = Numerics.max(xN, x);
				y0 = Numerics.min(y0, y);
				yN = Numerics.max(yN, y);
				z0 = Numerics.min(z0, z);
				zN = Numerics.max(zN, z);
			}
		}
		int ncx = xN-x0+1;
		int ncy = yN-y0+1;
		int ncz = zN-z0+1;

		System.out.println("compute average");
		float[][][] average = null;
		if (crop) {
			average = new float[ncx][ncy][ncz];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				int xyz = x+nx*y+nx*ny*z;
				
				average[x-x0][y-y0][z-z0] = 0.0f;
				for (int n=0;n<nsubj;n++) average[x-x0][y-y0][z-z0] += levelsets[n][xyz]/(float)nsubj;
			}
		} else {
			average = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				
				average[x][y][z] = 0.0f;
				for (int n=0;n<nsubj;n++) average[x][y][z] += levelsets[n][xyz]/(float)nsubj;
			}
		}
		// compute mean volume to normalize the average to the same size
		double meanvol = 0.0;
		float mindist = 1e6f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			for (int n=0;n<nsubj;n++) if (levelsets[n][xyz]<0) meanvol += 1.0/nsubj;
			if (average[x][y][z]<mindist) mindist = average[x][y][z];
		}
		System.out.println("mean volume (voxels): "+meanvol);
		System.out.println("minimum avg. distance: "+mindist);
		
		// find appropriate threshold to have correct volume; should use a fast marching approach!
		BinaryHeap4D	heap = new BinaryHeap4D(ncx*ncy+ncy*ncz+ncz*ncx, BinaryHeap4D.MINTREE);

		if (crop) {
		} else {
			// 6-neighborhood: pre-compute the index offsets
			int[] xoff = new int[]{1, -1, 0, 0, 0, 0};
			int[] yoff = new int[]{0, 0, 1, -1, 0, 0};
			int[] zoff = new int[]{0, 0, 0, 0, 1, -1};
			
			double vol = 0.0;
			boolean[][][] label = new boolean[nx][ny][nz];
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				if (average[x][y][z]<0 || average[x][y][z]==mindist) {
					vol++;
					label[x][y][z] = true;
				
					// check for neighbors outside
					for (int k = 0; k<6 ; k++) {
						if (average[x+xoff[k]][y+yoff[k]][z+zoff[k]]>=0) {
							// add to the heap (multiple instances are OK)
							heap.addValue(average[x+xoff[k]][y+yoff[k]][z+zoff[k]],x+xoff[k],y+yoff[k],z+zoff[k],1);
						}
					}
				} else {
					label[x][y][z] = false;
				}
			}
			// run until the volume exceeds the mean volume
			float threshold = 0.0f;
			while (heap.isNotEmpty() && vol<meanvol) {
				threshold = heap.getFirst();
				int x = heap.getFirstX();
 				int y = heap.getFirstY();
 				int z = heap.getFirstZ();
 				heap.removeFirst();	
 				if (label[x][y][z]==false) {
					vol++;
					label[x][y][z] = true;
					// add neighbors
					for (int k = 0; k<6; k++) {
						if (label[x+xoff[k]][y+yoff[k]][z+zoff[k]]==false) {
							heap.addValue(average[x+xoff[k]][y+yoff[k]][z+zoff[k]],x+xoff[k],y+yoff[k],z+zoff[k],1);
						}
					}
				}
			}
			System.out.println("Distance offset: "+threshold);
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				average[x][y][z] -= threshold;
			}
		}
		
		System.out.println("create shape matrix");
		double[][] shape = new double[ncx*ncy*ncz][nsubj];
		for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
			int xyz = x+nx*y+nx*ny*z;
			int xyzc = x-x0 + ncx*(y-y0) + ncx*ncy*(z-z0);
			
			for (int n=0;n<nsubj;n++) {
				if (crop) 	shape[xyzc][n] = levelsets[n][xyz]-average[x-x0][y-y0][z-z0];
				else 		shape[xyzc][n] = levelsets[n][xyz]-average[x][y][z];
			}
		}
		levelsets = null;
		
		// compute PCA components
		System.out.println("perform SVD");
		Matrix M = new Matrix(shape);
		SingularValueDecomposition svd = M.svd();
		
		// select the number of components to keep from the variations
		int npca = nsubj;
		for (int n=1;n<npca;n++) {
			if (Numerics.abs(svd.getSingularValues()[n])<ratioParam.getValue().floatValue()*Numerics.abs(svd.getSingularValues()[0])) npca = n;
		}
		System.out.println("number of components used:"+npca);
		float[] eigenval = new float[npca];
		Matrix U = svd.getU();
		System.out.println("eigenvalues:");
		for (int n=0;n<npca;n++) {
			eigenval[n] = (float)Numerics.abs(svd.getSingularValues()[n]);
			System.out.print(eigenval[n]+", ");
			eigenValues.add(eigenval[n]);
		}
		
		// normalize the vectors by eigenvalues and output
		float[][][][] components = null;
		
		if (crop) {
			components = new float[ncx][ncy][ncz][npca];
			for (int x=0;x<ncx;x++) for (int y=0;y<ncy;y++) for (int z=0;z<ncz;z++) {
				int xyzc = x + ncx*y + ncx*ncy*z;
				for (int n=0;n<npca;n++) {
					components[x][y][z][n] = (float)(U.get(xyzc, n));
				}
			}
		} else {
			components = new float[nx][ny][nz][npca];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				int xyzc = x-x0 + ncx*(y-y0) + ncx*ncy*(z-z0);
				for (int n=0;n<npca;n++) {
					components[x][y][z][n] = (float)(U.get(xyzc, n));
				}
			}
		}
		
		byte[][][][] segments = null;
		float scale = scaleParam.getValue().floatValue();
		if (crop) {
			segments = new byte[ncx][ncy][ncz][npca+1];
			for (int x=0;x<ncx;x++) for (int y=0;y<ncy;y++) for (int z=0;z<ncz;z++) {
				if (average[x][y][z]<0) segments[x][y][z][0] = 3;
				for (int n=0;n<npca;n++) {
					if (average[x][y][z]+scale*eigenval[n]*components[x][y][z][n]<0 
						&& average[x][y][z]-scale*eigenval[n]*components[x][y][z][n]<0)
						segments[x][y][z][n+1] = 3;
					else if (average[x][y][z]+scale*eigenval[n]*components[x][y][z][n]<0)
						segments[x][y][z][n+1] = 2;
					else if (average[x][y][z]-scale*eigenval[n]*components[x][y][z][n]<0)
						segments[x][y][z][n+1] = 1;
					else
						segments[x][y][z][n+1] = 0;
				}
			}
		} else {
			segments = new byte[nx][ny][nz][npca+1];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				if (average[x][y][z]<0) segments[x][y][z][0] = 3;
				for (int n=0;n<npca;n++) {
					if (average[x][y][z]+scale*eigenval[n]*components[x][y][z][n]<0 
						&& average[x][y][z]-scale*eigenval[n]*components[x][y][z][n]<0)
						segments[x][y][z][n+1] = 3;
					else if (average[x][y][z]+scale*eigenval[n]*components[x][y][z][n]<0)
						segments[x][y][z][n+1] = 2;
					else if (average[x][y][z]-scale*eigenval[n]*components[x][y][z][n]<0)
						segments[x][y][z][n+1] = 1;
					else
						segments[x][y][z][n+1] = 0;
				}
			}
		}
	
		ImageDataFloat meanData = new ImageDataFloat(average);		
		meanData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		meanData.setName(nameParam.getValue()+"_avg");
		meanImage.setValue(meanData);
		meanData = null;
		average = null;
		
		ImageDataFloat compData = new ImageDataFloat(components);		
		compData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		compData.setName(nameParam.getValue()+"_pca");
		pcaImage.setValue(compData);
		compData = null;
		components = null;
		
		ImageDataUByte segData = new ImageDataUByte(segments);		
		segData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		segData.setName(nameParam.getValue()+"_seg");
		segImage.setValue(segData);
		segData = null;
		segments = null;
	}


}
