package de.mpg.cbs.jist.shape;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamNumberCollection;
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
public class JistShapeLevelsetFullPCA extends ProcessingAlgorithm {

	// jist containers
	private ParamVolumeCollection intensImages;
	private ParamVolumeCollection lvlImages;
	
	private ParamDouble ratioParam;
	private ParamDouble distanceParam;
	private ParamBoolean cropParam;
	private ParamString nameParam;
	private ParamDouble scaleParam;

	private ParamVolume 	smeanImage;
	private ParamVolume 	spcaImage;
	private ParamVolume 	ameanImage;
	private ParamVolume 	apcaImage;
	
	private ParamVolume 	segImage;
	private ParamVolume 	imgImage;
	private	ParamNumberCollection	eigenValues;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(lvlImages = new ParamVolumeCollection("Levelset Image collection"));
		inputParams.add(intensImages = new ParamVolumeCollection("Intensity Image collection"));
		
		inputParams.add(ratioParam = new ParamDouble("Ratio of first to last kept eigenvalue", 0.0f, 1.0f, 0.1f));
		inputParams.add(distanceParam = new ParamDouble("Distance to levelset (mm)", 0.0f, 100.0f, 5.0f));
		inputParams.add(cropParam = new ParamBoolean("Output cropped atlas", false));
		inputParams.add(nameParam = new ParamString("Atlas name"));
		inputParams.add(scaleParam = new ParamDouble("Eigenmode visualization scaling (x stdev)", 0.0f, 10.0f, 1.0f));
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Shape.devel");
		inputParams.setLabel("Full Levelset PCA");
		inputParams.setName("FullLevelsetPCA");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Create an active appearance model (indep shape) from segmentations and a shape average.");
		
		info.setVersion("1.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(smeanImage = new ParamVolume("Average Shape",VoxelType.FLOAT));
		outputParams.add(spcaImage = new ParamVolume("PCA Shape Component",VoxelType.FLOAT, -1, -1, -1, -1));
		outputParams.add(ameanImage = new ParamVolume("Average Intensity",VoxelType.FLOAT));
		outputParams.add(apcaImage = new ParamVolume("PCA Component Intensity",VoxelType.FLOAT, -1, -1, -1, -1));
		outputParams.add(eigenValues = new ParamNumberCollection("Eigenvalues"));
		outputParams.add(segImage = new ParamVolume("Segment Image",VoxelType.UBYTE, -1, -1, -1, -1));
		outputParams.add(imgImage = new ParamVolume("PCA-shifted Image",VoxelType.FLOAT, -1, -1, -1, -1));
		
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
		
		float[][] intensities = new float[nsubj][nx*ny*nz];
		float[][] levelsets = new float[nsubj][nx*ny*nz];
		for (int n=0;n<nsubj;n++) {
			System.out.println("image "+(n+1));
			lvlImg = new ImageDataFloat(lvlImages.getImageDataList().get(n));
			ImageDataFloat intensImg = new ImageDataFloat(intensImages.getImageDataList().get(n));
			float[][][] tmpl = lvlImg.toArray3d();
			float[][][] tmpi = intensImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				intensities[n][xyz] = tmpi[x][y][z];
				levelsets[n][xyz] = tmpl[x][y][z];
			}
		}
		
		// crop the data to a reasonable size: data size + 3 voxels
		System.out.println("compute cropping parameters");
		boolean crop = cropParam.getValue().booleanValue();
		float maxdist = distanceParam.getValue().floatValue()/Numerics.max(rx,ry,rz);
		int x0 = nx, xN = 0;
		int y0 = ny, yN = 0;
		int z0 = nz, zN = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			boolean keep=false;
			for (int n=0;n<nsubj;n++) if (levelsets[n][xyz] < maxdist+3.0f) keep = true;
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

		System.out.println("compute shape average");
		float[][][] saverage = null;
		if (crop) {
			saverage = new float[ncx][ncy][ncz];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				int xyz = x+nx*y+nx*ny*z;
				
				saverage[x-x0][y-y0][z-z0] = 0.0f;
				for (int n=0;n<nsubj;n++) saverage[x-x0][y-y0][z-z0] += levelsets[n][xyz]/(float)nsubj;
			}
		} else {
			saverage = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				
				saverage[x][y][z] = 0.0f;
				for (int n=0;n<nsubj;n++) saverage[x][y][z] += levelsets[n][xyz]/(float)nsubj;
			}
		}
		// compute mean volume to normalize the average to the same size
		double meanvol = 0.0;
		float mindist = 1e6f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			for (int n=0;n<nsubj;n++) if (levelsets[n][xyz]<0) meanvol += 1.0/nsubj;
			if (saverage[x][y][z]<mindist) mindist = saverage[x][y][z];
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
				if (saverage[x][y][z]<0 || saverage[x][y][z]==mindist) {
					vol++;
					label[x][y][z] = true;
				
					// check for neighbors outside
					for (int k = 0; k<6 ; k++) {
						if (saverage[x+xoff[k]][y+yoff[k]][z+zoff[k]]>=0) {
							// add to the heap (multiple instances are OK)
							heap.addValue(saverage[x+xoff[k]][y+yoff[k]][z+zoff[k]],x+xoff[k],y+yoff[k],z+zoff[k],1);
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
							heap.addValue(saverage[x+xoff[k]][y+yoff[k]][z+zoff[k]],x+xoff[k],y+yoff[k],z+zoff[k],1);
						}
					}
				}
			}
			System.out.println("Distance offset: "+threshold);
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				saverage[x][y][z] -= threshold;
			}
		}
		
		// project intensities onto average shape
		System.out.println("create pca matrix");
		double[][] pcamatrix = new double[2*ncx*ncy*ncz][nsubj];

		System.out.println("create shape matrix");
		for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
			int xyz = x+nx*y+nx*ny*z;
			int xyzc = x-x0 + ncx*(y-y0) + ncx*ncy*(z-z0);
			
			for (int n=0;n<nsubj;n++) {
				if (crop) 	pcamatrix[xyzc][n] = levelsets[n][xyz]-saverage[x-x0][y-y0][z-z0];
				else 		pcamatrix[xyzc][n] = levelsets[n][xyz]-saverage[x][y][z];
			}
		}

		System.out.println("project intensities on average shape");
		
		float[] target = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			target[xyz] = saverage[x][y][z];
		}
		//saverage = null;
		
		boolean[] mask = new boolean[nx*ny*nz];
		//float maxls = 25.0f;
		//for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
		for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
			int xyz = x+nx*y+nx*ny*z;
			mask[xyz] = true;
			/*
			for (int n=0;n<nsubj;n++) if (levelsets[n][xyz]>maxdist) mask[xyz] = false;
			if (target[xyz]>maxls) mask[xyz] = false;
			*/
		}
		
		for (int n=0;n<nsubj;n++) {
			SmoothGdm evolve = new SmoothGdm(levelsets[n], target, nx, ny, nz, rx, ry, rz, mask, 0.5f, 0.2f, "no",null);
			//evolve.changeNarrowBandSize(20.0f);
			evolve.addIntensityToMap(intensities[n]);
			//evolve.evolveNarrowBandMapping(500, 0.001f);
			evolve.evolveCompleteMapping(500, 0.001f);
			
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				int xyz = x+nx*y+nx*ny*z;
				int xyzc = x-x0 + ncx*(y-y0) + ncx*ncy*(z-z0) + ncx*ncy*ncz;
				// only keep the intensities within a good distance
				if (target[xyz]<=maxdist) {
					pcamatrix[xyzc][n] = evolve.getIntensity()[xyz];
				}
			}
		}
		intensities = null;
		levelsets = null;
		target = null;
					
		System.out.println("compute intensity average");
		float[][][] amean = null;
		if (crop) {
			amean = new float[ncx][ncy][ncz];
			for (int x=0;x<ncx;x++) for (int y=0;y<ncy;y++) for (int z=0;z<ncz;z++) {
				int xyzc = x + ncx*y + ncx*ncy*z + ncx*ncy*ncz;
				
				amean[x][y][z] = 0.0f;
				for (int n=0;n<nsubj;n++) amean[x-x0][y-y0][z-z0] += pcamatrix[xyzc][n]/(float)nsubj;
			}
		} else {
			amean = new float[nx][ny][nz];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				int xyzc = x-x0 + ncx*(y-y0) + ncx*ncy*(z-z0) + ncx*ncy*ncz;
			
				amean[x][y][z] = 0.0f;
				for (int n=0;n<nsubj;n++) amean[x][y][z] += pcamatrix[xyzc][n]/(float)nsubj;
			}
		}
		
		System.out.println("create shape and appearance matrix");
		for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
			int xyz = x+nx*y+nx*ny*z;
			int xyzc = x-x0 + ncx*(y-y0) + ncx*ncy*(z-z0) + ncx*ncy*ncz;
			
			for (int n=0;n<nsubj;n++) {
				if (crop) 	pcamatrix[xyzc][n] -= amean[x-x0][y-y0][z-z0];
				else 		pcamatrix[xyzc][n] -= amean[x][y][z];
			}
		}
		
		// compute PCA components
		System.out.println("perform SVD");
		Matrix M = new Matrix(pcamatrix);
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
		pcamatrix = null;
		
		// normalize the vectors by eigenvalues and output
		float[][][][] scomponents = null;
		
		if (crop) {
			scomponents = new float[ncx][ncy][ncz][npca];
			for (int x=0;x<ncx;x++) for (int y=0;y<ncy;y++) for (int z=0;z<ncz;z++) {
				int xyzc = x + ncx*y + ncx*ncy*z;
				for (int n=0;n<npca;n++) {
					scomponents[x][y][z][n] = (float)(U.get(xyzc, n));
				}
			}
		} else {
			scomponents = new float[nx][ny][nz][npca];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				int xyzc = x-x0 + ncx*(y-y0) + ncx*ncy*(z-z0);
				for (int n=0;n<npca;n++) {
					scomponents[x][y][z][n] = (float)(U.get(xyzc, n));
				}
			}
		}
		
		float[][][][] acomponents = null;
		
		if (crop) {
			acomponents = new float[ncx][ncy][ncz][npca];
			for (int x=0;x<ncx;x++) for (int y=0;y<ncy;y++) for (int z=0;z<ncz;z++) {
				int xyzc = x + ncx*y + ncx*ncy*z + ncx*ncy*ncz;
				for (int n=0;n<npca;n++) {
					acomponents[x][y][z][n] = (float)(U.get(xyzc, n));
				}
			}
		} else {
			acomponents = new float[nx][ny][nz][npca];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				int xyzc = x-x0 + ncx*(y-y0) + ncx*ncy*(z-z0) + ncx*ncy*ncz;
				for (int n=0;n<npca;n++) {
					acomponents[x][y][z][n] = (float)(U.get(xyzc, n));
				}
			}
		}
		U = null;
		
		byte[][][][] segments = null;
		float scale = scaleParam.getValue().floatValue();
		if (crop) {
			segments = new byte[ncx][ncy][ncz][npca+1];
			for (int x=0;x<ncx;x++) for (int y=0;y<ncy;y++) for (int z=0;z<ncz;z++) {
				if (saverage[x][y][z]<0) segments[x][y][z][0] = 3;
				for (int n=0;n<npca;n++) {
					if (saverage[x][y][z]+scale*eigenval[n]*scomponents[x][y][z][n]<0 
						&& saverage[x][y][z]-scale*eigenval[n]*scomponents[x][y][z][n]<0)
						segments[x][y][z][n+1] = 3;
					else if (saverage[x][y][z]+scale*eigenval[n]*scomponents[x][y][z][n]<0)
						segments[x][y][z][n+1] = 2;
					else if (saverage[x][y][z]-scale*eigenval[n]*scomponents[x][y][z][n]<0)
						segments[x][y][z][n+1] = 1;
					else
						segments[x][y][z][n+1] = 0;
				}
			}
		} else {
			segments = new byte[nx][ny][nz][npca+1];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				if (saverage[x][y][z]<0) segments[x][y][z][0] = 3;
				for (int n=0;n<npca;n++) {
					if (saverage[x][y][z]+scale*eigenval[n]*scomponents[x][y][z][n]<0 
						&& saverage[x][y][z]-scale*eigenval[n]*scomponents[x][y][z][n]<0)
						segments[x][y][z][n+1] = 3;
					else if (saverage[x][y][z]+scale*eigenval[n]*scomponents[x][y][z][n]<0)
						segments[x][y][z][n+1] = 2;
					else if (saverage[x][y][z]-scale*eigenval[n]*scomponents[x][y][z][n]<0)
						segments[x][y][z][n+1] = 1;
					else
						segments[x][y][z][n+1] = 0;
				}
			}
		}

		float[][][][] imgs = null;
		if (crop) {
			imgs = new float[ncx][ncy][ncz][npca+1];
			for (int x=0;x<ncx;x++) for (int y=0;y<ncy;y++) for (int z=0;z<ncz;z++) {
				imgs[x][y][z][0] = amean[x][y][z];
				for (int n=0;n<npca;n++) {
					imgs[x][y][z][n+1] = amean[x][y][z] + scale*eigenval[n]*acomponents[x][y][z][n];
				}
			}
		} else {
			imgs = new float[nx][ny][nz][npca+1];
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				imgs[x][y][z][0] = amean[x][y][z];
				for (int n=0;n<npca;n++) {
					imgs[x][y][z][n+1] = amean[x][y][z] + scale*eigenval[n]*acomponents[x][y][z][n];
				}
			}
		}
		
		
		ImageDataFloat smeanData = new ImageDataFloat(saverage);		
		smeanData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		smeanData.setName(nameParam.getValue()+"_savg");
		smeanImage.setValue(smeanData);
		smeanData = null;
		saverage = null;
		
		ImageDataFloat scompData = new ImageDataFloat(scomponents);		
		scompData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		scompData.setName(nameParam.getValue()+"_spca");
		spcaImage.setValue(scompData);
		scompData = null;
		scomponents = null;
		
		ImageDataFloat ameanData = new ImageDataFloat(amean);		
		ameanData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		ameanData.setName(nameParam.getValue()+"_aavg");
		ameanImage.setValue(ameanData);
		ameanData = null;
		amean = null;
		
		ImageDataFloat acompData = new ImageDataFloat(acomponents);		
		acompData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		acompData.setName(nameParam.getValue()+"_apca");
		apcaImage.setValue(acompData);
		acompData = null;
		acomponents = null;
		
		ImageDataUByte segData = new ImageDataUByte(segments);		
		segData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		segData.setName(nameParam.getValue()+"_seg");
		segImage.setValue(segData);
		segData = null;
		segments = null;
		
		ImageDataFloat imgData = new ImageDataFloat(imgs);		
		imgData.setHeader(lvlImages.getImageDataList().get(0).getHeader());
		imgData.setName(nameParam.getValue()+"_img");
		imgImage.setValue(imgData);
		imgData = null;
		imgs = null;
	}


}
