package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
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


/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityTVDenoising extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputImage;
	private ParamOption ngbParam;
	private ParamFloat ratioParam;
	private ParamInteger binParam;
	private ParamBoolean adjustParam;
	private ParamOption histParam;
	private ParamBoolean wrappedParam;
	
	private ParamVolume denoiseImage;
	private ParamVolume edgeImage;
	private ParamVolume parcelImage;
	private ParamVolume histImage;
	
	// parameters
	private		static final String[]	ngbTypes = {"6C","18C","26C"};
	private		String		ngbType = "6C";
	
	private		static final String[]	histTypes = {"Full_histogram", "KIThreshold","ExpnormKI","NonnegKI"};
	private		String		histType = "Full_histogram";
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		
		inputParams.add(binParam = new ParamInteger("Histogram bins", 10, 10000, 100));
		inputParams.add(histParam = new ParamOption("Noise level estimation", histTypes));
		histParam.setValue(histType);
		
		inputParams.add(ngbParam = new ParamOption("Neighborhood connectivity", ngbTypes));
		ngbParam.setValue(ngbType);
		
		inputParams.add(ratioParam = new ParamFloat("Scaling ratio", 0.0f, 1.0f, 0.05f));
		inputParams.add(adjustParam = new ParamBoolean("Two-level denoising", false));
		inputParams.add(wrappedParam = new ParamBoolean("Wrap intensities", false));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity.devel");
		inputParams.setLabel("TV Denoising");
		inputParams.setName("TVDenoising");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Total-variation image denoising.");
		
		info.setVersion("3.0.2");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(denoiseImage = new ParamVolume("Denoised Image",VoxelType.FLOAT));
		outputParams.add(edgeImage = new ParamVolume("Edge Image",VoxelType.FLOAT));
		outputParams.add(parcelImage = new ParamVolume("Parcellated Image",VoxelType.INT));
		outputParams.add(histImage = new ParamVolume("Histogram Image",VoxelType.UBYTE,-1,-1,-1,-1));
		
		outputParams.setName("tv denoise images");
		outputParams.setLabel("tv denoise images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataFloat inImg = new ImageDataFloat(inputImage.getImageData());
		float[][][] image = inImg.toArray3d();
		int nx = inImg.getRows();
		int ny = inImg.getCols();
		int nz = inImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = inImg.getHeader().getDimResolutions()[0];
		float ry = inImg.getHeader().getDimResolutions()[1];
		float rz = inImg.getHeader().getDimResolutions()[2];
		inImg = null;
		
		// basic mask for zero values
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			//if (image[x][y][z]==0) mask[x][y][z] = false;
			//else mask[x][y][z]=true;
			mask[x][y][z]=true;
		}
		
		// rescale image to [0,1] for convenience
		float INF = 1e9f;
		float Imin = INF, Imax = -INF;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
			if (image[x][y][z]>Imax) Imax = image[x][y][z];
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			image[x][y][z] = (image[x][y][z]-Imin)/(Imax-Imin);
		}	
		// main algorithm
		
		// first estimate the noise level
		int ngb = 1; // 6C neighborhood
		if (ngbParam.getValue().equals("18C")) ngb = 2;
		else if (ngbParam.getValue().equals("26C")) ngb = 3;
		
		System.out.println("build local difference histogram");
		
		Histogram hist = new Histogram(binParam.getValue().intValue());
		
		hist.buildFromDifferences(image, mask, ngb, nx, ny, nz);
		
		//System.out.println(hist.printHistogram());
		
		float stdev = 0.0f;
		int threshold = 0;
		if (histParam.getValue().equals("KIThreshold")) {
			threshold = hist.computeKIThresholdExhaustive();
			// get noise model from threshold
			stdev = hist.getPartialMean(0,threshold);
		} else if (histParam.getValue().equals("ExpnormKI")) {
			threshold = hist.computeExpNormalKIThresholdExhaustive();
			// get noise model from threshold
			stdev = hist.getPartialMean(0,threshold);
		} else if (histParam.getValue().equals("NonnegKI")) {
			threshold = hist.computeNonnegKIThresholdExhaustive();
			// get noise model from threshold
			stdev = hist.getPartialMean(0,threshold);
		} else {
			stdev = hist.mean();
		}
		
		System.out.println("estimated noise stdev: "+(stdev*(Imax-Imin))+" (relative: "+stdev+")");
		
		System.out.println("start TV algorithm");
		
		// apply to TV algorithm
		TotalVariation algo = new TotalVariation(image,mask,nx,ny,nz, ratioParam.getValue().floatValue(), 0.125f, 0.00001f, 500);
		
		if (wrappedParam.getValue().booleanValue())
		    algo.solveWrapped();
		else 
            algo.denoiseImage(stdev, adjustParam.getValue().booleanValue());
		
		// output
		
		// return denoised image, histogram and threshold, estimated edges? (locations where gradient>stdev)
		
		byte[][] histo = hist.plotLogHistogram(threshold);
		
		float[][][] tvimg;
		if (wrappedParam.getValue().booleanValue()) tvimg = algo.exportResultWrapped();
		else tvimg = algo.exportResult();

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tvimg[x][y][z] = Imin + tvimg[x][y][z]*(Imax-Imin);
		}	
		
		float[][][] edges = new float[nx][ny][nz];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			edges[x][y][z] = 0.0f;
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if (i*i+j*j+l*l>0 && i*i+j*j+l*l<=ngb && mask[x+i][y+j][z+l]) {
					edges[x][y][z] = Numerics.max(edges[x][y][z],Numerics.abs(tvimg[x][y][z]-tvimg[x+i][y+j][z+l]));
				}
			}
		}
		
		int[][][] parcel = new int[nx][ny][nz];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			parcel[x][y][z] = 1+Numerics.round((tvimg[x][y][z]-Imin)/(Imax-Imin)/stdev);
		}
		//parcel = ObjectLabeling.connected6Object3D(parcel, nx, ny, nz);

		ImageDataFloat bufferData = new ImageDataFloat(tvimg);
		bufferData.setHeader(inputImage.getImageData().getHeader());
		bufferData.setName(inputImage.getImageData().getName()+"_tvdenoise");
		denoiseImage.setValue(bufferData);
		bufferData = null;
		tvimg = null;
		
		bufferData = new ImageDataFloat(edges);
		bufferData.setHeader(inputImage.getImageData().getHeader());
		bufferData.setName(inputImage.getImageData().getName()+"_tvedges");
		edgeImage.setValue(bufferData);
		bufferData = null;
		edges = null;
		
		ImageDataInt intbufferData = new ImageDataInt(parcel);
		intbufferData.setHeader(inputImage.getImageData().getHeader());
		intbufferData.setName(inputImage.getImageData().getName()+"_tvseg");
		parcelImage.setValue(intbufferData);
		intbufferData = null;
		parcel = null;
		
		ImageDataUByte histData = new ImageDataUByte(histo);		
		histData.setName(inputImage.getImageData().getName()+"_tvhist");
		histImage.setValue(histData);
		histData = null;
		histo = null;
	}


}
