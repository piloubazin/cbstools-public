package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
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


/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityMp2rageMasking extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inv2Image;
	private ParamVolume t1mapImage;
	private ParamVolume isoImage;
	private ParamOption distribParam;
	private ParamOption maskingParam;
	
	private ParamVolume probaImage;
	private ParamVolume maskImage;
	private ParamVolume t1maskImage;
	private ParamVolume isomaskImage;
	
	// parameters
	private		static final String[]	distribs = {"exponential","half-normal"};
	private		String		distrib = "half-normal";
	
	private		static final byte	HNORM = 1;
	private		static final byte	EXP = 2;
	
	private ParamBoolean	skip0Param;
	private ParamBoolean	noiterParam;
	
	private		static final String[]	maskings = {"binary","proba"};
	private		String		masking = "binary";
	
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inv2Image = new ParamVolume("Second inversion (Inv2) Image"));
		
		inputParams.add(t1mapImage = new ParamVolume("Quantitative T1 Map (T1_Images) Image"));
		inputParams.add(isoImage = new ParamVolume("T1-weighted (UNI) Image"));
		t1mapImage.setMandatory(false);
		isoImage.setMandatory(false);
		
		inputParams.add(distribParam = new ParamOption("Background distribution", distribs));
		distribParam.setDescription("Model distribution for background noise (default is half-normal, exponential is more stringent).");
		distribParam.setValue(distrib);
		
		inputParams.add(skip0Param = new ParamBoolean("Skip zero values", true));
		inputParams.add(noiterParam = new ParamBoolean("non-iterative estimate", false));
		
		inputParams.add(maskingParam = new ParamOption("Masking method", maskings));
		maskingParam.setDescription("Whether to use a binary threshold or a weighted average based on the probability.");
		maskingParam.setValue(masking);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity.MP2RAGE");
		inputParams.setLabel("MP2RAGE Background Masking");
		inputParams.setName("MP2RAGEBackgroundMasking");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Estimate a background signal mask for a MP2RAGE dataset.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(probaImage = new ParamVolume("Signal Proba Image",VoxelType.FLOAT));
		outputParams.add(maskImage = new ParamVolume("Signal Mask Image",VoxelType.UBYTE));
		outputParams.add(t1maskImage = new ParamVolume("Masked T1_Map Image",VoxelType.FLOAT));
		outputParams.add(isomaskImage = new ParamVolume("Masked T1-weighted Image",VoxelType.FLOAT));
		
		t1maskImage.setMandatory(false);
		isomaskImage.setMandatory(false);
		
		outputParams.setName("masking images");
		outputParams.setLabel("masking images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		int inv2 = 0;
		int nimg = 1;
		
		int t1map = -1, iso = -1, pv = -1;
		if (t1mapImage.getImageData() != null) { t1map = nimg; nimg++; }
		if (isoImage.getImageData() != null) { iso = nimg; nimg++; }
		
		ImageDataFloat invImg = new ImageDataFloat(inv2Image.getImageData());
		ImageDataFloat t1Img = null, isoImg = null;
		
		if (t1map>-1) t1Img = new ImageDataFloat(t1mapImage.getImageData());
		if (iso>-1) isoImg = new ImageDataFloat(isoImage.getImageData());
		
		float[][][] buffer = invImg.toArray3d();
		int nx = invImg.getRows();
		int ny = invImg.getCols();
		int nz = invImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = invImg.getHeader().getDimResolutions()[0];
		float ry = invImg.getHeader().getDimResolutions()[1];
		float rz = invImg.getHeader().getDimResolutions()[2];
		
		float[][] image = new float[nimg][nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[inv2][xyz] = buffer[x][y][z];
		}
		invImg = null;
		buffer = null;
		
		if (t1map>-1) {
			buffer = t1Img.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[t1map][xyz] = buffer[x][y][z];
			}
			t1Img = null;
			buffer = null;
		}			
		if (iso>-1) {
			buffer = isoImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[iso][xyz] = buffer[x][y][z];
			}
			isoImg = null;
			buffer = null;
		}			
		
		// main algorithm

		// assume  an exponential + outlier distributions for the inv2 image
		
		// use default parameters: 
		float maxdiff = 0.0001f;
		int itermax = 20;
		
		// min max for t1 map fixing
		float max = -1e10f;
		float min = 1e10f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (image[inv2][xyz]>max) max = image[inv2][xyz];
			if (image[inv2][xyz]<min) min = image[inv2][xyz];
		}
		// re-normalize
		for (int xyz=0;xyz<nxyz;xyz++) {
			image[inv2][xyz] = (image[inv2][xyz]-min)/(max-min);
		}
		System.out.println("image range: "+min+", "+max);

		boolean skip0 = skip0Param.getValue().booleanValue();
		
		double mean = 0.0f;
		double den = 0.0;
		for (int xyz=0;xyz<nxyz;xyz++) if (!skip0 || image[inv2][xyz]!=0) {
			mean += image[inv2][xyz];
			den++;
		}
		mean /= den;
		System.out.println("mean parameters: "+mean);
		
		double pi = Math.PI;
		
		byte model = HNORM;
		if (distribParam.getValue().equals("exponential")) model = EXP;
		
		// re-normalized probability map
		float[] proba = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (model==HNORM)
				// half-normal model
				proba[xyz] = (float)(2.0/(mean*pi)*Math.exp(-image[inv2][xyz]*image[inv2][xyz]/(pi*mean*mean)));
			else if (model==EXP)
				// exponential model
				proba[xyz] = (float)(Math.exp(-image[inv2][xyz]/mean)/mean);
			// bg vs. outlier test : proba is p(xyz is background)
			proba[xyz] = proba[xyz]/(1.0f+proba[xyz]);
		}
		// loop
		if (noiterParam.getValue().booleanValue()) itermax = 0;
		double diff = 1.0;
		double prev;
		for (int t=0;t<itermax && diff>maxdiff;t++) {
			System.out.println("iteration "+(t+1));
			prev = mean;
			mean = 0.0;
			den = 0.0;
			for (int xyz=0;xyz<nxyz;xyz++) if (!skip0 || image[inv2][xyz]!=0) {
				mean += proba[xyz]*image[inv2][xyz];
				den += proba[xyz];
			}
			mean /= den;
			System.out.println("mean parameters: "+mean);
			diff = Numerics.abs(prev-mean)/(0.5f*(prev+mean));
			System.out.println("diff parameters: "+diff);
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (model==HNORM)
					// half-normal model
					proba[xyz] = (float)(2.0/(mean*pi)*Math.exp(-image[inv2][xyz]*image[inv2][xyz]/(pi*mean*mean)));
				else if (model==EXP)
					// exponential model
					proba[xyz] = (float)(Math.exp(-image[inv2][xyz]/mean)/mean);
				// bg vs. outlier test : proba is p(xyz is background)
				proba[xyz] = proba[xyz]/(1.0f+proba[xyz]);
			}
		}
		
		byte[][][] mask = new byte[nx][ny][nz];
		float[][][] prob = new float[nx][ny][nz];
			
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (proba[xyz]<=0.5f) mask[x][y][z] = 1;
			else mask[x][y][z] = 0;
			prob[x][y][z] = 1.0f-proba[xyz];
		}
		proba = null;

		ImageDataFloat bufferData = new ImageDataFloat(prob);
		bufferData.setHeader(inv2Image.getImageData().getHeader());
		bufferData.setName(inv2Image.getImageData().getName()+"_proba");
		probaImage.setValue(bufferData);
		bufferData = null;
		
		ImageDataUByte maskData = new ImageDataUByte(mask);		
		maskData.setHeader(inv2Image.getImageData().getHeader());
		maskData.setName(inv2Image.getImageData().getName()+"_mask");
		maskImage.setValue(maskData);
		maskData = null;
		
		if (t1map>-1) {
			buffer = new float[nx][ny][nz];
			if (maskingParam.getValue().equals("proba")) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					buffer[x][y][z] = prob[x][y][z]*image[t1map][xyz];
				}
			} else if (maskingParam.getValue().equals("dilated")) {
				byte[][][] dilated = ObjectMorphology.dilateObject(mask, nx, ny, nz, 2, 2, 2);
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (dilated[x][y][z]>0) buffer[x][y][z] = image[t1map][xyz];
					else buffer[x][y][z] = 0.0f;
				}
				dilated = null;
			} else {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (mask[x][y][z]>0) buffer[x][y][z] = image[t1map][xyz];
					else buffer[x][y][z] = 0.0f;
				}
			}
			bufferData = new ImageDataFloat(buffer);
			bufferData.setHeader(t1mapImage.getImageData().getHeader());
			bufferData.setName(t1mapImage.getImageData().getName()+"_masked");
			t1maskImage.setValue(bufferData);
			bufferData = null;
			buffer = null;
		}
		
		if (iso>-1) {
			buffer = new float[nx][ny][nz];
			if (maskingParam.getValue().equals("proba")) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					buffer[x][y][z] = prob[x][y][z]*image[iso][xyz];
				}
			} else if (maskingParam.getValue().equals("dilated")) {
				byte[][][] dilated = ObjectMorphology.dilateObject(mask, nx, ny, nz, 2, 2, 2);
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (dilated[x][y][z]>0) buffer[x][y][z] = image[iso][xyz];
					else buffer[x][y][z] = 0.0f;
				}
				dilated = null;
			} else {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (mask[x][y][z]>0) buffer[x][y][z] = image[iso][xyz];
					else buffer[x][y][z] = 0.0f;
				}
			}
			bufferData = new ImageDataFloat(buffer);
			bufferData.setHeader(isoImage.getImageData().getHeader());
			bufferData.setName(isoImage.getImageData().getName()+"_masked");
			isomaskImage.setValue(bufferData);
			bufferData = null;
			buffer = null;
		}
		mask = null;
		
	}


}
