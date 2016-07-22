package de.mpg.cbs.jist.fmri;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
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

import org.apache.commons.math3.stat.descriptive.rank.*;
import org.apache.commons.math3.util.*;
import org.apache.commons.math3.special.*;

import gov.nih.mipav.model.file.FileInfoBase;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.model.structures.ModelStorageBase;

/*
 * @author Pierre-Louis Bazin
 */
public class JistFmriPartialCorrelationBoundaries extends ProcessingAlgorithm {

	// jist containers
	
	// input 
	private ParamVolume 	dataImage;
	//private ParamFloat 		extentParam;
	
	//private ParamOption outParam;
	//private static final String[] outTypes = {"max_boundary","resampled"}; 
	
	private static final byte X=0;
	private static final byte Y=1;
	private static final byte Z=2;
	
	// ouput
	private ParamVolume boundaryImage;
		
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(dataImage = new ParamVolume("fMRI Data Image",VoxelType.FLOAT,-1,-1,-1,-1));
		dataImage.setLoadAndSaveOnValidate(false);
		
		//inputParams.add(ngbParam = new ParamOption("Neighborhood type", ngbTypes));
		//inputParams.add(extentParam = new ParamFloat("Neighborhood distance (mm)", 0.0f, 30.0f, 2.0f));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("fMRI");
		inputParams.setLabel("fMRI Partial Correlation Boundaries");
		inputParams.setName("fMRIPartialCorrelationBoundaries");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Estimates local correlation boundaries for fMRI data.");
		
		info.setVersion("3.0.9");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(boundaryImage = new ParamVolume("Boundary Map",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setLabel("fMRI Correlation Boundaries");
		outputParams.setName("fMRICorrelationBoundaries");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat	dataImg = new ImageDataFloat(dataImage.getImageData());
		int nx = dataImg.getRows();
		int ny = dataImg.getCols();
		int nz = dataImg.getSlices();
		int nt = dataImg.getComponents();
		float rx = dataImg.getHeader().getDimResolutions()[0];
		float ry = dataImg.getHeader().getDimResolutions()[1];
		float rz = dataImg.getHeader().getDimResolutions()[2];
		
		String dataName = dataImg.getName();
		ImageHeader dataHeader = dataImg.getHeader();
		
		float[][][][] data = dataImg.toArray4d();
		dataImg.dispose();
		dataImage.dispose();
		
		System.out.println("fmri data loaded");

		// center and norm
		float[][][] center = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			double mean = 0.0;
			for (int t=0;t<nt;t++) mean += data[x][y][z][t];
			center[x][y][z] = (float)(mean/nt);
		}
		float[][][] norm = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			for (int t=0;t<nt;t++) data[x][y][z][t] -= center[x][y][z];
			double var = 0.0;
			for (int t=0;t<nt;t++) var += data[x][y][z][t]*data[x][y][z][t];
			norm[x][y][z] = (float)FastMath.sqrt(var/(nt-1));
		}
		
		// mask from zero norm
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			mask[x][y][z] = (norm[x][y][z]>0);
		}
		
		float[][][][] boundaries = new float[nx][ny][nz][26];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			for (int n=0;n<26;n++) {
				if (mask[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]] && mask[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]]) {
					// partial correlation formula?
					double corrXY = 0.0;
					for (int t=0;t<nt;t++) corrXY += data[x][y][z][t]*data[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]][t];
					corrXY /= nt*norm[x][y][z]*norm[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]];

					double corrXZ = 0.0;
					for (int t=0;t<nt;t++) corrXZ += data[x][y][z][t]*data[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]][t];
					corrXZ /= nt*norm[x][y][z]*norm[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]];
					
					double corrYZ = 0.0;
					for (int t=0;t<nt;t++) corrYZ += data[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]][t]*data[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]][t];
					corrYZ /= nt*norm[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]]*norm[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]];
					
					boundaries[x][y][z][n] = (float)((corrXY-corrXZ*corrYZ)/(FastMath.sqrt((1.0-corrXZ*corrXZ)*(1.0-corrYZ*corrYZ))));
				}
			}
		}
		
		System.out.println("output..");
		
		// output
		ImageDataFloat boundaryImg = new ImageDataFloat(boundaries);		
		boundaryImg.setHeader(dataHeader);
		boundaryImg.setName(dataName + "_pcbound");
		boundaryImage.setValue(boundaryImg);
	}

}
