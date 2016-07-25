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
	
	private ParamOption corrParam;
	private static final String[] corrTypes = {"absolute","positive","raw"}; 
	
	private static final byte X=0;
	private static final byte Y=1;
	private static final byte Z=2;
	
	private static final byte RAW=30;
	private static final byte ABS=31;
	private static final byte POS=32;
	private static final byte NEG=33;
	
	// ouput
	private ParamVolume pboundaryImage;
	private ParamVolume sboundaryImage;
	private ParamVolume vboundaryImage;
		
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(dataImage = new ParamVolume("fMRI Data Image",VoxelType.FLOAT,-1,-1,-1,-1));
		dataImage.setLoadAndSaveOnValidate(false);
		
		inputParams.add(corrParam = new ParamOption("Correlation type", corrTypes));
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
		outputParams.add(pboundaryImage = new ParamVolume("Partial Boundary Map",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(sboundaryImage = new ParamVolume("Superres Boundary Map",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(vboundaryImage = new ParamVolume("Boundary Vector Map",VoxelType.FLOAT,-1,-1,-1,-1));
		
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

		byte corrtype;
		if (corrParam.getValue().equals("absolute")) corrtype = ABS;
		else if (corrParam.getValue().equals("positive")) corrtype = POS;
		else if (corrParam.getValue().equals("negative")) corrtype = NEG;
		else corrtype = RAW; 
		
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
			norm[x][y][z] = (float)FastMath.sqrt(var/(nt-1.0));
		}
		
		// mask from zero norm
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			mask[x][y][z] = (norm[x][y][z]>0);
		}
		
		float[][][][] partialboundaries = new float[nx][ny][nz][26];
		float[][][] superboundaries = new float[2*nx-1][2*ny-1][2*nz-1];
		float[][][][] boundaryvector = new float[nx][ny][nz][3];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			for (int n=0;n<26;n++) {
				if (mask[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]] && mask[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]]) {
					// partial correlation formula?
					double corrXY = 0.0;
					for (int t=0;t<nt;t++) corrXY += data[x][y][z][t]*data[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]][t];
					corrXY /= nt*norm[x][y][z]*norm[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]];

					if (corrtype==ABS) corrXY = Numerics.abs(corrXY);
					else if (corrtype==POS) corrXY = Numerics.max(corrXY,0.0);
					else if (corrtype==NEG) corrXY = -Numerics.min(corrXY,0.0);
					
					double corrYZ = 0.0;
					for (int t=0;t<nt;t++) corrYZ += data[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]][t]*data[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]][t];
					corrYZ /= nt*norm[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]]*norm[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]];
					
					if (corrtype==ABS) corrYZ = Numerics.abs(corrYZ);
					else if (corrtype==POS) corrYZ = Numerics.max(corrYZ,0.0);
					else if (corrtype==NEG) corrYZ = -Numerics.min(corrYZ,0.0);
					
					double corrZX = 0.0;
					for (int t=0;t<nt;t++) corrZX += data[x][y][z][t]*data[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]][t];
					corrZX /= nt*norm[x][y][z]*norm[x-Ngb.x[n]][y-Ngb.y[n]][z-Ngb.z[n]];
					
					if (corrtype==ABS) corrZX = Numerics.abs(corrZX);
					else if (corrtype==POS) corrZX = Numerics.max(corrZX,0.0);
					else if (corrtype==NEG) corrZX= -Numerics.min(corrZX,0.0);
					
					double pcorr = (corrXY-corrZX*corrYZ)/FastMath.sqrt((1.0-corrZX*corrZX)*(1.0-corrYZ*corrYZ));
					partialboundaries[x][y][z][n] = (float)(pcorr);
					
					superboundaries[2*x-Ngb.x[n]][2*y-Ngb.y[n]][2*z-Ngb.z[n]] = Numerics.max(superboundaries[2*x-Ngb.x[n]][2*y-Ngb.y[n]][2*z-Ngb.z[n]], (float)pcorr);
					
					float[] v = Ngb.directionVector(n);
					boundaryvector[x][y][z][X] += (float)pcorr*v[X];
					boundaryvector[x][y][z][Y] += (float)pcorr*v[Y];
					boundaryvector[x][y][z][Z] += (float)pcorr*v[Z];
				}
			}
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			superboundaries[2*x][2*y][2*z] = superboundaries[2*x-Ngb.x[0]][2*y-Ngb.y[0]][2*z-Ngb.z[0]];
			for (int n=1;n<26;n++) {
				superboundaries[2*x][2*y][2*z] = Numerics.min(superboundaries[2*x][2*y][2*z], superboundaries[2*x-Ngb.x[n]][2*y-Ngb.y[n]][2*z-Ngb.z[n]]);
			}
		}
		System.out.println("output..");
		
		// output
		ImageDataFloat pboundaryImg = new ImageDataFloat(partialboundaries);		
		pboundaryImg.setHeader(dataHeader);
		pboundaryImg.setName(dataName + "_pbound");
		pboundaryImage.setValue(pboundaryImg);
		ImageDataFloat sboundaryImg = new ImageDataFloat(superboundaries);		
		sboundaryImg.setHeader(dataHeader);
		sboundaryImg.setName(dataName + "_sbound");
		sboundaryImage.setValue(sboundaryImg);
		ImageDataFloat vboundaryImg = new ImageDataFloat(boundaryvector);		
		vboundaryImg.setHeader(dataHeader);
		vboundaryImg.setName(dataName + "_pvec");
		vboundaryImage.setValue(vboundaryImg);
	}

}
