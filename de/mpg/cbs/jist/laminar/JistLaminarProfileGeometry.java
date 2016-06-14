package de.mpg.cbs.jist.laminar;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamPointDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
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

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarProfileGeometry extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume 	layersImage;
	private ParamDouble		extraParam;
	private	ParamDouble		smoothingParam;
	private ParamOption		calcParam;
	private ParamOption		regParam;
	
	private static final String[] calcTypes = {"thickness", 
						    "curvedness", 
						    "shape_index",
						    "mean_curvature", 
						    "gauss_curvature",
						    "profile_length", 
						    "profile_curvature", 
						    "profile_torsion",
						    "profile_direction"};
						    
	private static final String[] regTypes = {"none","Gaussian"};
	
	private ParamVolume calcImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private int	nx, ny, nz;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(layersImage = new ParamVolume("Profile Surface Image",null,-1,-1,-1,-1));
		
		inputParams.add(calcParam = new ParamOption("computed measure", calcTypes));
		
		inputParams.add(regParam = new ParamOption("regularization", regTypes));
		inputParams.add(smoothingParam=new ParamDouble("smoothing parameter",0,5,0.3));
		
		inputParams.add(extraParam=new ParamDouble("outside extension (mm)",0.0,10.0,0.0));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis");
		inputParams.setLabel("Profile Geometry");
		inputParams.setName("ProfileGeometry");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute various geometric quantities for a cortical layers.");
		
		info.setVersion("3.0.7");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(calcImage = new ParamVolume("Result",VoxelType.FLOAT,-1,-1,-1,-1));

		outputParams.setName("profile geometry");
		outputParams.setLabel("profile geometry");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays : TO DO
		
		ImageDataFloat	layersImg = new ImageDataFloat(layersImage.getImageData());
		
		nx = layersImg.getRows();
		ny = layersImg.getCols();
		nz = layersImg.getSlices();
		int nlayers = layersImg.getComponents()-1;
		int nxyz = nx*ny*nz;
		
		// note: we assume isotropic data
		float rx = layersImg.getHeader().getDimResolutions()[0];
		float ry = layersImg.getHeader().getDimResolutions()[1];
		float rz = layersImg.getHeader().getDimResolutions()[2];
		
		float[][] layers = new float[nlayers+1][nxyz];
		if (nlayers==0) {
			float[][][] buffer3 = layersImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				layers[0][xyz] = buffer3[x][y][z];
			}
			buffer3 = null;
		} else {
			float[][][][] buffer4 = layersImg.toArray4d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
				int xyz = x+nx*y+nx*ny*z;
				layers[l][xyz] = buffer4[x][y][z][l];
			}
			buffer4 = null;
		}
		layersImg = null;
		
		// main algorithm
		String imgname = layersImage.getImageData().getName();

		// create a mask for regions of value 0 on center image (kinda arbitrary..)
		boolean[] ctxmask = new boolean[nxyz];
		float extra = extraParam.getValue().floatValue()/rx;
		for (int xyz=0;xyz<nxyz;xyz++) {
			ctxmask[xyz] = (layers[0][xyz] >= -extra && layers[nlayers][xyz] <= extra);
		}
		
		// main algorithm: thickness & curve parameters
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		float[][][] calc = null;
		float[][][][] calc4d = null;
		if (calcParam.getValue().equals("profile_direction")) calc4d = new float[nx][ny][nz][3];
		else calc = new float[nx][ny][nz];

		// if surface processing: compute the smoothing
		/* in all cases
		if (calcParam.getValue().equals("mean_curvature") 
			|| calcParam.getValue().equals("gauss_curvature") 
			|| calcParam.getValue().equals("shape_index") 
			|| calcParam.getValue().equals("curvedness") ) {
		*/
		if (regParam.getValue().equals("Gaussian")) {
			float kscale = smoothingParam.getValue().floatValue();
			float[][] kernel = ImageFilters.separableGaussianKernel(kscale, kscale, kscale);
			int ks = (kernel[0].length-1)/2;
			
			for (int l=0;l<=nlayers;l++) {
				layers[l] = ImageFilters.separableConvolution(layers[l], nx, ny, nz, kernel, ks, ks, ks);
			}
		}
		//}
		
		for (int x=1; x<nx-1; x++) for (int y=1; y<ny-1; y++) for (int z=1; z<nz-1; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				if (calcParam.getValue().equals("thickness")) {
					// basic thickness definition
					calc[x][y][z] = (layers[0][xyz]-layers[nlayers][xyz])*rx;
				} else
				if (calcParam.getValue().equals("mean_curvature") 
					|| calcParam.getValue().equals("gauss_curvature") 
					|| calcParam.getValue().equals("shape_index") 
					|| calcParam.getValue().equals("curvedness")  
					|| calcParam.getValue().equals("profile_direction") ) {
				
					// search for closest surface to the point
					int lbest = 0;
					for (int l=1;l<=nlayers;l++) {
						if (Numerics.abs(layers[l][xyz])<Numerics.abs(layers[lbest][xyz]) ) {
							lbest = l;
						}
					}
					// compute the corresponding measure
					if (calcParam.getValue().equals("mean_curvature")) {
						calc[x][y][z] = ImageGeometry.meanCurvatureValue(layers[lbest], xyz, nx, ny, nz);
					} else if (calcParam.getValue().equals("gauss_curvature")) {
						calc[x][y][z] = ImageGeometry.gaussCurvatureValue(layers[lbest], xyz, nx, ny, nz);
					} else if (calcParam.getValue().equals("shape_index")) {
						calc[x][y][z] = ImageGeometry.shapeIndexValue(layers[lbest], xyz, nx, ny, nz);
					} else if (calcParam.getValue().equals("curvedness")) {
						calc[x][y][z] = ImageGeometry.curvednessValue(layers[lbest], xyz, nx, ny, nz);
					}else if (calcParam.getValue().equals("profile_direction")) {
						calc4d[x][y][z][X] = 0.5f*(layers[lbest][xyz+1]-layers[lbest][xyz-1]);
						calc4d[x][y][z][Y] = 0.5f*(layers[lbest][xyz+nx]-layers[lbest][xyz-nx]);
						calc4d[x][y][z][Z] = 0.5f*(layers[lbest][xyz+nx*ny]-layers[lbest][xyz-nx*ny]);
						float norm = (float)FastMath.sqrt(calc4d[x][y][z][X]*calc4d[x][y][z][X]+calc4d[x][y][z][Y]*calc4d[x][y][z][Y]+calc4d[x][y][z][Z]*calc4d[x][y][z][Z]);
						if (norm>0) {
							calc4d[x][y][z][X] /= norm; 
							calc4d[x][y][z][Y] /= norm; 
							calc4d[x][y][z][Z] /= norm; 
						}
					}
				} else {
					// first compute the profile	
					float err = profile.computeTrajectory(layers, x, y, z);
					if (calcParam.getValue().equals("profile_length")) {
						calc[x][y][z] = profile.computeLength();
					} else if (calcParam.getValue().equals("profile_curvature")) {
						calc[x][y][z] = profile.computeCurvature();
					} else if (calcParam.getValue().equals("profile_torsion")) {
						calc[x][y][z] = profile.computeTorsion();
					}
				}
			}
		}

		// output
		ImageDataFloat calcData;
		if (calc!=null) calcData = new ImageDataFloat(calc);		
		else calcData = new ImageDataFloat(calc4d);		
		calcData.setHeader(layersImage.getImageData().getHeader());
		calcData.setName(imgname+"_"+calcParam.getValue());
		calcImage.setValue(calcData);
		calcData = null;
		calc = null;
		calc4d = null;
				
		return;			
	}

}
