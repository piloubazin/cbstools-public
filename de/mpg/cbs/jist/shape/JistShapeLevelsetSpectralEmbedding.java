package de.mpg.cbs.jist.shape;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamMatrix;
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

import Jama.Matrix;
import Jama.EigenvalueDecomposition;
import Jama.QRDecomposition;

/*
 * @author Pierre-Louis Bazin
 */
public class JistShapeLevelsetSpectralEmbedding extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputImage;
	
	private ParamFloat boundParam;
	private ParamInteger connectParam;
	private ParamInteger dimParam;
	private ParamOption featureParam;
	private static final String[] featureTypes = {"distance","gradient","inverse_gradient"}; 
	private ParamOption methodParam;
	private static final String[] methodTypes = {"full_decomposition","partial_QR","power_iteration"}; 
	
	private ParamVolume embeddingImage;
	//private ParamVolume matrixImage;
	private ParamMatrix eigenvalMtx;
	
	// numerical quantities
	private static final	float	INVSQRT2 = (float)(1.0/FastMath.sqrt(2.0));
	private static final	float	INVSQRT3 = (float)(1.0/FastMath.sqrt(3.0));
	private static final	float	SQRT2 = (float)FastMath.sqrt(2.0);
	private static final	float	SQRT3 = (float)FastMath.sqrt(3.0);

	// direction labeling		
	public	static	final	byte	X = 0;
	public	static	final	byte	Y = 1;
	public	static	final	byte	Z = 2;

	// feature choice labeling		
	public	static	final	byte	DIST = 100;
	public	static	final	byte	GRAD = 101;
	public	static	final	byte	INVGRAD = 102;
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;

	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputImage = new ParamVolume("Shape image"));
		inputParams.add(featureParam = new ParamOption("Weighting type", featureTypes));
		
		inputParams.add(boundParam = new ParamFloat("Boundary thickness: ", 0.0f, 100.0f, 0.0f));
		inputParams.add(connectParam = new ParamInteger("Connectivity", 6, 26, 6));
		inputParams.add(dimParam = new ParamInteger("Embedding dimension", 0, 100, 5));
		inputParams.add(methodParam = new ParamOption("Computation method", methodTypes));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Shape.devel");
		inputParams.setLabel("Level set Spectral Embeding");
		inputParams.setName("LevelsetSpectralEmbedding");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Eliza Orasanu", "","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Create a spectral representation for a levelset surface");
		
		info.setVersion("3.0.9");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(embeddingImage = new ParamVolume("Embedding Image",VoxelType.FLOAT,-1,-1,-1,-1));
		//outputParams.add(matrixImage = new ParamVolume("Matrix Image",VoxelType.FLOAT));
		outputParams.add(eigenvalMtx = new ParamMatrix("Eigenvalues", 5, 1));

		outputParams.setName("embedding images");
		outputParams.setLabel("embedding images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		System.out.println("input");
		ImageDataFloat	inputImg = new ImageDataFloat(inputImage.getImageData());
		
		int nx = inputImg.getRows();
		int ny = inputImg.getCols();
		int nz = inputImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = inputImg.getHeader().getDimResolutions()[0];
		float ry = inputImg.getHeader().getDimResolutions()[1];
		float rz = inputImg.getHeader().getDimResolutions()[2];
		
		float[][][] input = inputImg.toArray3d();
		
		byte mode;
		if (featureParam.getValue().equals("distance")) mode = DIST;
		else if (featureParam.getValue().equals("gradient")) mode = GRAD;
		else if (featureParam.getValue().equals("inverse_gradient")) mode = INVGRAD;
		else mode = 0;
		
		float boundary = boundParam.getValue().floatValue();
		int connectivity = connectParam.getValue().intValue();
		int dimensions = dimParam.getValue().intValue();
		
		int datasize=0;
		int[][][] location = new int[nx][ny][nz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (Numerics.abs(input[x][y][z])<=boundary) {
				location[x][y][z] = datasize;
				datasize++;
			} else {
				location[x][y][z] = -1;
			}
		}
		System.out.println("data matrix size: "+datasize+" x "+datasize);
				
		// 1. build a gigantic matrix
		double[][] laplacian = new double[datasize][datasize];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (Numerics.abs(input[x][y][z])<=boundary) {
				int loc = location[x][y][z];
				for (int n=0;n<connectivity;n++) {
					int ngb = location[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]];
					if (ngb!=-1) {
						if (mode==DIST) {
							if (n<6) laplacian[loc][ngb] = 1.0;
							if (n<18) laplacian[loc][ngb] = INVSQRT2;
							if (n<26) laplacian[loc][ngb] = INVSQRT3;
						} else if (mode==GRAD) {
							double weight = Numerics.max(1.0-Numerics.abs(input[x][y][z]-input[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]]),0.0);
							if (n<6) laplacian[loc][ngb] = weight;
							if (n<18) laplacian[loc][ngb] = weight*INVSQRT2;
							if (n<26) laplacian[loc][ngb] = weight*INVSQRT3;
						} else if (mode==INVGRAD) {
							double weight = 1.0/Numerics.max(1e-3,Numerics.abs(input[x][y][z]-input[x+Ngb.x[n]][y+Ngb.y[n]][z+Ngb.z[n]]));
							if (n<6) laplacian[loc][ngb] = weight;
							if (n<18) laplacian[loc][ngb] = weight*INVSQRT2;
							if (n<26) laplacian[loc][ngb] = weight*INVSQRT3;
						}
					}
				}
			}
		}
		// L = G^-1 (D-W), D_ii = \sum_j W_ij; G_ii = 1/\phi_i
		for (int n=0;n<datasize;n++) {
			double sum = 0.0;
			for (int m=0;m<datasize;m++) { 
				sum += laplacian[n][m];
				laplacian[n][m] *= -1.0;
			}
			laplacian[n][n] += sum;
		}
		// G^-1
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (Numerics.abs(input[x][y][z])<=boundary) {
				int loc = location[x][y][z];
				laplacian[loc][loc] *= input[x][y][z];
			}
		}
		Matrix laplacianmatrix = new Matrix(laplacian);
		
		// 2. Get the eigenvalues
		Matrix eigenvectors = null;
		Matrix eigenvalues = null;
		if (methodParam.getValue().equals("full_decomposition")) {
			EigenvalueDecomposition eigen = new EigenvalueDecomposition(laplacianmatrix);
			
			eigenvectors = new Matrix(datasize,dimensions);
			eigenvalues = new Matrix(dimensions,1);
			for (int i = 0; i < datasize; i++) {
				for (int j = 0; j < dimensions; j++) {
					eigenvectors.set(i,j, eigen.getV().get(i,j));
				}
			}
			for (int j = 0; j < dimensions; j++) {
				eigenvalues.set(j,1, eigen.getRealEigenvalues()[j]);
			}
		} else if (methodParam.getValue().equals("partial_QR")) {
			// initial basis: random vector
			Matrix basis = Matrix.random(datasize,dimensions);
			Matrix qbasis = basis.qr().getQ();
			int tmax = 100;
			for (int t=0;t<tmax;t++) {
				// power iteration
				basis = laplacianmatrix.times(qbasis);
				// orthogonalization
				QRDecomposition qr = basis.qr();
				basis = qr.getQ();
				// difference
				double diff = 0.0;
				for (int i = 0; i < datasize; i++) {
					for (int j = 0; j < dimensions; j++) {
						diff = Numerics.max(diff, Numerics.abs(basis.get(i,j) - qbasis.get(i,j)));
					}
				}
				if (diff < 1e-3 || t==tmax-1) {
					t = tmax;
					eigenvectors = basis;
					eigenvalues = new Matrix(dimensions,1);
					Matrix r = qr.getR();
					for (int j = 0; j < dimensions; j++) {
						eigenvalues.set(j,1, r.get(j,j));
					}
				} else {
					qbasis = basis;
				}
			}
 		}
		// 3. Project back into image space
		float[][][][] embedding = new float[nx][ny][nz][dimensions];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (Numerics.abs(input[x][y][z])<=boundary) {
				int loc = location[x][y][z];
				for (int d=0;d<dimensions;d++) {
					embedding[x][y][z][d] = (float)eigenvectors.get(loc,d);
				}
			}
		}
		
		ImageDataFloat embeddingData = new ImageDataFloat(embedding);		
		embeddingData.setHeader(inputImage.getImageData().getHeader());
		embeddingData.setName(inputImage.getImageData().getName()+"_emb");
		embeddingImage.setValue(embeddingData);
		embeddingData = null;
		embedding = null;
		
		eigenvalMtx.setValue(eigenvalues);
		/*
		ImageDataFloat matrixData = new ImageDataFloat(matrix);		
		matrixData.setHeader(inputImage.getImageData().getHeader());
		matrixData.setName(inputImage.getImageData().getName()+"_mtx");
		matrixImage.setValue(matrixData);
		matrixData = null;
		matrix = null;
		*/
	}

}
