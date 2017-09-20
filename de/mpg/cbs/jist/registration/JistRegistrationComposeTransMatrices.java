package de.mpg.cbs.jist.registration;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

//import edu.jhu.ece.iacl.algorithms.PrinceGroupAuthors;
import edu.jhu.ece.iacl.algorithms.registration.RegistrationUtilities;
import edu.jhu.ece.iacl.algorithms.registration.RegistrationUtilities.InterpolationType;
import edu.jhu.ece.iacl.jist.io.ArrayDoubleMtxReaderWriter;
import edu.jhu.ece.iacl.jist.io.ArrayDoubleReaderWriter;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.io.ImageDataReaderWriter;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFileCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamMatrix;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamModel;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.utility.JistLogger;


/*
 * @author Min Chen (mchen55@jhu.edu)
 * @author Pierre-Louis Bazin (bazin@cbs.mpg.de)
 *
 */
public class JistRegistrationComposeTransMatrices extends ProcessingAlgorithm{
	//Inputs
	public ParamMatrix inParamTransform;
	public ParamMatrix inParamLeftMatrix,inParamRightMatrix;
	public ParamModel test;
	public ParamBoolean invertTransform;
	public ParamBoolean invertLeft;
	public ParamBoolean invertRight;
	
	//Outputs
	public ParamMatrix outParamTransform;

	
	//Internal Variables
	int XN, YN, ZN; 
	int chN = 3;
	File dir;
	ArrayDoubleMtxReaderWriter mtxRW;

	
	//Other Variables
	private static final String revnum = RegistrationUtilities.getVersion();
	private static final String shortDescription = "Multiplies a transformation matrices by a left or right matrix.";
	private static final String longDescription = "Multiplies a transformation matrices by a left or right matrix, with or without inversion.";


	protected void createInputParameters(ParamCollection inputParams) {
		//Set initial matrices as identity.(You do need two separate ones or setting one will affect the other).
		Matrix identityLeftMtx = new Matrix(4,4);
		Matrix identityRightMtx = new Matrix(4,4);
		for(int i = 0; i < 4; i++){
			identityLeftMtx.set(i,i,1);
			identityRightMtx.set(i,i,1);
		}
		
		inputParams.add(inParamLeftMatrix=new ParamMatrix("Matrix to Left Multiply By (Optional)",identityLeftMtx));
		inputParams.add(invertLeft=new ParamBoolean("Invert left matrix?",false));
		
		inputParams.add(inParamTransform=new ParamMatrix("Transformations Matrix to Multiply",4,4));
		inputParams.add(invertTransform=new ParamBoolean("Invert transformation matrix?",false));
		
		inputParams.add(inParamRightMatrix=new ParamMatrix("Matrix to Right Multiply By (Optional)",identityRightMtx));
		inputParams.add(invertRight=new ParamBoolean("Invert right matrix?",false));
		
		//inputParams.add(inParamMatrix = new ParamMatrix("Transformation Matrix", 4, 4));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration");
		inputParams.setLabel("Compose Transformation Matrices");
		inputParams.setName("ComposeTransMatrix");


		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Min Chen", "","iacl.ece.jhu.edu"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription(shortDescription);
		info.setLongDescription(shortDescription + longDescription);
		info.setVersion(revnum);
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.BETA);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outParamTransform=new ParamMatrix("Multiplied Transformation",4,4));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ExecuteWrapper wrapper=new ExecuteWrapper();
		monitor.observe(wrapper);
		mtxRW = new ArrayDoubleMtxReaderWriter();
		dir = new File(this.getOutputDirectory()+File.separator+edu.jhu.ece.iacl.jist.utility.FileUtil.forceSafeFilename(this.getAlgorithmName()));
		try{
			if(!dir.isDirectory()){
				(new File(dir.getCanonicalPath())).mkdir();
			}
		}catch(IOException e){
			e.printStackTrace();
		}
		wrapper.execute(this);
	}


	protected class ExecuteWrapper extends AbstractCalculation{
		public void execute(ProcessingAlgorithm alg){
			this.setLabel("combine Volumes");
			
			Matrix finalMatrix;
			Matrix leftMatrix = inParamLeftMatrix.getValue();
			Matrix rightMatrix = inParamRightMatrix.getValue();
			Matrix transMatrix = inParamTransform.getValue();
			
			//set as identify if null, check for inversion
			if(leftMatrix == null) leftMatrix = Matrix.identity(4,4);
			else if (invertLeft.getValue().booleanValue()) leftMatrix = leftMatrix.inverse();
			
			if(rightMatrix == null) rightMatrix = Matrix.identity(4,4);
			else if (invertRight.getValue().booleanValue()) rightMatrix = rightMatrix.inverse();
			
			if (invertTransform.getValue().booleanValue()) transMatrix = transMatrix.inverse();
			
			// compute the final matrix
			finalMatrix = leftMatrix.times(transMatrix.times(rightMatrix));
			
			outParamTransform.setValue(finalMatrix);
				
		}
	}
		
}