package de.mpg.cbs.jist.registration;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.jhu.ece.iacl.jist.io.ImageDataReaderWriter;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.PipeAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.factory.ParamFactory;
import edu.jhu.ece.iacl.jist.pipeline.parameter.*;
import edu.jhu.ece.iacl.jist.pipeline.parser.ShellScriptXmlParser;
import edu.jhu.ece.iacl.jist.pipeline.parser.ShellScriptXmlParser.ParamTypes;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistRegistrationEmbeddedSyN extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume sourceImage;
	private ParamVolume targetImage;
	
	private	ParamFile		synScriptParam;
	private ParamInteger 	coarseItParam;
	private ParamInteger 	medItParam;
	private ParamInteger 	fineItParam;
	private ParamBoolean	affineParam;
	private ParamOption		costParam;
	private ParamOption		interpolationParam;
	
	private static final String[] costTypes = {"Cross Correlation","Mutual Information"};
	private static final String[] interpolationTypes = {"Linear", "Nearest Neighbor", "BSpline"};
	
	private ParamVolume deformedImage;
	private ParamVolume mappingImage;
	private ParamVolume inverseMappingImage;
		
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;


	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(sourceImage = new ParamVolume("Source Image"));
		inputParams.add(targetImage = new ParamVolume("Target Image"));
		
		inputParams.add(synScriptParam = new ParamFile("SyN script"));
		
		inputParams.add(coarseItParam = new ParamInteger("coarse level iterations", 0, 1000, 40));
		inputParams.add(medItParam = new ParamInteger("medium level iterations", 0, 1000, 50));
		inputParams.add(fineItParam = new ParamInteger("fine level iterations", 0, 1000, 40));
		inputParams.add(affineParam = new ParamBoolean("run affine first", true));
		inputParams.add(costParam = new ParamOption("cost function",costTypes));
		inputParams.add(interpolationParam = new ParamOption("interpolation method", interpolationTypes));

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration");
		inputParams.setLabel("Embedded SyN");
		inputParams.setName("EmbeddedSyN");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Brian B. Avants", "",""));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new Citation("Avants BB, Epstein CL, Grossman M, Gee JC, "
								+"Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain, "
								+"Med Image Anal. 2008 Feb;12(1):26-41"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Provides a simple wrapper around the SyN algorithm from ANTs.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(deformedImage = new ParamVolume("Deformed source",VoxelType.FLOAT));
		outputParams.add(mappingImage = new ParamVolume("Mapping function",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(inverseMappingImage = new ParamVolume("Inverse Mapping function",VoxelType.FLOAT,-1,-1,-1,-1));
			
		outputParams.setName("deformed images");
		outputParams.setLabel("deformed images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays (we assume source and target are in the same space)
		
		ImageDataFloat	sImg = new ImageDataFloat(sourceImage.getImageData());
		ImageDataFloat	tImg = new ImageDataFloat(targetImage.getImageData());
		
		int nsx = sImg.getRows();
		int nsy = sImg.getCols();
		int nsz = sImg.getSlices();
		ImageHeader sHeader = new ImageHeader();
		//sHeader.setAxisOrientation(sImg.getHeader().getAxisOrientation());
		sHeader.setDimResolutions(sImg.getHeader().getDimResolutions());
		//sHeader.setImageOrientation(sImg.getHeader().getImageOrientation());
		//sHeader.setOrigin(sImg.getHeader().getOrigin());
		//sHeader.setSliceThickness(sImg.getHeader().getSliceThickness());
		//sHeader.setUnitsOfMeasure(sImg.getHeader().getUnitsOfMeasure());
		
		int ntx = tImg.getRows();
		int nty = tImg.getCols();
		int ntz = tImg.getSlices();
		ImageHeader tHeader = new ImageHeader();
		//tHeader.setAxisOrientation(tImg.getHeader().getAxisOrientation());
		tHeader.setDimResolutions(tImg.getHeader().getDimResolutions());
		//tHeader.setImageOrientation(tImg.getHeader().getImageOrientation());
		//tHeader.setOrigin(tImg.getHeader().getOrigin());
		//tHeader.setSliceThickness(tImg.getHeader().getSliceThickness());
		//tHeader.setUnitsOfMeasure(tImg.getHeader().getUnitsOfMeasure());
		
		float[][][] source = sImg.toArray3d();
		sImg = null;
		
		float[][][] target = tImg.toArray3d();
		tImg = null;
		
		BasicInfo.displayMessage("Init deformations\n");

		// init coordinates for SyN processing
		float[][][][] scoord = new float[nsx][nsy][nsz][3];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			scoord[x][y][z][X] = x;		
			scoord[x][y][z][Y] = y;		
			scoord[x][y][z][Z] = z;		
		}
		float[][][][] tcoord = new float[ntx][nty][ntz][3];
		for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			tcoord[x][y][z][X] = x;		
			tcoord[x][y][z][Y] = y;		
			tcoord[x][y][z][Z] = z;		
		}
		
		try {
			ImageDataReaderWriter readerWriter = new ImageDataReaderWriter();
			File destdir = new File(this.getOutputDirectory().getCanonicalFile()+File.separator+this.getAlgorithmName());
			
			// create the directory if it doesn't exist
			if(!destdir.isDirectory()){
				(new File(destdir.getCanonicalPath())).mkdir();
			}
			// write the temporary files to be processed
			String inputSourceLoc = destdir.getAbsolutePath() + File.separator + "input_source.nii";
			File inputSource = new File(inputSourceLoc);
			ImageDataFloat sourceData = new ImageDataFloat(source);
			sourceData.setHeader(sHeader);
			readerWriter.write(sourceData, inputSource);
						
			String inputTargetLoc = destdir.getAbsolutePath() + File.separator + "input_target.nii";
			File inputTarget = new File(inputTargetLoc);
			ImageDataFloat targetData = new ImageDataFloat(target);
			targetData.setHeader(tHeader);
			readerWriter.write(targetData, inputTarget);
						
			String inputSCoordLoc = destdir.getAbsolutePath() + File.separator + "input_scoord.nii";
			File inputSCoord = new File(inputSCoordLoc);
			ImageDataFloat scoordData = new ImageDataFloat(scoord);
			scoordData.setHeader(sHeader);
			readerWriter.write(scoordData, inputSCoord);
			
			String inputTCoordLoc = destdir.getAbsolutePath() + File.separator + "input_tcoord.nii";
			File inputTCoord = new File(inputTCoordLoc);
			ImageDataFloat tcoordData = new ImageDataFloat(tcoord);
			tcoordData.setHeader(tHeader);
			readerWriter.write(tcoordData, inputTCoord);
			
			String outputImageLoc = destdir.getAbsolutePath() + File.separator + "output_image.nii";
			File outputImage = new File(outputImageLoc);
			
			String outputMappingLoc = destdir.getAbsolutePath() + File.separator + "output_map.nii";
			File outputMapping = new File(outputMappingLoc);
			
			String outputInvMappingLoc = destdir.getAbsolutePath() + File.separator + "output_invmap.nii";
			File outputInvMapping = new File(outputInvMappingLoc);
			
			//build command array, size = 1 + # of inputs + # of outputs
			String[] comArray = new String[15];
			comArray[0] = synScriptParam.getValue().getAbsolutePath();
			// first command is shell argument
			comArray[1] = destdir.getAbsolutePath();
			comArray[2] = inputTarget.getAbsolutePath();
			comArray[3] = inputSource.getAbsolutePath();
			comArray[4] = inputTCoord.getAbsolutePath();
			comArray[5] = inputSCoord.getAbsolutePath();
			comArray[6] = Integer.toString(coarseItParam.getValue().intValue());
			comArray[7] = Integer.toString(medItParam.getValue().intValue());
			comArray[8] = Integer.toString(fineItParam.getValue().intValue());
			comArray[9] = Boolean.toString(affineParam.getValue().booleanValue());
			comArray[10] = costParam.getValue();
			comArray[11] = interpolationParam.getValue();
			comArray[12] = outputImage.getAbsolutePath();
			comArray[13] = outputMapping.getAbsolutePath();
			comArray[14] = outputInvMapping.getAbsolutePath();
			
			System.out.println("Current Command:");
			for(int i=0; i<comArray.length; i++){
				System.out.println(comArray[i]);
			}
			System.out.println("\n");
			
			// run the script
			Runtime rt = Runtime.getRuntime();
			Process pr = rt.exec(comArray);

			// collect the command line output
			BufferedReader inputOut = new BufferedReader(new InputStreamReader(pr.getInputStream()));
			BufferedReader inputError = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
			String lineOut=null;
			String lineError=null;
			while((lineOut=inputOut.readLine()) != null || (lineError=inputError.readLine()) != null) {
				System.out.println(lineOut);
				System.err.println(lineError);
			}
			int exitVal = pr.waitFor();
			
			if(exitVal==0){
				System.out.println("Executed successfully");
			}else{
				System.out.println("Executed with error(s)");
			}
			
			// output
			ImageData deformData = readerWriter.read(outputImage);
			deformData.setHeader(sourceImage.getImageData().getHeader());
			deformData.setName(sourceImage.getImageData().getName()+"_def");
			deformedImage.setValue(deformData);
			deformData = null;	
			
			ImageData mappingData = readerWriter.read(outputMapping);
			mappingData.setHeader(targetImage.getImageData().getHeader());
			mappingData.setName(sourceImage.getImageData().getName()+"_map");
			mappingImage.setValue(mappingData);
			mappingData = null;
			
			ImageData invmappingData = readerWriter.read(outputInvMapping);
			invmappingData.setHeader(sourceImage.getImageData().getHeader());
			invmappingData.setName(sourceImage.getImageData().getName()+"_invmap");
			inverseMappingImage.setValue(invmappingData);
			invmappingData = null;
			
			//removes all files created for extension compatibility
			System.out.format("Deleting: " + inputSource.getAbsolutePath() + "\n");
			System.out.format("Deleting: " + inputTarget.getAbsolutePath() + "\n");
			System.out.format("Deleting: " + inputSCoord.getAbsolutePath() + "\n");
			System.out.format("Deleting: " + inputTCoord.getAbsolutePath() + "\n");
			System.out.format("Deleting: " + outputImage.getAbsolutePath() + "\n");
			System.out.format("Deleting: " + outputMapping.getAbsolutePath() + "\n");
			System.out.format("Deleting: " + outputInvMapping.getAbsolutePath() + "\n");
			
			ImageDataReaderWriter.deleteImageFile(inputSource);
			ImageDataReaderWriter.deleteImageFile(inputTarget);
			ImageDataReaderWriter.deleteImageFile(inputSCoord);
			ImageDataReaderWriter.deleteImageFile(inputTCoord);
			ImageDataReaderWriter.deleteImageFile(outputImage);
			ImageDataReaderWriter.deleteImageFile(outputMapping);
			ImageDataReaderWriter.deleteImageFile(outputInvMapping);
			
		} catch(Exception e) {
			System.out.println("Executed with error(s)");
			e.printStackTrace();
		}
		
	}

}
