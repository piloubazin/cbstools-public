package de.mpg.cbs.jist.laminar;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamPointDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
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
import de.mpg.cbs.core.laminar.LaminarProfileGeometry;

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarProfileGeometry extends ProcessingAlgorithm {

    private LaminarProfileGeometry algorithm;
    
	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume 	layersImage;
	private ParamFloat		extraParam;
	private	ParamFloat		smoothingParam;
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
						    
	private static final String[] suffixTypes = {"th", 
						    "cur", 
						    "si",
						    "mc", 
						    "gc",
						    "pl", 
						    "pc", 
						    "pt",
						    "pd"};
						    
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
		inputParams.add(smoothingParam=new ParamFloat("smoothing parameter",0,5,0.3f));
		
		inputParams.add(extraParam=new ParamFloat("outside extension (mm)",0.0f,10.0f,0.0f));
		
		algorithm = new LaminarProfileGeometry();
		
		inputParams.setPackage(algorithm.getPackage());
		inputParams.setCategory(algorithm.getCategory());
		inputParams.setLabel(algorithm.getLabel());
		inputParams.setName(algorithm.getName());

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(References.getAuthor(algorithm.getAlgorithmAuthors()[0]));
		info.setAffiliation(algorithm.getAffiliation());
		info.setDescription(algorithm.getDescription());
		
		info.setVersion(algorithm.getVersion());
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
		
		// i/o variables
		String name = Interface.getName(layersImage);
		ImageHeader header = Interface.getHeader(layersImage);
		int[] dims = Interface.getDimensions4D(layersImage);
		float[] res = Interface.getResolutions(layersImage);
		
		// main algorithm
		algorithm = new LaminarProfileGeometry();

		algorithm.setLayerSurfaceImage(Interface.getFloatImage4D(layersImage));
		
		algorithm.setComputedMeasure(calcParam.getValue());
		algorithm.setRegularization(regParam.getValue());
		algorithm.setSmoothing(smoothingParam.getValue().floatValue());
		algorithm.setOutsideExtension_mm(extraParam.getValue().floatValue());
		
		algorithm.setLayersDimensions(dims);
		algorithm.setLayersResolutions(res);

		algorithm.execute();
		
		// output
		String suffix = "_lgeom"+suffixTypes[calcParam.getIndex()];
		
		if (algorithm.isResult4D()) {
            Interface.setFloatImage4D(algorithm.getResultImage(), dims, 3, calcImage, name+suffix, header);
        } else {
           Interface.setFloatImage3D(algorithm.getResultImage(), dims, calcImage, name+suffix, header);
        }

	}

}
