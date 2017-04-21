package de.mpg.cbs.jist.laminar;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurfaceCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import de.mpg.cbs.core.laminar.*;

import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3f;

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarProfileMeshing extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 			layersImage;
	private ParamSurface 			inputSurface;

	private ParamSurfaceCollection 	sampledSurfaces;
	
	private LaminarProfileMeshing algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputSurface=new ParamSurface("Input Surface"));
		inputParams.add(layersImage = new ParamVolume("Profile Surface Image (4D)"));
		
		algorithm = new LaminarProfileMeshing();
		
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
		outputParams.add(sampledSurfaces=new ParamSurfaceCollection("Sampled Surfaces"));
		
		outputParams.setName("sampled surfaces");
		outputParams.setLabel("sampled surfaces");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.getName(inputSurface);
		ImageHeader header = Interface.getHeader(layersImage);
		int[] dims = Interface.getDimensions(layersImage);
		float[] res = Interface.getResolutions(layersImage);
		int nlayers = Interface.getComponents(layersImage);
		
		// main algorithm
		algorithm = new LaminarProfileMeshing();
		
		algorithm.setProfileSurfaceImage(Interface.getFloatImage4D(layersImage));
		algorithm.setInputSurfacePoints(Interface.getSurfacePoints(inputSurface));
		algorithm.setInputSurfaceTriangles(Interface.getSurfaceTriangles(inputSurface));
		
		algorithm.setDimensions(dims[0],dims[1],dims[2],nlayers);
		algorithm.setResolutions(res);

		algorithm.execute();

		for (int n=0;n<nlayers;n++) {
			ParamSurface surf = new ParamSurface("Sampled Surface"+n);
			Interface.setSurface(algorithm.getSampledSurfacePoints(n),algorithm.getSampledSurfaceTriangles(n),surf,name+"_layer"+n);
			sampledSurfaces.add(surf);
		}
	}

}
