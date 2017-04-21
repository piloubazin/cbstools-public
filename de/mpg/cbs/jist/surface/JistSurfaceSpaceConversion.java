package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;

import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.graphics.smooth.SurfaceInflate;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.SurfaceToMask;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import javax.vecmath.Point3f;

/*
 * @author Pierre-Louis Bazin and Christine L Tardif
 */
public class JistSurfaceSpaceConversion extends ProcessingAlgorithm {

	// jist containers
	ParamSurface	inputSurface;
	ParamVolume		referenceVolume;
	//ParamBoolean 	mipavTransform;
	
	ParamSurface	outputSurface;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputSurface = new ParamSurface("Input Surface"));
		inputParams.add(referenceVolume = new ParamVolume("Reference Volume"));
		//inputParams.add(mipavTransform = new ParamBoolean("Align mesh to MIPAV image space", true));		
				
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setLabel("Surface Space Conversion");
		inputParams.setName("SurfaceSpaceConversion");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Transforms the surface mesh into a different space.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputSurface = new ParamSurface("Converted Surface"));

		outputParams.setName("converted surface mesh");
		outputParams.setLabel("converted surface mesh");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		EmbeddedSurface surf = inputSurface.getSurface();
		
		ImageDataFloat	refVol = new ImageDataFloat(referenceVolume.getImageData());
		int nx = refVol.getRows();
		int ny = refVol.getCols();
		int nz = refVol.getSlices();
		float rx = refVol.getHeader().getDimResolutions()[0];
		float ry = refVol.getHeader().getDimResolutions()[1];
		float rz = refVol.getHeader().getDimResolutions()[2];
		
		surf.scaleVertices(new float[]{1.0f, -1.0f, -1.0f});
		Point3f pt = new Point3f(0, (ny-1)*ry, (nz-1)*rz);
		surf.translate(pt);
		
		//output
		outputSurface.setValue(surf);
	}


}
