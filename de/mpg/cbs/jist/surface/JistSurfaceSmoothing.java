package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
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

import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;
import org.apache.commons.math3.util.FastMath;

public class JistSurfaceSmoothing  extends ProcessingAlgorithm{
	ParamSurface 	inputSurface;
	ParamInteger 	iterationsParam;
	ParamFloat		lambdaParam;
	
	ParamSurface 	outputSurface;
	    
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputSurface = new ParamSurface("Input Surface"));
		inputParams.add(iterationsParam = new ParamInteger("Number of smoothing iterations", 0, 100, 0));
		inputParams.add(lambdaParam = new ParamFloat("Lambda", 0.0f, 1.0f, 0.75f));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setName("SmoothSurfaceMesh");
		inputParams.setLabel("Smooth Surface Mesh");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Smooths a surface mesh by taking the average coordinates of a vertices neighbours");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputSurface=new ParamSurface("Output Surface"));
	}
	
	protected void execute(CalculationMonitor monitor) {
		
		EmbeddedSurface surf = inputSurface.getSurface();
		EmbeddedSurface smoothsurf = inputSurface.getSurface();

		surf.buildAllTables();
		float lambda = lambdaParam.getValue().floatValue();
		
		int[][] neighbors = surf.getNeighborVertexVertexTable();
		
		for (int i=0; i<iterationsParam.getValue().intValue(); i++) {
			for(int v=0; v<surf.getVertexCount(); v++) {
			
				Point3f pt = surf.getVertex(v);
				pt.scale(1.0f-lambda);			
				int N = neighbors[v].length;
				for (int n=0; n<N; n++) {
					Point3f tmp = surf.getVertex(neighbors[v][n]);
					tmp.scale(lambda/N);
					pt.add(tmp);
				}
				smoothsurf.setVertex(v, pt);
			}
			surf = smoothsurf.clone();
		}
		
		smoothsurf.setName(surf.getName()+"_smooth");
		outputSurface.setValue(smoothsurf);
		
	}


}
