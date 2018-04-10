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

import org.apache.commons.math3.stat.descriptive.rank.*;
import org.apache.commons.math3.util.*;
import org.apache.commons.math3.special.*;

import edu.jhu.ece.iacl.algorithms.graphics.intersector.SurfaceIntersector;
import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.ResampleLevelSet;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.SurfaceToMask;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.topology.TopologyCorrection;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.*;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;


public class JistSurfaceMeshDataZscore  extends ProcessingAlgorithm{
	private ParamSurface 	inputSurface;
	//ParamFloat		fwhmParam;
	private ParamOption    statParam;
	private		static final String[]         statTypes = {"standard", "robust"};
	
	private ParamSurface 	outputSurface;
	    
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputSurface = new ParamSurface("Input Surface"));
		//inputParams.add(fwhmParam = new ParamFloat("FWHM", 0.0f, 20.0f, 5.0f));
		inputParams.add(statParam = new ParamOption("statistics type: ", statTypes));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces.devel");
		inputParams.setName("ZscoreSurfaceMeshData");
		inputParams.setLabel("Z-score Surface Mesh Data");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Z-score the data mapped onto a surface mesh.");
		
		info.setVersion("3.1.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputSurface=new ParamSurface("Output Surface"));
	}
	
	protected void execute(CalculationMonitor monitor) {
		
		EmbeddedSurface surf = inputSurface.getSurface();
		
		double[][] sdata = surf.getVertexData();
		
		// compute the mean, stdev or med, iqr
		double[] center = new double[sdata[0].length];
		double[] scale = new double[sdata[0].length];
		for (int d=0; d<sdata[0].length; d++) {
		    center[d] = 0.0;
		    scale[d] = 0.0;
		    
		    double nval = 0.0;
            for(int v=0; v<surf.getVertexCount(); v++) if (sdata[v][d]!=0.0) {
                nval++;
            }
            
		    if (statParam.getValue().equals("standard")) {
		        for(int v=0; v<surf.getVertexCount(); v++) if (sdata[v][d]!=0.0) {
                    center[d] += sdata[v][d];
                }
                if (nval>0) center[d] /= nval;
                for(int v=0; v<surf.getVertexCount(); v++) if (sdata[v][d]!=0.0) {
                    scale[d] += (sdata[v][d]-center[d])*(sdata[v][d]-center[d]);
                }
                if (nval>1) scale[d] /= (nval-1.0);
                else scale[d] = 1.0;
                scale[d] = FastMath.sqrt(scale[d]);
            } else if (statParam.getValue().equals("robust")) {
                double[] samples = new double[(int)nval];
                int n=0;
                for(int v=0; v<surf.getVertexCount(); v++) if (sdata[v][d]!=0.0) {
                    samples[n] = sdata[v][d];
                    n++;
                }
                Percentile measure = new Percentile();
				center[d] = measure.evaluate(samples, 50.0);
				scale[d] = measure.evaluate(samples, 75.0) - measure.evaluate(samples, 25.0);
			}
			System.out.println("Dimension "+d+"("+nval+"; "+statParam.getValue()+"): ");
			System.out.println("center "+center[d]);
			System.out.println("scale "+scale[d]);
			
			for(int v=0; v<surf.getVertexCount(); v++) if (sdata[v][d]!=0.0) {
			    sdata[v][d] = (sdata[v][d]-center[d])/scale[d];
			}
		}
		
		surf.setVertexData(sdata);
		surf = surf.clone();
		
		surf.setName(surf.getName()+"_zscore");
		outputSurface.setValue(surf);
		
	}


}
