package de.mpg.cbs.jist.utilities;

import javax.vecmath.Point3i;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamPointInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;


/*
 * Copy Header
 * @author Pilou Bazin
 *
 */
public class JistUtilitiesCopyHeader extends ProcessingAlgorithm{
	private ParamVolume inputVol;
	private ParamVolume refVol;

	private ParamVolume outputVol;

	private static final String cvsversion = "$Revision: 1.0f $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");
	private static final String shortDescription = "Copy the header from one image into another one. No change is made on the data.";
	private static final String longDescription = "";


	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputVol = new ParamVolume("Image", null,-1,-1,-1,-1));
		inputParams.add(refVol = new ParamVolume("Reference"));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Utilities");
		inputParams.setLabel("Copy Header");
		inputParams.setName("CopyHeader");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription(shortDescription);
		info.setLongDescription(shortDescription + longDescription);
		info.setVersion("3.0");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {

		outputParams.add(outputVol = new ParamVolume("Changed Image", null,-1,-1,-1,-1));
		
	}


	protected void execute(CalculationMonitor monitor) {
		
		ImageData vol=inputVol.getImageData();
		ImageData ref=refVol.getImageData();
		
		ImageData res = vol.clone();
		res.setHeader(ref.getHeader());
		res.setName(res.getName()+"_hdr");
		outputVol.setValue(res);
	}
}
