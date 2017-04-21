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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
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
public class JistUtilitiesRenameImage extends ProcessingAlgorithm{
	private ParamVolume inputVol;
	private ParamVolume refVol;
	
	private ParamString prefixParam;
	private ParamString suffixParam;
	private	ParamBoolean removeParam;
	private	ParamInteger keptParam;
	
	private ParamVolume outputVol;

	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputVol = new ParamVolume("Image", null,-1,-1,-1,-1));
		inputParams.add(refVol = new ParamVolume("Naming reference"));
		
		inputParams.add(prefixParam = new ParamString("add as prefix"));
		inputParams.add(suffixParam = new ParamString("add as suffix"));
		inputParams.add(removeParam = new ParamBoolean("remove suffixes", false));
		inputParams.add(keptParam = new ParamInteger("kept suffixes", 0,100,1));
		keptParam.setDescription("number of underscore-separated suffixes kept in the reference.") ;
		refVol.setMandatory(false);
		prefixParam.setMandatory(false);
		suffixParam.setMandatory(false);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Utilities");
		inputParams.setLabel("Rename Image");
		inputParams.setName("RenameImage");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Renames an image. No change is made on the data.");
		info.setLongDescription("Renames an image. No change is made on the data. Suffixes separated by _ can be removed.");
		info.setVersion("3.0.5");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	/*
	 * Create output parameters for Image Cropper: cubic volume and extents
	 */
	protected void createOutputParameters(ParamCollection outputParams) {

		outputParams.add(outputVol = new ParamVolume("Renamed Image", null, -1, -1, -1, -1));
		
	}


	/*
	 * Execute Image Cropper
	 */
	protected void execute(CalculationMonitor monitor) {
		
		ImageData vol=inputVol.getImageData();
		ImageData ref;
		if (refVol.getImageData()!=null) ref=refVol.getImageData();
		else ref = vol;
		
		ImageData res = vol.clone();
		
		String newname = ref.getName();
		System.out.println("new name basis: "+newname);
		if (removeParam.getValue().booleanValue()) {
			int idx = newname.indexOf("_");
			for (int n=0;n<keptParam.getValue().intValue();n++) if (idx!=-1) {
				System.out.print(".");
				idx = newname.indexOf("_", idx+1);
			}
			if (idx>-1) newname = newname.substring(0, idx);
		}
		System.out.println("remove all but "+keptParam.getValue().intValue()+"suffixes: "+newname);
		newname = prefixParam.getValue()+newname;
		newname += suffixParam.getValue();
		System.out.println("final name: "+newname);
		
		res.setName(newname);
		outputVol.setValue(res);
	}
}
