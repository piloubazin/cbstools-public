/*
 *
 */
package de.mpg.cbs.jist.utilities;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFileCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurfaceCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile.DialogType;

import java.io.File;
import java.util.Vector;
import java.util.Collections;
import java.io.IOException;
import org.apache.commons.io.FileUtils;

import javax.swing.ProgressMonitor;
/*
 * @author Pierre-Louis Bazin (bazin@cbs.mpg.de)
 */
public class JistUtilitiesCopyData extends ProcessingAlgorithm{
	ParamVolumeCollection invol;
	ParamSurfaceCollection insurf;
	ParamFileCollection infile;
	ParamFile outdir;
	
	protected void createInputParameters(ParamCollection inputParams) {
		setRunningInSeparateProcess(true);
		
		inputParams.add(outdir=new ParamFile("Output directory", DialogType.DIRECTORY));
		
		inputParams.add(invol=new ParamVolumeCollection("Data Volumes (opt)",null,-1,-1,-1,-1));
		invol.setMandatory(false);		
		inputParams.add(insurf=new ParamSurfaceCollection("Data Surfaces (opt)"));
		insurf.setMandatory(false);		
		inputParams.add(infile=new ParamFileCollection("Data Files (opt)"));
		infile.setMandatory(false);		
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Utilities");
		inputParams.setLabel("Copy Data");
		inputParams.setName("Copy_Data");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("");
		info.setAffiliation("");
		info.setDescription("Copies the input data to a specified output directory");
		info.setVersion("3.0");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {

	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		// check if directory exists?
		
		// write to directory; if collection is sufficiently big, pass the parameter
		System.out.println("write image data to: "+outdir.getURI());
		/*
		if (invol.getValue()!=null) {
			for (ParamVolume volume : invol.getParamVolumeList()) {
				volume.setLoadAndSaveOnValidate(false);
				volume.getImageData();
				volume.writeAndFreeNow(outdir.getValue());
			}
		}
		
		if (insurf.getValue()!=null) {
			for (ParamFile surface : insurf.getParameters()) {
				surface.setLoadAndSaveOnValidate(false);
				surface.writeAndFreeNow(this);
			}
		}
		*/
		if (invol.getValue()!=null) {
			for (ParamFile file : invol.getParameters()) {
				try {
					FileUtils.copyFileToDirectory(file.getValue(), outdir.getValue());
				} catch (IOException ioe) {
					System.err.println(ioe.getMessage());
				}
			}
		}
		if (insurf.getValue()!=null) {
			for (ParamFile file : insurf.getParameters()) {
				try {
					FileUtils.copyFileToDirectory(file.getValue(), outdir.getValue());
				} catch (IOException ioe) {
					System.err.println(ioe.getMessage());
				}
			}
		}
		if (infile.getValue()!=null) {
			for (ParamFile file : infile.getParameters()) {
				try {
					FileUtils.copyFileToDirectory(file.getValue(), outdir.getValue());
				} catch (IOException ioe) {
					System.err.println(ioe.getMessage());
				}
			}
		}
	}
	
	

}
