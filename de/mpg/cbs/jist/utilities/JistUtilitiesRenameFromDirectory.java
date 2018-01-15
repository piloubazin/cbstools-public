// 
// Decompiled by Procyon v0.5.30
// 

package de.mpg.cbs.jist.utilities;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipavWrapper;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamModel;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;

public class JistUtilitiesRenameFromDirectory extends ProcessingAlgorithm
{
    private ParamVolume inputVol;
    private ParamInteger pathParam;
    private ParamOption dirParam;
    private ParamVolume outputVol;
    private static final String[] dirTypes;
    private static final String shortDescription = "Renames an image with parts of its path (e.g. its directory). No change is made on the data.";
    private static final String longDescription = "";
    
    protected void createInputParameters(final ParamCollection inputParams) {
        inputParams.add((ParamModel)(this.inputVol = new ParamVolume("Image", (VoxelType)null, -1, -1, -1, -1)));
        inputParams.add((ParamModel)(this.pathParam = new ParamInteger("directory index (1:first, etc.) ", 1, 100, 1)));
        inputParams.add((ParamModel)(this.dirParam = new ParamOption("count index from 1 = ", JistUtilitiesRenameFromDirectory.dirTypes)));
        inputParams.setPackage("CBS Tools");
        inputParams.setCategory("Utilities");
        inputParams.setLabel("Rename From Directory");
        inputParams.setName("RenameFromDirectory");
        final AlgorithmInformation info = this.getAlgorithmInformation();
        info.add(new AlgorithmInformation.AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de", "http://www.cbs.mpg.de/"));
        info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
        info.setDescription("Renames an image with parts of its path (e.g. its directory). No change is made on the data.");
        info.setLongDescription("Renames an image with parts of its path (e.g. its directory). No change is made on the data.");
        info.setVersion("3.0.2");
        info.setEditable(false);
        info.setStatus(DevelopmentStatus.RC);
    }
    
    protected void createOutputParameters(final ParamCollection outputParams) {
        outputParams.add((ParamModel)(this.outputVol = new ParamVolume("Renamed Image", (VoxelType)null, -1, -1, -1, -1)));
    }
    
    protected void execute(final CalculationMonitor monitor) {
        final ImageDataMipavWrapper input = (ImageDataMipavWrapper)this.inputVol.getImageData();
        final String path = input.getModelImageDirect().getImageDirectory();
        System.out.println("directory to parse: " + path);
        final String[] subdirs = path.split("/");
        String prefix = "";
        final int n = this.pathParam.getValue().intValue();
        if (this.dirParam.getValue().equals("root_directory")) {
            prefix = subdirs[n];
        }
        else {
            prefix = subdirs[subdirs.length - n];
        }
        final ImageData res = this.inputVol.getImageData();
        System.out.println("new name prefix: " + prefix);
        String newname = res.getName();
        System.out.println("image name: " + newname);
        newname = prefix + "_" + newname;
        System.out.println("full name: " + newname);
        newname = newname.replaceAll("\\.", "_");
        System.out.println("remove all dots: " + newname);
        res.setName(newname);
        this.outputVol.setValue(res);
    }
    
    static {
        dirTypes = new String[] { "root_directory", "last_directory" };
    }
}
