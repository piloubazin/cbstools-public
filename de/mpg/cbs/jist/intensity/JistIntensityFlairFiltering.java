// 
// Decompiled by Procyon v0.5.30
// 

package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import de.mpg.cbs.utilities.Numerics;
import de.mpg.cbs.utilities.Interface;
import de.mpg.cbs.libraries.ImageStatistics;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamModel;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;

public class JistIntensityFlairFiltering extends ProcessingAlgorithm
{
    private ParamVolume inputImage;
    private ParamOption algoParam;
    private ParamDouble inratioParam;
    private ParamDouble outratioParam;
    private ParamVolume priorImage;
    private ParamInteger priorParam;
    private static final String[] algotypes;
    private String algotype;
    private ParamVolume probaImage;
    private ParamVolume filterImage;
    
    public JistIntensityFlairFiltering() {
        this.algotype = "prior-based-normalization";
    }
    
    protected void createInputParameters(final ParamCollection inputParams) {
        inputParams.add((ParamModel)(this.inputImage = new ParamVolume("Input Image")));
        inputParams.add((ParamModel)(this.algoParam = new ParamOption("Algorithm", JistIntensityFlairFiltering.algotypes)));
        this.algoParam.setValue(this.algotype);
        inputParams.add((ParamModel)(this.inratioParam = new ParamDouble("Inlier ratio", 0.0, 1.0, 0.01)));
        inputParams.add((ParamModel)(this.outratioParam = new ParamDouble("Outlier ratio", 0.0, 1.0, 0.001)));
        inputParams.add((ParamModel)(this.priorImage = new ParamVolume("Prior Segmentation")));
        this.priorImage.setMandatory(false);
        inputParams.add((ParamModel)(this.priorParam = new ParamInteger("Lesion label", 0, 1000, 10)));
        inputParams.setPackage("CBS Tools");
        inputParams.setCategory("Intensity.FLAIR");
        inputParams.setLabel("Flair Filtering");
        inputParams.setName("FlairFiltering");
        final AlgorithmInformation info = this.getAlgorithmInformation();
        info.add(new AlgorithmInformation.AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de", "http://www.cbs.mpg.de/"));
        info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
        info.setDescription("Filter a FLAIR image to highlight hyperintensities (either based on given parameters or a prior segmentation)");
        info.setVersion("3.0.1");
        info.setStatus(DevelopmentStatus.RC);
        info.setEditable(false);
    }
    
    protected void createOutputParameters(final ParamCollection outputParams) {
        outputParams.add((ParamModel)(this.probaImage = new ParamVolume("Proba Image", VoxelType.FLOAT)));
        outputParams.add((ParamModel)(this.filterImage = new ParamVolume("Filtered Image", VoxelType.FLOAT)));
        outputParams.setName("flair filter images");
        outputParams.setLabel("flair filter images");
    }
    
    protected void execute(final CalculationMonitor monitor) {
        final ImageDataFloat inImg = new ImageDataFloat(this.inputImage.getImageData());
        final float[][][] image = inImg.toArray3d();
        final int nx = inImg.getRows();
        final int ny = inImg.getCols();
        final int nz = inImg.getSlices();
        final float rx = inImg.getHeader().getDimResolutions()[0];
        final float ry = inImg.getHeader().getDimResolutions()[1];
        final float rz = inImg.getHeader().getDimResolutions()[2];
        final boolean[][][] mask = new boolean[nx][ny][nz];
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int z = 0; z < nz; ++z) {
                    mask[x][y][z] = (image[x][y][z] > 0.0f);
                }
            }
        }
        final float[][][] proba = new float[nx][ny][nz];
        final float[][][] filter = new float[nx][ny][nz];
        if (this.algoParam.getValue().equals("inlier-outlier-proba")) {
            final float Iin = ImageStatistics.robustMaximum(image, mask, this.inratioParam.getValue().floatValue(), 4, nx, ny, nz);
            final float Iout = ImageStatistics.robustMaximum(image, mask, this.outratioParam.getValue().floatValue(), 4, nx, ny, nz);
            Interface.displayMessage("image inlier / outliers: " + Iin + ", " + Iout + "\n");
            final float mean = Iout - Iin;
            final float offset = Iout;
            float avg = 0.0f;
            float sum = 0.0f;
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        proba[x2][y2][z2] = (float)Math.exp(-0.5f * Numerics.square(Numerics.max(offset - image[x2][y2][z2], 0.0f) / mean));
                        proba[x2][y2][z2] = Numerics.square(proba[x2][y2][z2]);
                        avg += proba[x2][y2][z2] * image[x2][y2][z2];
                        sum += proba[x2][y2][z2];
                    }
                }
            }
            avg /= sum;
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        filter[x2][y2][z2] = Numerics.max(2.0f * proba[x2][y2][z2] - 1.0f, 0.0f) * (offset + mean) + (1.0f - Numerics.max(2.0f * proba[x2][y2][z2] - 1.0f, 0.0f)) * image[x2][y2][z2];
                    }
                }
            }
        }
        else if (this.algoParam.getValue().equals("prior-based-normalization")) {
            final ImageDataInt priorImg = new ImageDataInt(this.priorImage.getImageData());
            final int[][][] prior = priorImg.toArray3d();
            final int lesionlb = this.priorParam.getValue().intValue();
            int braincount = 0;
            int lesioncount = 0;
            final int bglb = prior[0][0][0];
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        if (prior[x2][y2][z2] > bglb) {
                            ++braincount;
                        }
                        if (prior[x2][y2][z2] == lesionlb) {
                            ++lesioncount;
                        }
                    }
                }
            }
            final float ratio = lesioncount / braincount;
            Interface.displayMessage("lesion ratio: " + ratio + "\n");
            final float Iout2 = ImageStatistics.robustMaximum(image, mask, ratio, 4, nx, ny, nz);
            Interface.displayMessage("image outliers: " + Iout2 + "\n");
            for (int x3 = 0; x3 < nx; ++x3) {
                for (int y3 = 0; y3 < ny; ++y3) {
                    for (int z3 = 0; z3 < nz; ++z3) {
                        if (image[x3][y3][z3] >= Iout2) {
                            proba[x3][y3][z3] = 1.0f;
                        }
                        else {
                            proba[x3][y3][z3] = 0.0f;
                        }
                        filter[x3][y3][z3] = Numerics.bounded(image[x3][y3][z3], 0.0f, Iout2);
                    }
                }
            }
        }
        ImageDataFloat bufferData = new ImageDataFloat(proba);
        bufferData.setHeader(this.inputImage.getImageData().getHeader());
        bufferData.setName(this.inputImage.getImageData().getName() + "_proba");
        this.probaImage.setValue((ImageData)bufferData);
        bufferData = null;
        bufferData = new ImageDataFloat(filter);
        bufferData.setHeader(this.inputImage.getImageData().getHeader());
        bufferData.setName(this.inputImage.getImageData().getName() + "_filter");
        this.filterImage.setValue((ImageData)bufferData);
        bufferData = null;
    }
    
    static {
        algotypes = new String[] { "inlier-outlier-proba", "prior-based-normalization" };
    }
}
