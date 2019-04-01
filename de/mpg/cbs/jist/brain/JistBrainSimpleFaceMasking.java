// 
// Decompiled by Procyon v0.5.30
// 

package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import de.mpg.cbs.libraries.Morphology;
import de.mpg.cbs.libraries.Gdm3d;
import de.mpg.cbs.libraries.BinaryTopology;
import de.mpg.cbs.methods.MinMaxFiltering;
import de.mpg.cbs.utilities.Interface;
import de.mpg.cbs.utilities.Numerics;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamModel;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;

public class JistBrainSimpleFaceMasking extends ProcessingAlgorithm
{
    private ParamVolume inputImage;
    private ParamVolume stripMaskImage;
    private static final byte HNORM = 101;
    private static final byte EXP = 102;
    private static final byte UNKNOWN = 0;
    private static final byte XPOS = 1;
    private static final byte XNEG = -1;
    private static final byte YPOS = 2;
    private static final byte YNEG = -2;
    private static final byte ZPOS = 3;
    private static final byte ZNEG = -3;
    private ParamBoolean skip0Param;
    private ParamVolume maskedInputImage;
    private ParamVolume faceMaskImage;
    
    protected void createInputParameters(final ParamCollection inputParams) {
        inputParams.add((ParamModel)(this.inputImage = new ParamVolume("Input Intensity Image")));
        inputParams.add((ParamModel)(this.stripMaskImage = new ParamVolume("Stripping Mask Image")));
        inputParams.add((ParamModel)(this.skip0Param = new ParamBoolean("Skip zero values", false)));
        inputParams.setPackage("CBS Tools");
        inputParams.setCategory("Brain Processing.devel");
        inputParams.setLabel("Simple Face Masking");
        inputParams.setName("SimpleFaceMasking");
        final AlgorithmInformation info = this.getAlgorithmInformation();
        info.add(new AlgorithmInformation.AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de", "http://www.cbs.mpg.de/"));
        info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
        info.setDescription("Mask face information for a dataset, given a skull stripping mask.");
        info.setVersion("3.1.4");
        info.setStatus(DevelopmentStatus.RC);
        info.setEditable(false);
    }
    
    protected void createOutputParameters(final ParamCollection outputParams) {
        outputParams.add((ParamModel)(this.maskedInputImage = new ParamVolume("Masked Input Image", VoxelType.FLOAT)));
        outputParams.add((ParamModel)(this.faceMaskImage = new ParamVolume("Face Mask Image", VoxelType.UBYTE)));
        
        outputParams.setName("face masked images");
        outputParams.setLabel("face masked images");
    }
    
    protected void execute(final CalculationMonitor monitor) {

        ImageDataFloat inputImg = new ImageDataFloat(this.inputImage.getImageData());
        ImageDataUByte maskImg = new ImageDataUByte(this.stripMaskImage.getImageData());
        
        final int nx = inputImg.getRows();
        final int ny = inputImg.getCols();
        final int nz = inputImg.getSlices();
        final int nxyz = nx * ny * nz;
        final float rx = inputImg.getHeader().getDimResolutions()[0];
        final float ry = inputImg.getHeader().getDimResolutions()[1];
        final float rz = inputImg.getHeader().getDimResolutions()[2];
        final ImageHeader.AxisOrientation orx = inputImg.getHeader().getAxisOrientation()[0];
        final ImageHeader.AxisOrientation ory = inputImg.getHeader().getAxisOrientation()[1];
        final ImageHeader.AxisOrientation orz = inputImg.getHeader().getAxisOrientation()[2];
        int facedir = 0;
        if (orx == ImageHeader.AxisOrientation.A2P_TYPE) {
            facedir = 1;
        }
        else if (orx == ImageHeader.AxisOrientation.P2A_TYPE) {
            facedir = -1;
        }
        else if (ory == ImageHeader.AxisOrientation.A2P_TYPE) {
            facedir = 2;
        }
        else if (ory == ImageHeader.AxisOrientation.P2A_TYPE) {
            facedir = -2;
        }
        else if (orz == ImageHeader.AxisOrientation.A2P_TYPE) {
            facedir = 3;
        }
        else if (orz == ImageHeader.AxisOrientation.P2A_TYPE) {
            facedir = -3;
        }
        int bodydir = 0;
        if (orx == ImageHeader.AxisOrientation.I2S_TYPE) {
            bodydir = 1;
        }
        else if (orx == ImageHeader.AxisOrientation.S2I_TYPE) {
            bodydir = -1;
        }
        else if (ory == ImageHeader.AxisOrientation.I2S_TYPE) {
            bodydir = 2;
        }
        else if (ory == ImageHeader.AxisOrientation.S2I_TYPE) {
            bodydir = -2;
        }
        else if (orz == ImageHeader.AxisOrientation.I2S_TYPE) {
            bodydir = 3;
        }
        else if (orz == ImageHeader.AxisOrientation.S2I_TYPE) {
            bodydir = -3;
        }
        System.out.println("face direction: " + facedir);
        System.out.println("body direction: " + bodydir);
        final float[] image = new float[nxyz];
        float[][][] buffer = inputImg.toArray3d();
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int z = 0; z < nz; ++z) {
                    final int xyz = x + nx * y + nx * ny * z;
                    image[xyz] = buffer[x][y][z];
                }
            }
        }
        inputImg = null;
        buffer = null;
 
        final byte[] brain = new byte[nxyz];
        byte[][][] buffer2 = maskImg.toArray3d();
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int z = 0; z < nz; ++z) {
                    final int xyz = x + nx * y + nx * ny * z;
                    brain[xyz] = buffer2[x][y][z];
                }
            }
        }
        maskImg = null;
        buffer2 = null;
        
        System.out.println("create outside mask");
        boolean[] outmask = new boolean[nxyz];
        for (int xyz9 = 0; xyz9 < nxyz; ++xyz9) {
            outmask[xyz9] = (brain[xyz9] <= 0.0f);
        }
        final float erosion = 5.0f;
        outmask = Morphology.erodeObject(outmask, nx, ny, nz, Numerics.ceil(erosion / rx), Numerics.ceil(erosion / ry), Numerics.ceil(erosion / rz));
        System.out.println("limit to face direction");
        final boolean[] dirmask = new boolean[nxyz];
        final float front = 20.0f;
        if (facedir == 1) {
            int xmin = nx;
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (!outmask[xyz10] && x2 < xmin) {
                            xmin = x2;
                        }
                    }
                }
            }
            xmin = Numerics.floor(xmin + front / rx);
            for (int y3 = 0; y3 < ny; ++y3) {
                for (int z3 = 0; z3 < nz; ++z3) {
                    for (int x3 = nx - 1; x3 > 0; --x3) {
                        final int xyz10 = x3 + nx * y3 + nx * ny * z3;
                        if (!outmask[xyz10]) {
                            for (int xi = 0; xi < Numerics.min(x3, xmin); ++xi) {
                                dirmask[xyz10 - x3 + xi] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (facedir == -1) {
            int xmax = 0;
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (!outmask[xyz10] && x2 > xmax) {
                            xmax = x2;
                        }
                    }
                }
            }
            xmax = Numerics.ceil(xmax - front / rx);
            for (int y3 = 0; y3 < ny; ++y3) {
                for (int z3 = 0; z3 < nz; ++z3) {
                    for (int x3 = 0; x3 < nx; ++x3) {
                        final int xyz10 = x3 + nx * y3 + nx * ny * z3;
                        if (!outmask[xyz10]) {
                            for (int xi = nx - 1; xi > Numerics.max(x3, xmax); --xi) {
                                dirmask[xyz10 - x3 + xi] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (facedir == 2) {
            int ymin = ny;
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (!outmask[xyz10] && y2 < ymin) {
                            ymin = y2;
                        }
                    }
                }
            }
            ymin = Numerics.floor(ymin + front / ry);
            for (int z4 = 0; z4 < nz; ++z4) {
                for (int x4 = 0; x4 < nx; ++x4) {
                    for (int y4 = ny - 1; y4 > 0; --y4) {
                        final int xyz10 = x4 + nx * y4 + nx * ny * z4;
                        if (!outmask[xyz10]) {
                            for (int yi = 0; yi < Numerics.min(y4, ymin); ++yi) {
                                dirmask[xyz10 - nx * y4 + nx * yi] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (facedir == -2) {
            int ymax = 0;
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (!outmask[xyz10] && y2 > ymax) {
                            ymax = y2;
                        }
                    }
                }
            }
            ymax = Numerics.ceil(ymax - front / ry);
            for (int z4 = 0; z4 < nz; ++z4) {
                for (int x4 = 0; x4 < nx; ++x4) {
                    for (int y4 = 0; y4 < ny; ++y4) {
                        final int xyz10 = x4 + nx * y4 + nx * ny * z4;
                        if (!outmask[xyz10]) {
                            for (int yi = ny - 1; yi > Numerics.max(y4, ymax); --yi) {
                                dirmask[xyz10 - nx * y4 + nx * yi] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (facedir == 3) {
            int zmin = nz;
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (!outmask[xyz10] && z2 < zmin) {
                            zmin = z2;
                        }
                    }
                }
            }
            zmin = Numerics.floor(zmin + front / rz);
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = nz - 1; z2 > 0; --z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (!outmask[xyz10]) {
                            for (int zi = 0; zi < Numerics.min(z2, zmin); ++zi) {
                                dirmask[xyz10 - nx * ny * z2 + nx * ny * zi] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (facedir == -3) {
            int zmax = 0;
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (!outmask[xyz10] && z2 > zmax) {
                            zmax = z2;
                        }
                    }
                }
            }
            zmax = Numerics.floor(zmax - front / rz);
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (!outmask[xyz10]) {
                            for (int zi = nz - 1; zi > Numerics.max(z2, zmax); --zi) {
                                dirmask[xyz10 - nx * ny * z2 + nx * ny * zi] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else {
            for (int x5 = 0; x5 < nx; ++x5) {
                for (int y3 = 0; y3 < ny; ++y3) {
                    for (int z3 = 0; z3 < nz; ++z3) {
                        final int xyz11 = x5 + nx * y3 + nx * ny * z3;
                        dirmask[xyz11] = true;
                    }
                }
            }
        }
        for (int x5 = 0; x5 < nx; ++x5) {
            for (int y3 = 0; y3 < ny; ++y3) {
                for (int z3 = 0; z3 < nz; ++z3) {
                    final int xyz11 = x5 + nx * y3 + nx * ny * z3;
                    dirmask[xyz11] = (outmask[xyz11] && dirmask[xyz11]);
                }
            }
        }
        if (bodydir == 1) {
            for (int y5 = 0; y5 < ny; ++y5) {
                for (int z4 = 0; z4 < nz; ++z4) {
                    for (int x4 = nx - 1; x4 > 0; --x4) {
                        final int xyz11 = x4 + nx * y5 + nx * ny * z4;
                        if (dirmask[xyz11]) {
                            for (int xi2 = 0; xi2 < x4 && outmask[xyz11 - x4 + xi2]; ++xi2) {
                                dirmask[xyz11 - x4 + xi2] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (bodydir == -1) {
            for (int y5 = 0; y5 < ny; ++y5) {
                for (int z4 = 0; z4 < nz; ++z4) {
                    for (int x4 = 0; x4 < nx; ++x4) {
                        final int xyz11 = x4 + nx * y5 + nx * ny * z4;
                        if (dirmask[xyz11]) {
                            for (int xi2 = nx - 1; xi2 > x4 && outmask[xyz11 - x4 + xi2]; --xi2) {
                                dirmask[xyz11 - x4 + xi2] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (bodydir == 2) {
            for (int z5 = 0; z5 < nz; ++z5) {
                for (int x2 = 0; x2 < nx; ++x2) {
                    for (int y2 = ny - 1; y2 > 0; --y2) {
                        final int xyz11 = x2 + nx * y2 + nx * ny * z5;
                        if (dirmask[xyz11]) {
                            for (int yi2 = 0; yi2 < y2 && outmask[xyz11 - nx * y2 + nx * yi2]; ++yi2) {
                                dirmask[xyz11 - nx * y2 + nx * yi2] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (bodydir == -2) {
            for (int z5 = 0; z5 < nz; ++z5) {
                for (int x2 = 0; x2 < nx; ++x2) {
                    for (int y2 = 0; y2 < ny; ++y2) {
                        final int xyz11 = x2 + nx * y2 + nx * ny * z5;
                        if (dirmask[xyz11]) {
                            for (int yi2 = ny - 1; yi2 > y2 && outmask[xyz11 - nx * y2 + nx * yi2]; --yi2) {
                                dirmask[xyz11 - nx * y2 + nx * yi2] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (bodydir == 3) {
            for (int x5 = 0; x5 < nx; ++x5) {
                for (int y3 = 0; y3 < ny; ++y3) {
                    for (int z3 = nz - 1; z3 > 0; --z3) {
                        final int xyz11 = x5 + nx * y3 + nx * ny * z3;
                        if (dirmask[xyz11]) {
                            for (int zi2 = 0; zi2 < z3 && outmask[xyz11 - nx * ny * z3 + nx * ny * zi2]; ++zi2) {
                                dirmask[xyz11 - nx * ny * z3 + nx * ny * zi2] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (bodydir == -3) {
            for (int x5 = 0; x5 < nx; ++x5) {
                for (int y3 = 0; y3 < ny; ++y3) {
                    for (int z3 = 0; z3 < nz; ++z3) {
                        final int xyz11 = x5 + nx * y3 + nx * ny * z3;
                        if (dirmask[xyz11]) {
                            for (int zi2 = nz - 1; zi2 > z3 && outmask[xyz11 - nx * ny * z3 + nx * ny * zi2]; --zi2) {
                                dirmask[xyz11 - nx * ny * z3 + nx * ny * zi2] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }
        outmask = dirmask;
        System.out.println("mask original images");
        buffer = new float[nx][ny][nz];
        for (int x5 = 0; x5 < nx; ++x5) {
            for (int y3 = 0; y3 < ny; ++y3) {
                for (int z3 = 0; z3 < nz; ++z3) {
                    final int xyz11 = x5 + nx * y3 + nx * ny * z3;
                    if (outmask[xyz11]) {
                        buffer[x5][y3][z3] = 0.0f;
                    }
                    else {
                        buffer[x5][y3][z3] = image[xyz11];
                    }
                }
            }
        }
        ImageDataFloat bufferData = new ImageDataFloat(buffer);
        bufferData.setHeader(this.inputImage.getImageData().getHeader());
        bufferData.setName(this.inputImage.getImageData().getName() + "_msk");
        this.maskedInputImage.setValue((ImageData)bufferData);
        bufferData = null;
        buffer = null;

        buffer2 = new byte[nx][ny][nz];
        for (int x2 = 0; x2 < nx; ++x2) {
            for (int y2 = 0; y2 < ny; ++y2) {
                for (int z2 = 0; z2 < nz; ++z2) {
                    final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                    if (outmask[xyz10]) {
                        buffer2[x2][y2][z2] = 0;
                    }
                    else {
                        buffer2[x2][y2][z2] = 1;
                    }
                }
            }
        }
        ImageDataUByte buffer2Data = new ImageDataUByte(buffer2);
        buffer2Data.setHeader(this.inputImage.getImageData().getHeader());
        buffer2Data.setName(this.inputImage.getImageData().getName() + "_fmsk");
        this.faceMaskImage.setValue((ImageData)buffer2Data);
        buffer2Data = null;
        buffer2 = null;
    }
}
