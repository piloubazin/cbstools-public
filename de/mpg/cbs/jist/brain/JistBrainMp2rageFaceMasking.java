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
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamModel;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;

public class JistBrainMp2rageFaceMasking extends ProcessingAlgorithm
{
    private ParamVolume inv1Image;
    private ParamVolume inv2Image;
    private ParamVolume phs1Image;
    private ParamVolume phs2Image;
    private ParamVolume t1mapImage;
    private ParamVolume isoImage;
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
    private ParamVolume inv1maskImage;
    private ParamVolume inv2maskImage;
    private ParamVolume phs1maskImage;
    private ParamVolume phs2maskImage;
    private ParamVolume t1maskImage;
    private ParamVolume isomaskImage;
    
    protected void createInputParameters(final ParamCollection inputParams) {
        inputParams.add((ParamModel)(this.inv1Image = new ParamVolume("First inversion (Inv1) Image (opt)")));
        inputParams.add((ParamModel)(this.inv2Image = new ParamVolume("Second inversion (Inv2) Image")));
        inputParams.add((ParamModel)(this.phs1Image = new ParamVolume("First inversion phase (Inv1_phs) Image (opt)")));
        inputParams.add((ParamModel)(this.phs2Image = new ParamVolume("Second inversion phase (Inv2_phs) Image (opt)")));
        inputParams.add((ParamModel)(this.t1mapImage = new ParamVolume("T1 Map (T1_Images) Image (opt)")));
        inputParams.add((ParamModel)(this.isoImage = new ParamVolume("T1-weighted (UNI) Image (opt)")));
        this.inv1Image.setMandatory(false);
        this.phs1Image.setMandatory(false);
        this.phs2Image.setMandatory(false);
        this.t1mapImage.setMandatory(false);
        this.isoImage.setMandatory(false);
        inputParams.add((ParamModel)(this.skip0Param = new ParamBoolean("Skip zero values", false)));
        inputParams.setPackage("CBS Tools");
        inputParams.setCategory("Brain Processing.devel");
        inputParams.setLabel("MP2RAGE Face Masking");
        inputParams.setName("Mp2rageFaceMasking");
        final AlgorithmInformation info = this.getAlgorithmInformation();
        info.add(new AlgorithmInformation.AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de", "http://www.cbs.mpg.de/"));
        info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
        info.setDescription("Mask face information for a MP2RAGE dataset. At least an INV2 image is required.");
        info.setVersion("3.0.2");
        info.setStatus(DevelopmentStatus.RC);
        info.setEditable(false);
    }
    
    protected void createOutputParameters(final ParamCollection outputParams) {
        outputParams.add((ParamModel)(this.inv1maskImage = new ParamVolume("Masked Inv1 Image", VoxelType.FLOAT)));
        outputParams.add((ParamModel)(this.inv2maskImage = new ParamVolume("Masked Inv2 Image", VoxelType.FLOAT)));
        outputParams.add((ParamModel)(this.phs1maskImage = new ParamVolume("Masked Inv1_phs Image", VoxelType.FLOAT)));
        outputParams.add((ParamModel)(this.phs2maskImage = new ParamVolume("Masked Inv2_phs Image", VoxelType.FLOAT)));
        outputParams.add((ParamModel)(this.t1maskImage = new ParamVolume("Masked T1 Map Image", VoxelType.FLOAT)));
        outputParams.add((ParamModel)(this.isomaskImage = new ParamVolume("Masked T1-weighted Image", VoxelType.FLOAT)));
        this.t1maskImage.setMandatory(false);
        this.isomaskImage.setMandatory(false);
        this.inv1maskImage.setMandatory(false);
        this.phs1maskImage.setMandatory(false);
        this.phs2maskImage.setMandatory(false);
        outputParams.setName("face masked images");
        outputParams.setLabel("face masked images");
    }
    
    protected void execute(final CalculationMonitor monitor) {
        final int inv2 = 0;
        int nimg = 1;
        int t1map = -1;
        int iso = -1;
        int inv3 = -1;
        int phs1 = -1;
        int phs2 = -1;
        if (this.t1mapImage.getImageData() != null) {
            t1map = nimg;
            ++nimg;
        }
        if (this.isoImage.getImageData() != null) {
            iso = nimg;
            ++nimg;
        }
        if (this.inv1Image.getImageData() != null) {
            inv3 = nimg;
            ++nimg;
        }
        if (this.phs1Image.getImageData() != null) {
            phs1 = nimg;
            ++nimg;
        }
        if (this.phs2Image.getImageData() != null) {
            phs2 = nimg;
            ++nimg;
        }
        ImageDataFloat invImg = new ImageDataFloat(this.inv2Image.getImageData());
        ImageDataFloat t1Img = null;
        ImageDataFloat isoImg = null;
        ImageDataFloat inv1Img = null;
        ImageDataFloat phs1Img = null;
        ImageDataFloat phs2Img = null;
        if (t1map > -1) {
            t1Img = new ImageDataFloat(this.t1mapImage.getImageData());
        }
        if (iso > -1) {
            isoImg = new ImageDataFloat(this.isoImage.getImageData());
        }
        if (inv3 > -1) {
            inv1Img = new ImageDataFloat(this.inv1Image.getImageData());
        }
        if (phs1 > -1) {
            phs1Img = new ImageDataFloat(this.phs1Image.getImageData());
        }
        if (phs2 > -1) {
            phs2Img = new ImageDataFloat(this.phs2Image.getImageData());
        }
        float[][][] buffer = invImg.toArray3d();
        final int nx = invImg.getRows();
        final int ny = invImg.getCols();
        final int nz = invImg.getSlices();
        final int nxyz = nx * ny * nz;
        final float rx = invImg.getHeader().getDimResolutions()[0];
        final float ry = invImg.getHeader().getDimResolutions()[1];
        final float rz = invImg.getHeader().getDimResolutions()[2];
        final ImageHeader.AxisOrientation orx = invImg.getHeader().getAxisOrientation()[0];
        final ImageHeader.AxisOrientation ory = invImg.getHeader().getAxisOrientation()[1];
        final ImageHeader.AxisOrientation orz = invImg.getHeader().getAxisOrientation()[2];
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
        final float[][] image = new float[nimg][nxyz];
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int z = 0; z < nz; ++z) {
                    final int xyz = x + nx * y + nx * ny * z;
                    image[inv2][xyz] = buffer[x][y][z];
                }
            }
        }
        invImg = null;
        buffer = null;
        if (t1map > -1) {
            buffer = t1Img.toArray3d();
            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {
                    for (int z = 0; z < nz; ++z) {
                        final int xyz = x + nx * y + nx * ny * z;
                        image[t1map][xyz] = buffer[x][y][z];
                    }
                }
            }
            t1Img = null;
            buffer = null;
        }
        if (iso > -1) {
            buffer = isoImg.toArray3d();
            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {
                    for (int z = 0; z < nz; ++z) {
                        final int xyz = x + nx * y + nx * ny * z;
                        image[iso][xyz] = buffer[x][y][z];
                    }
                }
            }
            isoImg = null;
            buffer = null;
        }
        if (inv3 > -1) {
            buffer = inv1Img.toArray3d();
            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {
                    for (int z = 0; z < nz; ++z) {
                        final int xyz = x + nx * y + nx * ny * z;
                        image[inv3][xyz] = buffer[x][y][z];
                    }
                }
            }
            inv1Img = null;
            buffer = null;
        }
        if (phs1 > -1) {
            buffer = phs1Img.toArray3d();
            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {
                    for (int z = 0; z < nz; ++z) {
                        final int xyz = x + nx * y + nx * ny * z;
                        image[phs1][xyz] = buffer[x][y][z];
                    }
                }
            }
            phs1Img = null;
            buffer = null;
        }
        if (phs2 > -1) {
            buffer = phs2Img.toArray3d();
            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {
                    for (int z = 0; z < nz; ++z) {
                        final int xyz = x + nx * y + nx * ny * z;
                        image[phs2][xyz] = buffer[x][y][z];
                    }
                }
            }
            phs2Img = null;
            buffer = null;
        }
        final float maxdiff = 1.0E-4f;
        final int itermax = 20;
        float max = -1.0E10f;
        float min = 1.0E10f;
        for (int xyz2 = 0; xyz2 < nxyz; ++xyz2) {
            if (image[inv2][xyz2] > max) {
                max = image[inv2][xyz2];
            }
            if (image[inv2][xyz2] < min) {
                min = image[inv2][xyz2];
            }
        }
        for (int xyz2 = 0; xyz2 < nxyz; ++xyz2) {
            image[inv2][xyz2] = (image[inv2][xyz2] - min) / (max - min);
        }
        System.out.println("image range: " + min + ", " + max);
        final boolean skip0 = this.skip0Param.getValue();
        double mean = 0.0;
        double den = 0.0;
        for (int xyz3 = 0; xyz3 < nxyz; ++xyz3) {
            if (!skip0 || image[inv2][xyz3] != 0.0f) {
                mean += image[inv2][xyz3];
                ++den;
            }
        }
        mean /= den;
        System.out.println("mean parameters: " + mean);
        final double pi = 3.141592653589793;
        final byte model = 102;
        final float[] proba = new float[nxyz];
        for (int xyz4 = 0; xyz4 < nxyz; ++xyz4) {
            if (model == 101) {
                proba[xyz4] = (float)(2.0 / (mean * pi) * Math.exp(-image[inv2][xyz4] * image[inv2][xyz4] / (pi * mean * mean)));
            }
            else if (model == 102) {
                proba[xyz4] = (float)(Math.exp(-image[inv2][xyz4] / mean) / mean);
            }
            proba[xyz4] /= 1.0f + proba[xyz4];
        }
        double diff = 1.0;
        for (int t = 0; t < itermax && diff > maxdiff; ++t) {
            System.out.println("iteration " + (t + 1));
            diff = mean;
            mean = 0.0;
            den = 0.0;
            for (int xyz5 = 0; xyz5 < nxyz; ++xyz5) {
                if (!skip0 || image[inv2][xyz5] != 0.0f) {
                    mean += proba[xyz5] * image[inv2][xyz5];
                    den += proba[xyz5];
                }
            }
            mean /= den;
            System.out.println("mean parameters: " + mean);
            diff = Numerics.abs(diff - mean);
            System.out.println("diff parameters: " + diff);
            for (int xyz5 = 0; xyz5 < nxyz; ++xyz5) {
                if (model == 101) {
                    proba[xyz5] = (float)(2.0 / (mean * pi) * Math.exp(-image[inv2][xyz5] * image[inv2][xyz5] / (pi * mean * mean)));
                }
                else if (model == 102) {
                    proba[xyz5] = (float)(Math.exp(-image[inv2][xyz5] / mean) / mean);
                }
                proba[xyz5] /= 1.0f + proba[xyz5];
            }
        }
        Interface.displayMessage("background-based skull masking");
        final MinMaxFiltering minmax = new MinMaxFiltering(proba, nx, ny, nz, rx, ry, rz);
        float[] brain = minmax.growRegion(new float[] { 0.0f }, new float[] { 0.9f }, new float[] { 0.9f }, 16.0f, 10, false);
        BinaryTopology topo = null;
        Gdm3d gdm = null;
        if (t1map > -1 || iso > -1) {
            float t1max = 0.0f;
            float isomax = 0.0f;
            final float pvmax = 0.0f;
            for (int xyz6 = 0; xyz6 < nxyz; ++xyz6) {
                if (t1map > -1 && image[t1map][xyz6] > t1max) {
                    t1max = image[t1map][xyz6];
                }
                if (iso > -1 && image[iso][xyz6] > isomax) {
                    isomax = image[iso][xyz6];
                }
            }
            final int[] mask = new int[nxyz];
            for (int xyz7 = 0; xyz7 < nxyz; ++xyz7) {
                if (brain[xyz7] > 0.0f) {
                    mask[xyz7] = 1;
                }
                else {
                    mask[xyz7] = 0;
                }
            }
            final float[] balloon = new float[nxyz];
            for (int xyz8 = 0; xyz8 < nxyz; ++xyz8) {
                float force = -0.1f;
                if (t1map > -1) {
                    force = Numerics.max(force, Numerics.bounded((image[t1map][xyz8] / t1max - 0.75f) / 0.25f, -1.0f, 1.0f));
                }
                if (iso > -1) {
                    force = Numerics.max(force, 1.0f - Numerics.bounded((image[iso][xyz8] / isomax - 0.25f) / 0.25f, -1.0f, 1.0f));
                }
                balloon[xyz8] = force;
            }
            topo = new BinaryTopology(mask, nx, ny, nz, rx, ry, rz, "wcs");
            topo.outsideSphericalTopology();
            gdm = new Gdm3d(topo.exportIntSegmentation(), nx, ny, nz, rx, ry, rz, (float[][])null, balloon, 0.0f, 0.1f, 0.9f, "wcs");
            gdm.evolveNarrowBand(100, 0.001f);
            brain = gdm.exportSegmentation();
        }
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
                        buffer[x5][y3][z3] = (max - min) * image[inv2][xyz11] + min;
                    }
                }
            }
        }
        ImageDataFloat bufferData = new ImageDataFloat(buffer);
        bufferData.setHeader(this.inv2Image.getImageData().getHeader());
        bufferData.setName(this.inv2Image.getImageData().getName() + "_msk");
        this.inv2maskImage.setValue((ImageData)bufferData);
        bufferData = null;
        buffer = null;
        if (t1map > -1) {
            buffer = new float[nx][ny][nz];
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (outmask[xyz10]) {
                            buffer[x2][y2][z2] = 0.0f;
                        }
                        else {
                            buffer[x2][y2][z2] = image[t1map][xyz10];
                        }
                    }
                }
            }
            bufferData = new ImageDataFloat(buffer);
            bufferData.setHeader(this.t1mapImage.getImageData().getHeader());
            bufferData.setName(this.t1mapImage.getImageData().getName() + "_msk");
            this.t1maskImage.setValue((ImageData)bufferData);
            bufferData = null;
            buffer = null;
        }
        if (iso > -1) {
            buffer = new float[nx][ny][nz];
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (outmask[xyz10]) {
                            buffer[x2][y2][z2] = 0.0f;
                        }
                        else {
                            buffer[x2][y2][z2] = image[iso][xyz10];
                        }
                    }
                }
            }
            bufferData = new ImageDataFloat(buffer);
            bufferData.setHeader(this.isoImage.getImageData().getHeader());
            bufferData.setName(this.isoImage.getImageData().getName() + "_msk");
            this.isomaskImage.setValue((ImageData)bufferData);
            bufferData = null;
            buffer = null;
        }
        if (inv3 > -1) {
            buffer = new float[nx][ny][nz];
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (outmask[xyz10]) {
                            buffer[x2][y2][z2] = 0.0f;
                        }
                        else {
                            buffer[x2][y2][z2] = image[inv3][xyz10];
                        }
                    }
                }
            }
            bufferData = new ImageDataFloat(buffer);
            bufferData.setHeader(this.inv1Image.getImageData().getHeader());
            bufferData.setName(this.inv1Image.getImageData().getName() + "_msk");
            this.inv1maskImage.setValue((ImageData)bufferData);
            bufferData = null;
            buffer = null;
        }
        if (phs1 > -1) {
            buffer = new float[nx][ny][nz];
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (outmask[xyz10]) {
                            buffer[x2][y2][z2] = 0.0f;
                        }
                        else {
                            buffer[x2][y2][z2] = image[phs1][xyz10];
                        }
                    }
                }
            }
            bufferData = new ImageDataFloat(buffer);
            bufferData.setHeader(this.phs1Image.getImageData().getHeader());
            bufferData.setName(this.phs1Image.getImageData().getName() + "_msk");
            this.phs1maskImage.setValue((ImageData)bufferData);
            bufferData = null;
            buffer = null;
        }
        if (phs2 > -1) {
            buffer = new float[nx][ny][nz];
            for (int x2 = 0; x2 < nx; ++x2) {
                for (int y2 = 0; y2 < ny; ++y2) {
                    for (int z2 = 0; z2 < nz; ++z2) {
                        final int xyz10 = x2 + nx * y2 + nx * ny * z2;
                        if (outmask[xyz10]) {
                            buffer[x2][y2][z2] = 0.0f;
                        }
                        else {
                            buffer[x2][y2][z2] = image[phs2][xyz10];
                        }
                    }
                }
            }
            bufferData = new ImageDataFloat(buffer);
            bufferData.setHeader(this.phs2Image.getImageData().getHeader());
            bufferData.setName(this.phs2Image.getImageData().getName() + "_msk");
            this.phs2maskImage.setValue((ImageData)bufferData);
            bufferData = null;
            buffer = null;
        }
    }
}
