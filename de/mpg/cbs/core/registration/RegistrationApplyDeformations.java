package de.mpg.cbs.core.registration;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class RegistrationApplyDeformations {

	// jist containers
	private float[] sourceImage;
	private float[] referenceImage = null;
	private float[] deformation1Image;
	private float[] deformation2Image = null;
	private float[] deformation3Image = null;
	private float[] deformation4Image = null;
	private String type1Option;
	private String type2Option = "none";
	private String type3Option = "none";
	private String type4Option = "none";
	private String interpOption = "nearest";
	private String padOption = "closest";
	
	public static final String[] types = {"none", "deformation(voxels)", "mapping(voxels)", "deformation(mm)", "mapping(mm)"};
	public static final String[] interp = {"nearest", "linear"};
	public static final String[] pads = {"closest", "zero", "min", "max"};
	
	private float[] deformedImage;
	
	private int nsx, nsy, nsz, nst, nsxyz;
	private float rsx, rsy, rsz;
	private int nrx, nry, nrz, nrt, nrxyz;
	private float rrx, rry, rrz;
	
	private int nd1x, nd1y, nd1z, nd1xyz;
	private float rd1x, rd1y, rd1z;
	private int nd2x, nd2y, nd2z, nd2xyz;
	private float rd2x, rd2y, rd2z;
	private int nd3x, nd3y, nd3z, nd3xyz;
	private float rd3x, rd3y, rd3z;
	private int nd4x, nd4y, nd4z, nd4xyz;
	private float rd4x, rd4y, rd4z;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;
	private static final byte T = 3;

	// create inputs
	public final void setImageToDeform(float[] val) { sourceImage = val; }
	//public final void setReferenceImage(float[] val) { referenceImage = val; }
	public final void setDeformationMapping1(float[] val) { deformation1Image = val; }
	public final void setDeformationType1(String val) { type1Option = val; }
	public final void setDeformationMapping2(float[] val) { deformation2Image = val; }
	public final void setDeformationType2(String val) { type2Option = val; }
	public final void setDeformationMapping3(float[] val) { deformation3Image = val; }
	public final void setDeformationType3(String val) { type3Option = val; }
	public final void setDeformationMapping4(float[] val) { deformation4Image = val; }
	public final void setDeformationType4(String val) { type4Option = val; }
	public final void setInterpolationType(String val) { interpOption = val; }
	public final void setImagePadding(String val) { padOption = val; }
		
	
	public final void setImageDimensions(int x, int y, int z) { nsx=x; nsy=y; nsz=z; nst=1; nsxyz=nsx*nsy*nsz; }
	public final void setImageDimensions(int x, int y, int z, int t) { nsx=x; nsy=y; nsz=z; nst=t; nsxyz=nsx*nsy*nsz; }
	public final void setImageDimensions(int[] dim) { nsx=dim[0]; nsy=dim[1]; nsz=dim[2]; nst=dim[3]; nsxyz=nsx*nsy*nsz; }
	
	public final void setImageResolutions(float x, float y, float z) { rsx=x; rsy=y; rsz=z; }
	public final void setImageResolutions(float[] res) { rsx=res[0]; rsy=res[1]; rsz=res[2]; }

	//public final void setReferenceDimensions(int x, int y, int z) { nrx=x; nry=y; nrz=z; nrt=t; nrxyz=nrx*nry*nrz; }
	//public final void setReferenceDimensions(int[] dim) { nrx=dim[0]; nry=dim[1]; nrz=dim[2]; nrt=dim[3]; nrxyz=nrx*nry*nrz; }
	
	//public final void setReferenceResolutions(float x, float y, float z) { rrx=x; rry=y; rrz=z; }
	//public final void setReferenceResolutions(float[] res) { rrx=res[0]; rry=res[1]; rrz=res[2]; }

	public final void setDeformation1Dimensions(int x, int y, int z) { nd1x=x; nd1y=y; nd1z=z; nd1xyz=nd1x*nd1y*nd1z; }
	public final void setDeformation1Dimensions(int[] dim) { nd1x=dim[0]; nd1y=dim[1]; nd1z=dim[2]; nd1xyz=nd1x*nd1y*nd1z; }
	
	public final void setDeformation1Resolutions(float x, float y, float z) { rd1x=x; rd1y=y; rd1z=z; }
	public final void setDeformation1Resolutions(float[] res) { rd1x=res[0]; rd1y=res[1]; rd1z=res[2]; }

	public final void setDeformation2Dimensions(int x, int y, int z) { nd2x=x; nd2y=y; nd2z=z; nd2xyz=nd2x*nd2y*nd2z; }
	public final void setDeformation2Dimensions(int[] dim) { nd2x=dim[0]; nd2y=dim[1]; nd2z=dim[2]; nd2xyz=nd2x*nd2y*nd2z; }
	
	public final void setDeformation2Resolutions(float x, float y, float z) { rd2x=x; rd2y=y; rd2z=z; }
	public final void setDeformation2Resolutions(float[] res) { rd2x=res[0]; rd2y=res[1]; rd2z=res[2]; }

	public final void setDeformation3Dimensions(int x, int y, int z) { nd3x=x; nd3y=y; nd3z=z; nd3xyz=nd3x*nd3y*nd3z; }
	public final void setDeformation3Dimensions(int[] dim) { nd3x=dim[0]; nd3y=dim[1]; nd3z=dim[2]; nd3xyz=nd3x*nd3y*nd3z; }
	
	public final void setDeformation3Resolutions(float x, float y, float z) { rd3x=x; rd3y=y; rd3z=z; }
	public final void setDeformation3Resolutions(float[] res) { rd3x=res[0]; rd3y=res[1]; rd3z=res[2]; }

	public final void setDeformation4Dimensions(int x, int y, int z) { nd4x=x; nd4y=y; nd4z=z; nd4xyz=nd4x*nd4y*nd4z; }
	public final void setDeformation4Dimensions(int[] dim) { nd4x=dim[0]; nd4y=dim[1]; nd4z=dim[2]; nd4xyz=nd4x*nd4y*nd4z; }
	
	public final void setDeformation4Resolutions(float x, float y, float z) { rd4x=x; rd4y=y; rd4z=z; }
	public final void setDeformation4Resolutions(float[] res) { rd4x=res[0]; rd4y=res[1]; rd4z=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Registration"; }
	public final String getLabel() { return "Apply Deformations"; }
	public final String getName() { return "ApplyDeformations"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Apply multiple deformations from non-linear transformations. Transformations can be a deformation field or a coordinate mapping."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.3"; };
	
	public final float[] getDeformedImage() { return deformedImage; }
	
	public void execute() {
				
        // deformation: in reference space
        System.out.println("load deformation 1");
        // scale to voxels if needed
        if (type1Option.endsWith("(mm)")) {
            System.out.println("normalize to resolution ("+rd1x+", "+rd1y+", "+rd1z+")");
            for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) for (int z=0;z<nd1z;z++) {
                int xyz = x + nd1x*y + nd1x*nd1y*z;
                deformation1Image[xyz + X*nd1xyz] /= rd1x; 	
                deformation1Image[xyz + Y*nd1xyz] /= rd1y; 	
                deformation1Image[xyz + Z*nd1xyz] /= rd1z; 	
            }
        }
        // turn into a mapping if needed
        if (type1Option.startsWith("deformation")) {
            for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) for (int z=0;z<nd1z;z++) {
                int xyz = x + nd1x*y + nd1x*nd1y*z;
                deformation1Image[xyz + X*nd1xyz] += x;
                deformation1Image[xyz + Y*nd1xyz] += y;
                deformation1Image[xyz + Z*nd1xyz] += z;
            }
        }
        // check for bad borders
        boolean[] boundary = new boolean[nd1x*nd1y*nd1z];
        boolean growBoundaries = false;
        for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) for (int z=0;z<nd1z;z++) {
            int xyz = x + nd1x*y + nd1x*nd1y*z;
            if (deformation1Image[xyz + X*nd1xyz]==0 && deformation1Image[xyz + Y*nd1xyz]==0 && deformation1Image[xyz + Z*nd1xyz]==0) {
                for (byte k=0;k<6;k++) {
                    if (x+Ngb.x[k]>=0 && x+Ngb.x[k]<nd1x && y+Ngb.y[k]>=0 && y+Ngb.y[k]<nd1y && z+Ngb.z[k]>=0 && z+Ngb.z[k]<nd1z) {
                        int ngb = Ngb.neighborIndex(k, xyz,nd1x,nd1y,nd1z);
                        if (deformation1Image[ngb + X*nd1xyz]!=0 || deformation1Image[ngb + Y*nd1xyz]!=0 || deformation1Image[ngb + Z*nd1xyz]!=0) {
                            growBoundaries = true;
                            boundary[ngb] = true;
                            k=6;
                        }
                    }
                }
            }
        }
        while (growBoundaries) {
            boolean[] changed = new boolean[nd1x*nd1y*nd1z];
            growBoundaries = false;
            for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) for (int z=0;z<nd1z;z++) {
                int xyz = x + nd1x*y + nd1x*nd1y*z;
                if (boundary[xyz]) {
                    for (byte k=0;k<6;k++) {
                        if (x+Ngb.x[k]>=0 && x+Ngb.x[k]<nd1x && y+Ngb.y[k]>=0 && y+Ngb.y[k]<nd1y && z+Ngb.z[k]>=0 && z+Ngb.z[k]<nd1z) {
                            int ngb = Ngb.neighborIndex(k, xyz,nd1x,nd1y,nd1z);
                            if (deformation1Image[ngb + X*nd1xyz]==0 && deformation1Image[ngb + Y*nd1xyz]==0 && deformation1Image[ngb + Z*nd1xyz]==0) {
                                deformation1Image[ngb + X*nd1xyz] = deformation1Image[xyz + X*nd1xyz];
                                deformation1Image[ngb + Y*nd1xyz] = deformation1Image[xyz + Y*nd1xyz];
                                deformation1Image[ngb + Z*nd1xyz] = deformation1Image[xyz + Z*nd1xyz];
                                growBoundaries = true;
                                changed[ngb] = true;
                            }
                        }
                    }
                }
            }
            boundary = changed;
        }
        
        // add to final deformation
        float[] deformation = deformation1Image;
        
        nrx = nd1x; nry = nd1y; nrz = nd1z;
        rrx = rd1x; rry = rd1y; rrz = rd1z;
        
        // second deformation if given
        if (!type2Option.equals("none") && deformation2Image!=null) {
            System.out.println("load deformation 2");
            // scale to voxels if needed
            if (type2Option.endsWith("(mm)")) {
                System.out.println("normalize to resolution ("+rd2x+", "+rd2y+", "+rd2z+")");
                for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
                    int xyz = x + nd2x*y + nd2x*nd2y*z;
                    deformation2Image[xyz + X*nd2xyz] /= rd2x; 	
                    deformation2Image[xyz + Y*nd2xyz] /= rd2y; 	
                    deformation2Image[xyz + Z*nd2xyz] /= rd2z; 	
                }
            }
            // turn into a mapping if needed
            if (type2Option.startsWith("deformation")) {
                for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
                    int xyz = x + nd2x*y + nd2x*nd2y*z;
                    deformation2Image[xyz + X*nd2xyz] += x;
                    deformation2Image[xyz + Y*nd2xyz] += y;
                    deformation2Image[xyz + Z*nd2xyz] += z;
                }
            }
            // check for bad borders
            boundary = new boolean[nd2x*nd2y*nd2z];
            growBoundaries = false;
            for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
                int xyz = x + nd2x*y + nd2x*nd2y*z;
                if (deformation2Image[xyz + X*nd2xyz]==0 && deformation2Image[xyz + Y*nd2xyz]==0 && deformation2Image[xyz + Z*nd2xyz]==0) {
                    for (byte k=0;k<6;k++) {
                        if (x+Ngb.x[k]>=0 && x+Ngb.x[k]<nd2x && y+Ngb.y[k]>=0 && y+Ngb.y[k]<nd2y && z+Ngb.z[k]>=0 && z+Ngb.z[k]<nd2z) {
                            int ngb = Ngb.neighborIndex(k, xyz,nd2x,nd2y,nd2z);
                            if (deformation2Image[ngb + X*nd2xyz]!=0 || deformation2Image[ngb + Y*nd2xyz]!=0 || deformation2Image[ngb + Z*nd2xyz]!=0) {
                                growBoundaries = true;
                                boundary[ngb] = true;
                                k=6;
                            }
                        }
                    }
                }
            }
            while (growBoundaries) {
                boolean[] changed = new boolean[nd2x*nd2y*nd2z];
                growBoundaries = false;
                for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
                    int xyz = x + nd2x*y + nd2x*nd2y*z;
                    if (boundary[xyz]) {
                        for (byte k=0;k<6;k++) {
                            if (x+Ngb.x[k]>=0 && x+Ngb.x[k]<nd2x && y+Ngb.y[k]>=0 && y+Ngb.y[k]<nd2y && z+Ngb.z[k]>=0 && z+Ngb.z[k]<nd2z) {
                                int ngb = Ngb.neighborIndex(k, xyz,nd2x,nd2y,nd2z);
                                if (deformation2Image[ngb + X*nd2xyz]==0 && deformation2Image[ngb + Y*nd2xyz]==0 && deformation2Image[ngb + Z*nd2xyz]==0) {
                                    deformation2Image[ngb + X*nd2xyz] = deformation2Image[xyz + X*nd2xyz];
                                    deformation2Image[ngb + Y*nd2xyz] = deformation2Image[xyz + Y*nd2xyz];
                                    deformation2Image[ngb + Z*nd2xyz] = deformation2Image[xyz + Z*nd2xyz];
                                    growBoundaries = true;
                                    changed[ngb] = true;
                                }
                            }
                        }
                    }
                }
                boundary = changed;
            }
        
            // compose the deformations: X' = def1(def2(X))
            System.out.println("compose deformations");
            float[] composed12 = new float[nd2x*nd2y*nd2z*3];
            for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
                int xyz = x + nd2x*y + nd2x*nd2y*z;
                composed12[xyz+X*nd2xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation2Image[xyz+X*nd2xyz], deformation2Image[xyz+Y*nd2xyz], deformation2Image[xyz+Z*nd2xyz], X, nd1x, nd1y, nd1z, 3);
                composed12[xyz+Y*nd2xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation2Image[xyz+X*nd2xyz], deformation2Image[xyz+Y*nd2xyz], deformation2Image[xyz+Z*nd2xyz], Y, nd1x, nd1y, nd1z, 3);
                composed12[xyz+Z*nd2xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation2Image[xyz+X*nd2xyz], deformation2Image[xyz+Y*nd2xyz], deformation2Image[xyz+Z*nd2xyz], Z, nd1x, nd1y, nd1z, 3);
            }
            deformation = composed12;
            deformation1Image = null;
            deformation2Image = null;

            nrx = nd2x; nry = nd2y; nrz = nd2z;
            rrx = rd2x; rry = rd2y; rrz = rd2z;

            // third deformation if given
            if (!type3Option.equals("none") && deformation3Image!=null) {
                System.out.println("load deformation 3");
                // scale to voxels if needed
                if (type3Option.endsWith("(mm)")) {
                    System.out.println("normalize to resolution ("+rd3x+", "+rd3y+", "+rd3z+")");
                    for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
                        int xyz = x + nd3x*y + nd3x*nd3y*z;
                        deformation3Image[xyz + X*nd3xyz] /= rd3x; 	
                        deformation3Image[xyz + Y*nd3xyz] /= rd3y; 	
                        deformation3Image[xyz + Z*nd3xyz] /= rd3z; 	
                    }
                }
                // turn into a mapping if needed
                if (type3Option.startsWith("deformation")) {
                    for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
                        int xyz = x + nd3x*y + nd3x*nd3y*z;
                        deformation3Image[xyz + X*nd3xyz] += x;
                        deformation3Image[xyz + Y*nd3xyz] += y;
                        deformation3Image[xyz + Z*nd3xyz] += z;
                    }
                }
                // check for bad borders
                boundary = new boolean[nd3x*nd3y*nd3z];
                growBoundaries = false;
                for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
                    int xyz = x + nd3x*y + nd3x*nd3y*z;
                    if (deformation3Image[xyz + X*nd3xyz]==0 && deformation3Image[xyz + Y*nd3xyz]==0 && deformation3Image[xyz + Z*nd3xyz]==0) {
                        for (byte k=0;k<6;k++) {
                            if (x+Ngb.x[k]>=0 && x+Ngb.x[k]<nd3x && y+Ngb.y[k]>=0 && y+Ngb.y[k]<nd3y && z+Ngb.z[k]>=0 && z+Ngb.z[k]<nd3z) {
                                int ngb = Ngb.neighborIndex(k, xyz,nd3x,nd3y,nd3z);
                                if (deformation3Image[ngb + X*nd3xyz]!=0 || deformation3Image[ngb + Y*nd3xyz]!=0 || deformation3Image[ngb + Z*nd3xyz]!=0) {
                                    growBoundaries = true;
                                    boundary[ngb] = true;
                                    k=6;
                                }
                            }
                        }
                    }
                }
                while (growBoundaries) {
                    boolean[] changed = new boolean[nd3x*nd3y*nd3z];
                    growBoundaries = false;
                    for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
                        int xyz = x + nd3x*y + nd3x*nd3y*z;
                        if (boundary[xyz]) {
                            for (byte k=0;k<6;k++) {
                                if (x+Ngb.x[k]>=0 && x+Ngb.x[k]<nd3x && y+Ngb.y[k]>=0 && y+Ngb.y[k]<nd3y && z+Ngb.z[k]>=0 && z+Ngb.z[k]<nd3z) {
                                    int ngb = Ngb.neighborIndex(k, xyz,nd3x,nd3y,nd3z);
                                    if (deformation3Image[ngb + X*nd3xyz]==0 && deformation3Image[ngb + Y*nd3xyz]==0 && deformation3Image[ngb + Z*nd3xyz]==0) {
                                        deformation3Image[ngb + X*nd3xyz] = deformation3Image[xyz + X*nd3xyz];
                                        deformation3Image[ngb + Y*nd3xyz] = deformation3Image[xyz + Y*nd3xyz];
                                        deformation3Image[ngb + Z*nd3xyz] = deformation3Image[xyz + Z*nd3xyz];
                                        growBoundaries = true;
                                        changed[ngb] = true;
                                    }
                                }
                            }
                        }
                    }
                    boundary = changed;
                }
            
                // compose the deformations: X' = def1(def2(def3(X)))
                System.out.println("compose deformations");
                float[] composed123 = new float[nd3x*nd3y*nd3z*3];
                for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
                    int xyz = x + nd3x*y + nd3x*nd3y*z;
                    composed123[xyz+X*nd3xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation3Image[xyz+X*nd3xyz], deformation3Image[xyz+Y*nd3xyz], deformation3Image[xyz+Z*nd3xyz], X, nd2x, nd2y, nd2z, 3);
                    composed123[xyz+Y*nd3xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation3Image[xyz+X*nd3xyz], deformation3Image[xyz+Y*nd3xyz], deformation3Image[xyz+Z*nd3xyz], Y, nd2x, nd2y, nd2z, 3);
                    composed123[xyz+Z*nd3xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation3Image[xyz+X*nd3xyz], deformation3Image[xyz+Y*nd3xyz], deformation3Image[xyz+Z*nd3xyz], Z, nd2x, nd2y, nd2z, 3);
                }
                deformation = composed123;
                deformation3Image = null;
                composed12 = null;
                
                nrx = nd3x; nry = nd3y; nrz = nd3z;
                rrx = rd3x; rry = rd3y; rrz = rd3z;
                
                if (!type4Option.equals("none") && deformation4Image!=null) {
                    System.out.println("load deformation 4");
                    // scale to voxels if needed
                    if (type4Option.endsWith("(mm)")) {
                        System.out.println("normalize to resolution ("+rd4x+", "+rd4y+", "+rd4z+")");
                        for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) for (int z=0;z<nd4z;z++) {
                            int xyz = x + nd4x*y + nd4x*nd4y*z;
                            deformation4Image[xyz + X*nd4xyz] /= rd4x; 	
                            deformation4Image[xyz + Y*nd4xyz] /= rd4y; 	
                            deformation4Image[xyz + Z*nd4xyz] /= rd4z; 	
                        }
                    }
                    // turn into a mapping if needed
                    if (type4Option.startsWith("deformation")) {
                        for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) for (int z=0;z<nd4z;z++) {
                            int xyz = x + nd4x*y + nd4x*nd4y*z;
                            deformation4Image[xyz + X*nd4xyz] += x;
                            deformation4Image[xyz + Y*nd4xyz] += y;
                            deformation4Image[xyz + Z*nd4xyz] += z;
                        }
                    }
                    // check for bad borders
                    boundary = new boolean[nd4x*nd4y*nd4z];
                    growBoundaries = false;
                    for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) for (int z=0;z<nd4z;z++) {
                        int xyz = x + nd4x*y + nd4x*nd4y*z;
                        if (deformation4Image[xyz + X*nd4xyz]==0 && deformation4Image[xyz + Y*nd4xyz]==0 && deformation4Image[xyz + Z*nd4xyz]==0) {
                            for (byte k=0;k<6;k++) {
                                if (x+Ngb.x[k]>=0 && x+Ngb.x[k]<nd4x && y+Ngb.y[k]>=0 && y+Ngb.y[k]<nd4y && z+Ngb.z[k]>=0 && z+Ngb.z[k]<nd4z) {
                                    int ngb = Ngb.neighborIndex(k, xyz,nd4x,nd4y,nd4z);
                                    if (deformation4Image[ngb + X*nd4xyz]!=0 || deformation4Image[ngb + Y*nd4xyz]!=0 || deformation4Image[ngb + Z*nd4xyz]!=0) {
                                        growBoundaries = true;
                                        boundary[ngb] = true;
                                        k=6;
                                    }
                                }
                            }
                        }
                    }
                    while (growBoundaries) {
                        boolean[] changed = new boolean[nd4x*nd4y*nd4z];
                        growBoundaries = false;
                        for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) for (int z=0;z<nd4z;z++) {
                            int xyz = x + nd4x*y + nd4x*nd4y*z;
                            if (boundary[xyz]) {
                                for (byte k=0;k<6;k++) {
                                    if (x+Ngb.x[k]>=0 && x+Ngb.x[k]<nd4x && y+Ngb.y[k]>=0 && y+Ngb.y[k]<nd4y && z+Ngb.z[k]>=0 && z+Ngb.z[k]<nd4z) {
                                        int ngb = Ngb.neighborIndex(k, xyz,nd4x,nd4y,nd4z);
                                        if (deformation4Image[ngb + X*nd4xyz]==0 && deformation4Image[ngb + Y*nd4xyz]==0 && deformation4Image[ngb + Z*nd4xyz]==0) {
                                            deformation4Image[ngb + X*nd4xyz] = deformation4Image[xyz + X*nd4xyz];
                                            deformation4Image[ngb + Y*nd4xyz] = deformation4Image[xyz + Y*nd4xyz];
                                            deformation4Image[ngb + Z*nd4xyz] = deformation4Image[xyz + Z*nd4xyz];
                                            growBoundaries = true;
                                            changed[ngb] = true;
                                        }
                                    }
                                }
                            }
                        }
                        boundary = changed;
                    }
                    // compose the deformations: X' = def1(def2(def3(X)))
                    System.out.println("compose deformations");
                    float[] composed1234 = new float[nd4x*nd4y*nd4z*3];
                    for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) for (int z=0;z<nd4z;z++) {
                        int xyz = x + nd4x*y + nd4x*nd4y*z;
                        composed1234[xyz+X*nd4xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation4Image[xyz+X*nd4xyz], deformation4Image[xyz+Y*nd4xyz], deformation4Image[xyz+Z*nd4xyz], X, nd3x, nd3y, nd3z, 3);
                        composed1234[xyz+Y*nd4xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation4Image[xyz+X*nd4xyz], deformation4Image[xyz+Y*nd4xyz], deformation4Image[xyz+Z*nd4xyz], Y, nd3x, nd3y, nd3z, 3);
                        composed1234[xyz+Z*nd4xyz] = ImageInterpolation.linearClosestInterpolation(deformation, deformation4Image[xyz+X*nd4xyz], deformation4Image[xyz+Y*nd4xyz], deformation4Image[xyz+Z*nd4xyz], Z, nd3x, nd3y, nd3z, 3);
                    }
                    deformation = composed1234;
                    deformation4Image = null;
                    composed123 = null;
                    
                    nrx = nd4x; nry = nd4y; nrz = nd4z;
                    rrx = rd4x; rry = rd4y; rrz = rd4z;
                }
            }
        }
        nrxyz = nrx*nry*nrz;
        
        System.out.println("output dimensions: "+nrx+" x "+nry+" x "+nrz+"("+nst+")");
        
        // new image
        System.out.println("deform image");
        float min = 1e10f, max = -1e10f;
        if (padOption.equals("min") || padOption.equals("max")) {
            for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) for (int t=0;t<nst;t++) {
                int xyzt = x + nsx*y + nsx*nsy*z + nsx*nsy*nsz*t;
                if (sourceImage[xyzt]<min) min = sourceImage[xyzt];
                if (sourceImage[xyzt]>max) max = sourceImage[xyzt];
            }
        }
        deformedImage = new float[nrx*nry*nrz*nst];
        if (interpOption.equals("nearest")) {
            for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
                int xyz = x + nrx*y + nrx*nry*z;
                for (int t=0;t<nst;t++) {
                    if (padOption.equals("closest"))
                        deformedImage[xyz + nrxyz*t] = ImageInterpolation.nearestNeighborClosestInterpolation(sourceImage, deformation[xyz+X*nrxyz], deformation[xyz+Y*nrxyz], deformation[xyz+Z*nrxyz], t, nsx, nsy, nsz, nst);
                    else if (padOption.equals("zero"))
                        deformedImage[xyz + nrxyz*t] = ImageInterpolation.nearestNeighborInterpolation(sourceImage, 0.0f, deformation[xyz+X*nrxyz], deformation[xyz+Y*nrxyz], deformation[xyz+Z*nrxyz], t, nsx, nsy, nsz, nst);
                    else if (padOption.equals("min"))
                        deformedImage[xyz + nrxyz*t] = ImageInterpolation.nearestNeighborInterpolation(sourceImage, min, deformation[xyz+X*nrxyz], deformation[xyz+Y*nrxyz], deformation[xyz+Z*nrxyz], t, nsx, nsy, nsz, nst);
                    else if (padOption.equals("max"))
                        deformedImage[xyz + nrxyz*t] = ImageInterpolation.nearestNeighborInterpolation(sourceImage, max, deformation[xyz+X*nrxyz], deformation[xyz+Y*nrxyz], deformation[xyz+Z*nrxyz], t, nsx, nsy, nsz, nst); 
                }
            }
        } else if (interpOption.equals("linear")) {
            for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
                int xyz = x + nrx*y + nrx*nry*z;
                for (int t=0;t<nst;t++) {
                    if (padOption.equals("closest"))
                        deformedImage[xyz + nrxyz*t] = ImageInterpolation.linearClosestInterpolation(sourceImage, deformation[xyz+X*nrxyz], deformation[xyz+Y*nrxyz], deformation[xyz+Z*nrxyz], t, nsx, nsy, nsz, nst);
                    else if (padOption.equals("zero"))
                        deformedImage[xyz + nrxyz*t] = ImageInterpolation.linearInterpolation(sourceImage, 0.0f, deformation[xyz+X*nrxyz], deformation[xyz+Y*nrxyz], deformation[xyz+Z*nrxyz], t, nsx, nsy, nsz, nst);
                    else if (padOption.equals("min"))
                        deformedImage[xyz + nrxyz*t] = ImageInterpolation.linearInterpolation(sourceImage, min, deformation[xyz+X*nrxyz], deformation[xyz+Y*nrxyz], deformation[xyz+Z*nrxyz], t, nsx, nsy, nsz, nst);
                    else if (padOption.equals("max"))
                        deformedImage[xyz + nrxyz*t] = ImageInterpolation.linearInterpolation(sourceImage, max, deformation[xyz+X*nrxyz], deformation[xyz+Y*nrxyz], deformation[xyz+Z*nrxyz], t, nsx, nsy, nsz, nst);
                }
            }
        }
        sourceImage = null;
        deformation = null;
    }
}
