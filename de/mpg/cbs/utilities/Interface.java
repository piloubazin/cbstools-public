package de.mpg.cbs.utilities;


import java.awt.*;
import java.awt.event.*;
import java.net.URL;

import javax.swing.*;

import javax.vecmath.*;

import java.io.*;
import java.net.*;
import java.util.*;
import java.util.jar.*;
import java.util.zip.*;

import gov.nih.mipav.view.Preferences;
import gov.nih.mipav.view.MipavUtil;
import gov.nih.mipav.view.ViewUserInterface;
import gov.nih.mipav.model.structures.*;
import gov.nih.mipav.model.file.*;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataDouble;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurfaceCollection;


// TODO: Auto-generated Javadoc
/**
 * A collection of utilities for I/O, GUI, interface with MIPAV, etc.
 * 
 * @author Pierre-Louis Bazin
 * @author Blake Lucas
 */
public class Interface {
//	protected static ViewUserInterface userInterface=null;
	/** The quiet. */
	//protected static final boolean quiet=true;

	/**
	 * Sets the quiet.
	 * 
	 * @param q the new quiet
	 */
	//public static final void setQuiet(boolean q){
	//	quiet=q;
	//}
	
	/**
	 * Checks if is quiet.
	 * 
	 * @return true, if is quiet
	 */
	//public static final boolean isQuiet(){
	//	return quiet;
	//}
	
	/**
	 * Singleton method to get User BasicInfo.
	 * WARNING: DO NOT USE THIS WITH JIST 
	 * 
	 * @return userInterface
	 */
	public static ViewUserInterface getUI(){
		return ViewUserInterface.getReference();		
	}
	
	/**
	 * Display message.
	 * 
	 * @param message the message
	 */
	public static final void displayMessage(String message){
		// GUI output
		//if (!isQuiet()) getUI().setGlobalDataText(message);
		// console output
		System.out.print(message);
		System.out.flush();
	}
	
	/**
	 * Display error.
	 * 
	 * @param message the message
	 */
	public static final void displayError(String message){
		// GUI output
		//if (!isQuiet()) MipavUtil.displayError(message);
		// console output
		System.err.print(message);
		System.err.flush();
	}
    
    /** Displays the Java Help dialog indexed directly to the section identified by the ID passed in. */
    //static HelpSet hs;
    
    /** The help broker. */
    //static HelpBroker helpBroker;
    
    /**
	 * Copy important file information between ModelImage structures,
	 * assuming all slices have same properties (uses only the first slice from
	 * the source).
	 * 
	 * @param image the image
	 * @param resultImage the result image
	 */
    public static final void updateFileInfo(ModelImage image, ModelImage resultImage) {
        FileInfoBase[] fileInfo;

        if (resultImage.getNDims() == 2) {
            fileInfo = resultImage.getFileInfo();
            
			fileInfo[0].setModality(image.getFileInfo()[0].getModality());
            fileInfo[0].setFileDirectory(image.getFileInfo()[0].getFileDirectory());
			fileInfo[0].setEndianess(image.getFileInfo()[0].getEndianess());
            fileInfo[0].setUnitsOfMeasure(image.getFileInfo()[0].getUnitsOfMeasure());
            fileInfo[0].setResolutions(image.getFileInfo()[0].getResolutions());
            fileInfo[0].setAxisOrientation(image.getFileInfo()[0].getAxisOrientation());
            fileInfo[0].setOrigin(image.getFileInfo()[0].getOrigin());
            fileInfo[0].setPixelPadValue(image.getFileInfo()[0].getPixelPadValue());
            fileInfo[0].setPhotometric(image.getFileInfo()[0].getPhotometric());
			
			fileInfo[0].setImageOrientation(image.getImageOrientation());
            
			fileInfo[0].setExtents(resultImage.getExtents());
            fileInfo[0].setMax(resultImage.getMax());
            fileInfo[0].setMin(resultImage.getMin());
            
        } else if (resultImage.getNDims() == 3) {
			//System.out.print("3:");
            fileInfo = resultImage.getFileInfo();

            for (int i = 0; i < resultImage.getExtents()[2]; i++) {
                fileInfo[i].setModality(image.getFileInfo()[0].getModality());
                fileInfo[i].setFileDirectory(image.getFileInfo()[0].getFileDirectory());
				fileInfo[i].setEndianess(image.getFileInfo()[0].getEndianess());
                fileInfo[i].setUnitsOfMeasure(image.getFileInfo()[0].getUnitsOfMeasure());
                fileInfo[i].setResolutions(image.getFileInfo()[0].getResolutions());
                fileInfo[i].setAxisOrientation(image.getFileInfo()[0].getAxisOrientation());
                fileInfo[i].setOrigin(image.getFileInfo()[0].getOrigin());
                fileInfo[i].setPixelPadValue(image.getFileInfo()[0].getPixelPadValue());
                fileInfo[i].setPhotometric(image.getFileInfo()[0].getPhotometric());
               
		        fileInfo[i].setImageOrientation(image.getImageOrientation());
				
				fileInfo[i].setExtents(resultImage.getExtents());
                fileInfo[i].setMax(resultImage.getMax());
                fileInfo[i].setMin(resultImage.getMin());
            }
        } else if (resultImage.getNDims() == 4) {
            //System.out.print("4:");
            fileInfo = resultImage.getFileInfo();

			int[] units = new int[4];
			float[] res = new float[4];
			for (int n=0;n<4;n++) {
				if (n<image.getNDims()) {
					units[n] = image.getFileInfo()[0].getUnitsOfMeasure()[n];
					res[n] = image.getFileInfo()[0].getResolutions()[n];
				} else {
					units[n] = image.getFileInfo()[0].getUnitsOfMeasure()[image.getNDims()-1];
					res[n] = image.getFileInfo()[0].getResolutions()[image.getNDims()-1];
				}					
			}
				
            for (int i = 0; i < (resultImage.getExtents()[2] * resultImage.getExtents()[3]); i++) {
                fileInfo[i].setModality(image.getFileInfo()[0].getModality());
                fileInfo[i].setFileDirectory(image.getFileInfo()[0].getFileDirectory());
                fileInfo[i].setEndianess(image.getFileInfo()[0].getEndianess());
                fileInfo[i].setAxisOrientation(image.getFileInfo()[0].getAxisOrientation());
                fileInfo[i].setOrigin(image.getFileInfo()[0].getOrigin());
                fileInfo[i].setPixelPadValue(image.getFileInfo()[0].getPixelPadValue());
                fileInfo[i].setPhotometric(image.getFileInfo()[0].getPhotometric());
				
				fileInfo[i].setUnitsOfMeasure(units);
                fileInfo[i].setResolutions(res);
                
				fileInfo[i].setImageOrientation(image.getImageOrientation());
                
				fileInfo[i].setExtents(resultImage.getExtents());
                fileInfo[i].setMax(resultImage.getMax());
                fileInfo[i].setMin(resultImage.getMin());
            }
        }
    }
    
	/**
	 * Copy important file information between ModelImage structures,
	 * assuming all slices have same properties (uses only the first slice from
	 * the source).
	 * 
	 * @param info the info
	 * @param resultImage the result image
	 */
    public static final void updateFileInfo(FileInfoBase info, ModelImage resultImage) {
        FileInfoBase[] fileInfo;

        if (resultImage.getNDims() == 2) {
            fileInfo = resultImage.getFileInfo();
            
			fileInfo[0].setModality(info.getModality());
            fileInfo[0].setFileDirectory(info.getFileDirectory());
			fileInfo[0].setEndianess(info.getEndianess());
            fileInfo[0].setUnitsOfMeasure(info.getUnitsOfMeasure());
            fileInfo[0].setResolutions(info.getResolutions());
            fileInfo[0].setAxisOrientation(info.getAxisOrientation());
            fileInfo[0].setOrigin(info.getOrigin());
            fileInfo[0].setPixelPadValue(info.getPixelPadValue());
            fileInfo[0].setPhotometric(info.getPhotometric());
			
			fileInfo[0].setImageOrientation(info.getImageOrientation());
            
			fileInfo[0].setExtents(resultImage.getExtents());
            fileInfo[0].setMax(resultImage.getMax());
            fileInfo[0].setMin(resultImage.getMin());
            
        } else if (resultImage.getNDims() == 3) {
			//System.out.print("3:");
            fileInfo = resultImage.getFileInfo();

            for (int i = 0; i < resultImage.getExtents()[2]; i++) {
                fileInfo[i].setModality(info.getModality());
                fileInfo[i].setFileDirectory(info.getFileDirectory());
				fileInfo[i].setEndianess(info.getEndianess());
                fileInfo[i].setUnitsOfMeasure(info.getUnitsOfMeasure());
                fileInfo[i].setResolutions(info.getResolutions());
                fileInfo[i].setAxisOrientation(info.getAxisOrientation());
                fileInfo[i].setOrigin(info.getOrigin());
                fileInfo[i].setPixelPadValue(info.getPixelPadValue());
                fileInfo[i].setPhotometric(info.getPhotometric());
               
		        fileInfo[i].setImageOrientation(info.getImageOrientation());
				
				fileInfo[i].setExtents(resultImage.getExtents());
                fileInfo[i].setMax(resultImage.getMax());
                fileInfo[i].setMin(resultImage.getMin());
            }
        } else if (resultImage.getNDims() == 4) {
            //System.out.print("4:");
            fileInfo = resultImage.getFileInfo();

			int[] units = new int[4];
			float[] res = new float[4];
			for (int n=0;n<4;n++) {
				units[n] = info.getUnitsOfMeasure()[n];
				res[n] = info.getResolutions()[n];
			}
				
            for (int i = 0; i < (resultImage.getExtents()[2] * resultImage.getExtents()[3]); i++) {
                fileInfo[i].setModality(info.getModality());
                fileInfo[i].setFileDirectory(info.getFileDirectory());
                fileInfo[i].setEndianess(info.getEndianess());
                fileInfo[i].setAxisOrientation(info.getAxisOrientation());
                fileInfo[i].setOrigin(info.getOrigin());
                fileInfo[i].setPixelPadValue(info.getPixelPadValue());
                fileInfo[i].setPhotometric(info.getPhotometric());
				
				fileInfo[i].setUnitsOfMeasure(units);
                fileInfo[i].setResolutions(res);
                
				fileInfo[i].setImageOrientation(info.getImageOrientation());
                
				fileInfo[i].setExtents(resultImage.getExtents());
                fileInfo[i].setMax(resultImage.getMax());
                fileInfo[i].setMin(resultImage.getMin());
            }
        }
    }
    
    /** return the maximum common base name for two strings (assumes the same beginning) */
    public static final String commonNameBase(String name1, String name2) {
    	int id=1;
    	while (name1.regionMatches(0, name2, 0, id)) id++;
    	
    	return name1.substring(0,id-1);
    }
    /** return the part of the base name not present in the second string (assumes the same beginning) */
    public static final String differentNameEnding(String name1, String name2) {
    	int id=1;
    	while (name1.regionMatches(0, name2, 0, id)) id++;
    	
    	return name1.substring(id-1);
    }
    
    // internal copy of MIPAV's internal orientation, etc. labels
    // hard-coded...
    public static final int AXIAL = 	0;	// FileInfoBase.AXIAL;
    public static final int CORONAL = 	1;	// FileInfoBase.CORONAL;
    public static final int SAGITTAL = 	2;	// FileInfoBase.SAGITTAL;
    
    public static final int R2L = 1; // FileInfoBase.ORI_R2L_TYPE;
    public static final int L2R = 2; // FileInfoBase.ORI_L2R_TYPE;
    public static final int P2A = 3; // FileInfoBase.ORI_P2A_TYPE;
    public static final int A2P = 4; // FileInfoBase.ORI_A2P_TYPE;
    public static final int I2S = 5; // FileInfoBase.ORI_I2S_TYPE;
    public static final int S2I = 6; // FileInfoBase.ORI_S2I_TYPE;
    
    // imported / adapted convenience classes from MIPAV, copied for self-containedness
    public static final boolean getBoolean(final StringTokenizer st) {
        final String str = st.nextToken();

        // returns true if str.equalsIgnoreCase( "true" )
        return new Boolean(str).booleanValue();
    }

    public static final float getFloat(final StringTokenizer st) {
        final String str = st.nextToken();

        return Float.parseFloat(str);
    }

    public static final double getDouble(final StringTokenizer st) {
        final String str = st.nextToken();

        return Double.parseDouble(str);
    }

    public static final int getInt(final StringTokenizer st) {
        final String str = st.nextToken();

        return Integer.parseInt(str);
    }

    // inputs
    public static final String getName(ParamVolume input) {
    	return input.getImageData().getName();
    }
    
    public static final ImageHeader getHeader(ParamVolume input) {
    	return input.getImageData().getHeader();
    }
    
    public static final int[] getDimensions(ParamVolume input) {
		int[] dims = new int[3];
		dims[0] = input.getImageData().getRows();
		dims[1] = input.getImageData().getCols();
		dims[2] = input.getImageData().getSlices();
		return dims;
	}
	
    public static final int[] getDimensions4D(ParamVolume input) {
		int[] dims = new int[4];
		dims[0] = input.getImageData().getRows();
		dims[1] = input.getImageData().getCols();
		dims[2] = input.getImageData().getSlices();
		dims[3] = input.getImageData().getComponents();
		return dims;
	}
	
    public static final int getComponents(ParamVolume input) {
    	return input.getImageData().getComponents();
	}
	
	public static final boolean isImage4D(ParamVolume input) {
		return (input.getImageData().getComponents()!=1);
	}
	
	public static final boolean isValid(ParamVolume input) {
		return (input.getImageData()!=null);
	}
	
    public static final float[] getResolutions(ParamVolume input) {
    	float[] res = new float[3];
    	res[0] = input.getImageData().getHeader().getDimResolutions()[0];
		res[1] = input.getImageData().getHeader().getDimResolutions()[1];
		res[2] = input.getImageData().getHeader().getDimResolutions()[2];
		return res;
	}
	
	public static final int[] getOrientations(ParamVolume input) {
		int[] ori = new int[4];
		ori[0] = input.getImageData().getHeader().getImageOrientation().ordinal();
		ori[1] = input.getImageData().getHeader().getAxisOrientation()[0].ordinal();
		ori[2] = input.getImageData().getHeader().getAxisOrientation()[1].ordinal();
		ori[3] = input.getImageData().getHeader().getAxisOrientation()[2].ordinal();
		return ori;
	}
	
	public static final String getName(ParamSurface input) {
    	return input.getSurface().getName();
    }
    
    public static final float[] getSurfacePoints(ParamSurface surface) {
		int npt = surface.getSurface().getVertexCount();
		float[] points = new float[3*npt];
		for(int i=0; i<npt; i++){
			Point3f p = surface.getSurface().getVertex(i);
			points[3*i+0] = p.x;
			points[3*i+1] = p.y;
			points[3*i+2] = p.z;
		}
		return points;
	}
	
	public static final int[] getSurfaceTriangles(ParamSurface surface) {
		int ntr = surface.getSurface().getFaceCount();
		int[] triangles = new int[3*ntr];
		for(int i=0; i<ntr; i++){
			int[] tri = surface.getSurface().getFaceVertexIds(i);
			triangles[3*i+0] = tri[0];
			triangles[3*i+1] = tri[1];
			triangles[3*i+2] = tri[2];
		}
		return triangles;
	}
	
	public static final void setSurface(float[] points, int[] triangles, ParamSurface container, String name) {
		int npt = points.length/3;
		Point3f[] pts = new Point3f[npt];
		for(int i=0; i<npt; i++){
			pts[i] = new Point3f(points[3*i],points[3*i+1],points[3*i+2]);
		}
		EmbeddedSurface surf = new EmbeddedSurface(pts, triangles);
		surf.setName(name);
		container.setValue(surf);
	}
	
	public static final ParamVolume getVolumeFromCollection(int n, ParamVolumeCollection input) {
	    return input.getParamVolume(n);
	}
	
	public static final int getVolumeCollectionSize(ParamVolumeCollection input) {
	    return input.size();
	}
	
    public static final float[] getFloatImage3D(ParamVolume input) {
		// import the image data into 1D arrays
		ImageDataFloat inImg = new ImageDataFloat(input.getImageData());
		float[][][] buffer = inImg.toArray3d();
		int nx = inImg.getRows();
		int ny = inImg.getCols();
		int nz = inImg.getSlices();
		//int nxyz = nx*ny*nz;
		
		float[] image = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[xyz] = buffer[x][y][z];
		}
		return image;
	}

    public static final float[] getFloatImageComponent4D(ParamVolume input, int c) {
		// import the image data into 1D arrays
		ImageDataFloat inImg = new ImageDataFloat(input.getImageData());
		float[][][][] buffer = inImg.toArray4d();
		int nx = inImg.getRows();
		int ny = inImg.getCols();
		int nz = inImg.getSlices();
		//int nxyz = nx*ny*nz;
		
		float[] image = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[xyz] = buffer[x][y][z][c];
		}
		return image;
	}

    public static final float[] getFloatImage4D(ParamVolume input) {
		// import the image data into 1D arrays
		ImageDataFloat inImg = new ImageDataFloat(input.getImageData());
		float[][][][] buffer = inImg.toArray4d();
		int nx = inImg.getRows();
		int ny = inImg.getCols();
		int nz = inImg.getSlices();
		int nc = inImg.getComponents();
		//int nxyz = nx*ny*nz*nc;
		
		float[] image = new float[nx*ny*nz*nc];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int c=0;c<nc;c++) {
			int xyz = x+nx*y+nx*ny*z+nx*ny*nz*c;
			image[xyz] = buffer[x][y][z][c];
		}
		return image;
	}

    public static final byte[] getUByteImage3D(ParamVolume input) {
		if (input.getImageData()!=null) {  	
			// import the image data into 1D arrays
			ImageDataUByte inImg = new ImageDataUByte(input.getImageData());
			byte[][][] buffer = inImg.toArray3d();
			int nx = inImg.getRows();
			int ny = inImg.getCols();
			int nz = inImg.getSlices();
			//int nxyz = nx*ny*nz;
			
			byte[] image = new byte[nx*ny*nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[xyz] = buffer[x][y][z];
			}
			return image;
		} else {
			return null;
		}
	}

    public static final int[] getIntegerImage3D(ParamVolume input) {
		// import the image data into 1D arrays
		ImageDataInt inImg = new ImageDataInt(input.getImageData());
		int[][][] buffer = inImg.toArray3d();
		int nx = inImg.getRows();
		int ny = inImg.getCols();
		int nz = inImg.getSlices();
		//int nxyz = nx*ny*nz;
		
		int[] image = new int[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[xyz] = buffer[x][y][z];
		}
		return image;
	}

   public static final int[] getIntegerImage4D(ParamVolume input) {
		// import the image data into 1D arrays
		ImageDataInt inImg = new ImageDataInt(input.getImageData());
		int[][][][] buffer = inImg.toArray4d();
		int nx = inImg.getRows();
		int ny = inImg.getCols();
		int nz = inImg.getSlices();
		int nc = inImg.getComponents();
		//int nxyz = nx*ny*nz*nc;
		
		int[] image = new int[nx*ny*nz*nc];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int c=0;c<nc;c++) {
			int xyz = x+nx*y+nx*ny*z+nx*ny*nz*c;
			image[xyz] = buffer[x][y][z][c];
		}
		return image;
	}

	// outputs
	public static final void setFloatImage3D(float[] image, int[] dimensions, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		float[][][] buffer = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			buffer[x][y][z] = image[xyz];
		}
		
		ImageDataFloat bufferData = new ImageDataFloat(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}
	
	// outputs
	public static final void setFloatImage4D(float[] image, int[] dimensions, int length, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		float[][][][] buffer = new float[nx][ny][nz][length];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<length;l++) {
			int xyz = x+nx*y+nx*ny*z+nx*ny*nz*l;
			buffer[x][y][z][l] = image[xyz];
		}
		
		ImageDataFloat bufferData = new ImageDataFloat(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}
	
	public static final void setUByteImage3D(byte[] image, int[] dimensions, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		byte[][][] buffer = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			buffer[x][y][z] = image[xyz];
		}
		
		ImageDataUByte bufferData = new ImageDataUByte(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}
	
	public static final void setUByteImage4D(byte[] image, int[] dimensions, int length, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		byte[][][][] buffer = new byte[nx][ny][nz][length];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<length;l++) {
			int xyz = x+nx*y+nx*ny*z+nx*ny*nz*l;
			buffer[x][y][z][l] = image[xyz];
		}
		
		ImageDataUByte bufferData = new ImageDataUByte(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}

	public static final void setByteImage3D(byte[] image, int[] dimensions, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		byte[][][] buffer = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			buffer[x][y][z] = image[xyz];
		}
		
		ImageDataByte bufferData = new ImageDataByte(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}
	
	public static final void setByteImage4D(byte[] image, int[] dimensions, int length, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		byte[][][][] buffer = new byte[nx][ny][nz][length];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<length;l++) {
			int xyz = x+nx*y+nx*ny*z+nx*ny*nz*l;
			buffer[x][y][z][l] = image[xyz];
		}
		
		ImageDataByte bufferData = new ImageDataByte(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}
	public static final void setIntegerImage3D(int[] image, int[] dimensions, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		int[][][] buffer = new int[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			buffer[x][y][z] = image[xyz];
		}
		
		ImageDataInt bufferData = new ImageDataInt(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}
	
	public static final void setIntegerImage4D(int[] image, int[] dimensions, int length, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		int[][][][] buffer = new int[nx][ny][nz][length];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<length;l++) {
			int xyz = x+nx*y+nx*ny*z+nx*ny*nz*l;
			buffer[x][y][z][l] = image[xyz];
		}
		
		ImageDataInt bufferData = new ImageDataInt(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}


	// outputs
	public static final void setFloatImage2D(float[] image, int[] dimensions, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		//int nxyz = nx*ny*nz;
		
		float[][] buffer = new float[nx][ny];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
			int xyz = x+nx*y;
			buffer[x][y] = image[xyz];
		}
		
		ImageDataFloat bufferData = new ImageDataFloat(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}
	// outputs
	public static final void setDoubleImage2D(double[] image, int[] dimensions, ParamVolume container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		//int nxyz = nx*ny*nz;
		
		double[][] buffer = new double[nx][ny];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
			int xyz = x+nx*y;
			buffer[x][y] = image[xyz];
		}
		
		ImageDataDouble bufferData = new ImageDataDouble(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.setValue(bufferData);
		
		return;
	}
	
	// outputs
	public static final void addFloatImage4D(float[] image, int[] dimensions, int length, ParamVolumeCollection container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		float[][][][] buffer = new float[nx][ny][nz][length];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<length;l++) {
			int xyz = x+nx*y+nx*ny*z+nx*ny*nz*l;
			buffer[x][y][z][l] = image[xyz];
		}
		
		ImageDataFloat bufferData = new ImageDataFloat(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.add(bufferData);
		
		return;
	}
	
	// outputs
	public static final void addFloatImage3D(float[] image, int[] dimensions, ParamVolumeCollection container, String name, ImageHeader header) {
		int nx = dimensions[0];
		int ny = dimensions[1];
		int nz = dimensions[2];
		//int nxyz = nx*ny*nz;
		
		float[][][] buffer = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			buffer[x][y][z] = image[xyz];
		}
		
		ImageDataFloat bufferData = new ImageDataFloat(buffer);
		bufferData.setHeader(header);
		bufferData.setName(name);
		container.add(bufferData);
		
		return;
	}
	

}
