package de.mpg.cbs.utilities;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.*;
// TODO: Auto-generated Javadoc

/**
 * Utility for Cropping images.
 * 
 * @author Navid Shiee
 */

public class CubicVolumeCropper extends AbstractCalculation{

	/*
	private	int x0;
	private int xN;
	private int y0;
	private int yN;
	private int z0;
	private int zN;
*/
	/**
	 * Instantiates a new cubic volume cropper.
	 */
	public CubicVolumeCropper (){
		super();
		/*
		x0=0;
		xN=0;
		y0=0;
		yN=0;
		z0=0;
		zN=0;
	*/
	}

	/**
	 * Providing the boundries.
	 * 
	 * @param _f the _f
	 * @param _val the _val
	 * @param _border the _border
	 * 
	 * @return the crop parameters
	 */
/*
	public int[] boundingBox() {
		int[] box = new int[6];
		box[0] = x0; box[1] = xN;
		box[2] = y0; box[3] = yN;
		box[4] = z0; box[5] = zN;
		return box;
	}    
*/
	/**
	 * This method finds the bounderies of image and adding given border to
	 * compute desired extents of image 
	 */

	public CropParameters findBoundaries(ImageData _f, double _val, int _border){
		int nx = _f.getRows();
		int ny = _f.getCols();
		int nz = _f.getSlices();
		int nc = _f.getComponents();
		int x0 = nx;
		int xN = 0;
		int y0 = ny;
		int yN = 0;
		int z0 = nz;
		int zN = 0;
		if (nc==1){
			for (int x=0; x<nx; x++){
				for (int y=0; y<ny; y++){
					for (int z=0; z<nz; z++){
						
						if (_f.getDouble(x, y, z)>_val) {
							
							if (x < x0) x0 = x;
							if (x > xN) xN = x;
							if (y < y0) y0 = y;
							if (y > yN) yN = y;
							if (z < z0) z0 = z;
							if (z > zN) zN = z;
						} 
					}
				}
			}
			incrementCompletedUnits();
		}
		else{

			for (int c=0; c<nc; c++ ){
				for (int x=0; x<nx; x++){
					for (int y=0; y<ny; y++){
						for (int z=0; z<nz; z++){
							if (_f.get(x,y,z,c).doubleValue()>_val) {
								if (x < x0) x0 = x;
								if (x > xN) xN = x;
								if (y < y0) y0 = y;
								if (y > yN) yN = y;
								if (z < z0) z0 = z;
								if (z > zN) zN = z;
							}
						}
					}
				}
				incrementCompletedUnits();
			}
		}
		x0 -= _border;
		xN += _border;
		y0 -= _border;
		yN += _border;
		z0 -= _border;
		zN += _border;   	
		return new CropParameters(_f.getRows(),_f.getCols(),_f.getSlices(),_f.getComponents(),x0,xN,y0,yN,z0,zN,_val);
	}


	/**
	 * Method assumes all input images have the same initail dimensions.
	 * Boundaries are computed with respect to extents of all images.
	 * 
	 * @param _f the _f
	 * @param _val the _val
	 * @param _border the _border
	 * 
	 * @return the crop parameters
	 */

	private CropParameters findBoundaries(List<ImageData> _f,double _val, int _border){

		int min_x0 = Integer.MAX_VALUE;
		int max_xN = 0;
		int min_y0 = Integer.MAX_VALUE;
		int max_yN = 0;
		int min_z0 = Integer.MAX_VALUE;
		int max_zN = 0;
		CropParameters p=null;
		for (ImageData vol:_f){
			p=findBoundaries(vol,_val,_border);
			min_x0= Math.min(min_x0, p.xmin);
			max_xN= Math.max(max_xN, p.xmax);
			min_y0= Math.min(min_y0, p.ymin);
			max_yN= Math.max(max_yN, p.ymax);
			min_z0= Math.min(min_z0, p.zmin);
			max_zN= Math.max(max_zN, p.zmax);
		}
		if(p!=null){
			return new CropParameters(p.rows,p.cols,p.slices,p.components,min_x0,max_xN,min_y0,max_yN,min_z0,max_zN,_val);
		} else return null;
	}


	/**
	 * This method make a cropped version of input image
	 * BECAUSE set(Number) is no longer supported, will not work for color volumes.
	 * 
	 * @param f the f
	 * @param params the params
	 * 
	 * @return the image data
	 */
	public  ImageData crop(ImageData f,CropParameters params){
		String name=f.getName();
		System.out.println(getClass().getCanonicalName()+"\t"+"Cropping Image...");
		int x0=params.xmin;
		int xN=params.xmax;
		int y0=params.ymin;
		int yN=params.ymax;
		int z0=params.zmin;
		int zN=params.zmax;
		double val=params.background;
		//System.out.println(getClass().getCanonicalName()+"\t"+"Original dimension :(" + f.getRows() + ", " + f.getCols() + ", " + f.getSlices() + ")" );
		//System.out.println(getClass().getCanonicalName()+"\t"+"Extents after adding border: (" + x0+ ", " + xN + ") (" +
			//	y0 + ", " + yN + ") (" + z0 + ", " + zN +")");
		//Cropped image dimensions:
		int mx = xN-x0+1;
		int my = yN-y0+1;
		int mz = zN-z0+1;
		int mc = f.getComponents();
		ImageDataMipav V;
		if (mc == 1){					
			V = new ImageDataMipav(name+"_cp",f.getType(),mx, my, mz);			
			for (int x=0; x<mx; x++){
				for (int y=0; y<my; y++){
					for (int z=0; z<mz; z++){
						if ( (x+x0 >= 0) && (x+x0 < f.getRows()) && (y+y0 >= 0) && (y+y0 < f.getCols()) && (z+z0 >= 0) && (z+z0 < f.getSlices()))
							V.set(x, y, z, f.getFloat(x+x0, y+y0, z+z0)) ;
						else
							V.set(x,y,z,params.background);
					}
				}
			}
		}
		else{			
			V = new ImageDataMipav(name+"_cp",f.getType(),mx, my, mz, mc);
			System.out.println("output volume: " + mx + " " + my + " " + mz + " " + mc);
			System.out.println("offset: " + x0 + " " + y0 + " " + z0 );
			for (int c=0;c<mc;c++){
				for (int x=0; x<mx; x++){
					for (int y=0; y<my; y++){
						for (int z=0; z<mz; z++){
							if ( (x+x0 >= 0) && (x+x0 < f.getRows()) && (y+y0 >= 0) && (y+y0 < f.getCols()) && (z+z0 >= 0) && (z+z0 < f.getSlices()))
								V.set(x, y, z, c, f.getDouble(x+x0, y+y0, z+z0, c)) ;
							else
								V.set(x,y,z,params.background);
							
						}
					}
				}
			}
		}		

//		V.setCropParameters(params);
		cropParams = params;
		V.setHeader(f.getHeader().clone());
		return V;

	}
	
	/** The crop params. */
	CropParameters cropParams;
	
	/**
	 * Gets the last crop params.
	 * 
	 * @return the last crop params
	 */
	public CropParameters getLastCropParams(){
		return cropParams;
	}
	
	/**
	 * Crop.
	 * 
	 * @param f the f
	 * @param val the val
	 * @param border the border
	 * 
	 * @return the image data
	 */
	public  ImageData crop(ImageData f,double val, int border){
		return crop(f,findBoundaries(f, val, border));
	}
	/*
	public  CubicVolume crop(CubicVolume f,double val, int border){

		String name=f.getName();
		System.out.println(getClass().getCanonicalName()+"\t"+"Cropping Image...");
		CropParameters params=findBounderies(f, val,border);
		System.out.println(getClass().getCanonicalName()+"\t"+"Original dimension :(" + f.getRows() + ", " + f.getCols() + ", " + f.getSlices() + ")" );
		System.out.println(getClass().getCanonicalName()+"\t"+"Extents after adding border: (" + x0+ ", " + xN + ") (" +y0 + ", " + yN + ") (" + z0 + ", " + zN +")");
		//Cropped image dimensions:
		int mx = xN-x0+1;
		int my = yN-y0+1;
		int mz = zN-z0+1;
		int mc = f.getComponents();
		CubicVolumeMipav V;
		if (mc == 1){
			V = new CubicVolumeMipav(f.getType(),mx, my, mz);
			for (int x=0; x<mx; x++){
				for (int y=0; y<my; y++){
					for (int z=0; z<mz; z++){
						V.set(x, y, z, f.get(x+x0, y+y0, z+z0)) ;
					}
				}
			}
			incrementCompletedUnits();
		}
		else{
			V = new CubicVolumeMipav(f.getType(),mx, my, mz, mc);
			for (int c=0;c<mc;c++){
				for (int x=0; x<mx; x++){
					for (int y=0; y<my; y++){
						for (int z=0; z<mz; z++){
							V.set(x, y, z, c, f.get(x+x0, y+y0, z+z0, c)) ;
						}
					}
				}
				incrementCompletedUnits();
			}
		}

		V.setName(name);
		V.setCropParameters(params);
		markCompleted();
		return V;

	}
	*/
	/**
	 * Crop.
	 * 
	 * @param imgs the imgs
	 * @param val the val
	 * @param border the border
	 * 
	 * @return the list< image data>
	 */
	public List<ImageData> crop(List<ImageData> imgs,double val, int border){
		CropParameters params=findBoundaries(imgs, val,border);
		return crop(imgs,params);
	}
	
	public List<ImageData> crop(List<ImageData> imgs,CropParameters params){
		ArrayList<ImageData> result=new ArrayList<ImageData>();
		setTotalUnits(imgs.size());
		for(ImageData img:imgs){
			incrementCompletedUnits();
			result.add(crop(img,params));
		}
		markCompleted();
		return result;
	}
	
	/**
	 * Uncrop.
	 * 
	 * @param imgs the imgs
	 * @param params the params
	 * 
	 * @return the list< image data>
	 */
	public List<ImageData> uncrop(List<ImageData> imgs,CropParameters params){
		ArrayList<ImageData> result=new ArrayList<ImageData>();
		setTotalUnits(imgs.size());
		for(ImageData img:imgs){
			incrementCompletedUnits();
			result.add(uncrop(img,params));
		}
		markCompleted();
		return result;
	}
	
	/**
	 * Uncrop.
	 * 
	 * @param vol the vol
	 * @param params the params
	 * 
	 * @return the image data
	 */
	public ImageData uncrop(ImageData vol,CropParameters params){
		
		//Cropped image dimensions:
		System.out.println(getClass().getCanonicalName()+"\t"+"UNCROP "+vol.getName());
		ImageDataMipav V;
		System.out.println(getClass().getCanonicalName()+"\t"+params);
		double val=0;
		if (params.components == 1){
			V = new ImageDataMipav(vol.getName()+"_uc",vol.getType(),params.rows, params.cols, params.slices);
			for (int x=0; x<params.rows; x++){
				for (int y=0; y<params.cols; y++){
					for (int z=0; z<params.slices; z++){
						if(x>=params.xmin&&y>=params.ymin&&z>=params.zmin&&
								x<params.xmax&&y<params.ymax&&z<params.zmax){
							V.set(x, y, z, val=vol.getDouble(x-params.xmin, y-params.ymin, z-params.zmin)) ;
						} else {
							V.set(x,y,z,params.background);
						}
					}
				}
			}
			incrementCompletedUnits();
		}
		else{
			V = new ImageDataMipav(vol.getName()+"_uc",vol.getType(),params.rows, params.cols, params.slices,params.components);
			for (int c=0;c<params.components;c++){
				for (int x=0; x<params.rows; x++){
					for (int y=0; y<params.cols; y++){
						for (int z=0; z<params.slices; z++){
							if(x>=params.xmin&&y>=params.ymin&&z>=params.zmin&&
									x<params.xmax&&y<params.ymax&&z<params.zmax){
								V.set(x, y, z,c, vol.getDouble(x-params.xmin, y-params.ymin, z-params.zmin,c)) ;
							} else {
								V.set(x,y,z,c,params.background);
							}
						}
					}
				}
				incrementCompletedUnits();
			}
		}
//		V.setCropParameters(params);
		cropParams=params;
		V.setHeader(vol.getHeader().clone());
		markCompleted();
		return V;

		//return uncrop(vol,params.xmin,params.rows-params.xmax,params.ymin,params.cols-params.ymax,params.zmin,params.slices-params.zmax,params.background);
	}
/*
	private CubicVolume uncrop(CubicVolume f, int upRow, int downRow,
			int leftColumn, int rightColumn,
			int upSlice, int downSlice, double val) {
		System.out.println(getClass().getCanonicalName()+"\t"+"Uncropping Image...");
		String name=f.getName();
		int nx = f.getRows();
		int ny = f.getCols();
		int nz = f.getSlices();
		int nc = f.getComponents();
		System.out.println(getClass().getCanonicalName()+"\t"+"Original dimension :(" + nx + ", " + ny + ", " + nz + ")" );


	
		int mx = nx + upRow + downRow;
		int my = ny + leftColumn + rightColumn;
		int mz = nz + upSlice + downSlice;
		System.out.println(getClass().getCanonicalName()+"\t"+"New Dimensions: (" + mx + ", " + my + ", " + mz + ")");        
		CubicVolume V = null;
		if (nc ==1){
			V=f.mimic(mx,my,mz,1);
			for (int x=0; x<mx; x++){
				for (int y=0; y<my; y++){
					for (int z=0; z<mz; z++){
						if ((upRow-1 < x)&&(x<nx+upRow)&&(leftColumn-1 < y)&&(y<ny+leftColumn)
								&&(upSlice-1 < z)&&(z<nz+upSlice)){
							V.set(x,y,z,f.get(x-upRow,y-leftColumn,z-upSlice));
						}
						else{
							V.set(x,y,z,val);
						}
					}
				}
			}
		}
		else {
			V = f.mimic(mx, my, mz, nc);
			for (int l=0; l<nc; l++){
				for (int x=0; x<mx; x++){
					for (int y=0; y<my; y++){
						for (int z=0; z<mz; z++){
							if ((upRow-1 < x)&&(x<nx+upRow)&&(leftColumn-1 < y)&&(y<ny+leftColumn)
									&&(upSlice-1 < z)&&(z<nz+upSlice)){
								V.set(x,y,z,l,f.get(x-upRow,y-leftColumn,z-upSlice,l));
							}
							else{
								V.set(x,y,z,l,val);
							}
						}
					}
				} 		
			}
		}
		V.setName(name);
		return V;	

	}
*/
}
