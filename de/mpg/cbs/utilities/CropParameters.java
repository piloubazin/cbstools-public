package de.mpg.cbs.utilities;

import javax.vecmath.Point3i;

// TODO: Auto-generated Javadoc
/**
 * The Class CropParameters.
 */
public class CropParameters {
	
	/** The xmin. */
	public int xmin;
    
    /** The xmax. */
    public int xmax;
    
    /** The ymin. */
    public int ymin;
    
    /** The ymax. */
    public int ymax;
    
    /** The zmin. */
    public int zmin;
    
    /** The zmax. */
    public int zmax;
    
    /** The rows. */
    public int rows;
    
    /** The cols. */
    public int cols;
    
    /** The slices. */
    public int slices;
    
    /** The components. */
    public int components;
    
    /** The background. */
    public double background;
	
	/**
	 * Instantiates a new crop parameters.
	 * 
	 * @param rows the rows
	 * @param cols the cols
	 * @param slices the slices
	 * @param components the components
	 * @param xmin the xmin
	 * @param xmax the xmax
	 * @param ymin the ymin
	 * @param ymax the ymax
	 * @param zmin the zmin
	 * @param zmax the zmax
	 * @param background the background
	 */
	public CropParameters(int rows,int cols,int slices,int components,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax,double background){
		this.rows=rows;
		this.cols=cols;
		this.slices=slices;
		this.components=components;
		this.xmin=xmin;
		this.xmax=xmax;
		this.ymin=ymin;
		this.ymax=ymax;
		this.zmin=zmin;
		this.zmax=zmax;
		this.background=background;
	}
	
	public CropParameters(Point3i offset, Point3i dims){
		rows=dims.x-offset.x-1;
		cols=dims.y-offset.y-1;
		slices=dims.z-offset.z-1;
		xmin = offset.x;
		ymin = offset.y;
		zmin = offset.z;
		xmax = dims.x+offset.x-1;
		ymax = dims.y+offset.y-1;
		zmax = dims.z+offset.z-1;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString(){
		return String.format("Corners: Lower (%d,%d,%d) Upper (%d,%d,%d)\nDimensions: Uncropped (%d,%d,%d) Cropped (%d,%d,%d)", xmin,ymin,zmin,xmax,ymax,zmax,rows,cols,slices,getCroppedRows(),getCroppedCols(),getCroppedSlices());
	}
	
	/**
	 * Gets the cropped rows.
	 * 
	 * @return the cropped rows
	 */
	public int getCroppedRows(){
		return xmax-xmin+1;
	}
	
	/**
	 * Gets the cropped cols.
	 * 
	 * @return the cropped cols
	 */
	public int getCroppedCols(){
		return ymax-ymin+1;
	}
	
	/**
	 * Gets the cropped slices.
	 * 
	 * @return the cropped slices
	 */
	public int getCroppedSlices(){
		return zmax-zmin+1;
	}
}
