package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;

import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This algorithm dilates the data on a levelset surface in the direction normal to the surface to a defined distance
 *
 *	@version    Feb 2013
 *	@author     Christine Lucas Tardif 
 *		
 *
 */
 
public class CorticalSurfaceDilation {
	
	private static int connectivity;
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	// data and membership buffers
	private 	float[]			dilData;
    private		float[]			data;
    private		boolean[]		datamask;
    private 	float[] 		lvlset;  			// level set functions
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	
	// new constructor
	public CorticalSurfaceDilation(float[] data_, boolean[] datamask_, float[] lvlset_, boolean[] mask_, 
			int nx_, int ny_, int nz_,
			float rx_, float ry_, float rz_,
			String connectivityType_) {
			
		data = data_;
		datamask = datamask_;
		lvlset = lvlset_;
		mask = mask_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		if (datamask == null) {
			datamask = new boolean[nx*ny*nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				datamask[xyz] = true;
			}
		}
		
		dilData = new float[nx*ny*nz];
		
		if (connectivityType_.equals("6")) {
			connectivity = 6;
			xoff = new int[]{1, -1, 0, 0, 0, 0};
			yoff = new int[]{0, 0, nx, -nx, 0, 0};
			zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		} else if (connectivityType_.equals("18")) {
			connectivity = 18;
			xoff = new int[]{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0};
			yoff = new int[]{0, 0, nx, -nx, 0, 0, nx, -nx, -nx, nx, 0, 0, 0, 0, nx, -nx, nx, -nx};
			zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny, 0, 0, 0, 0, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny,-nx*ny, nx*ny};
		} else if (connectivityType_.equals("26")) {
			connectivity = 26;
			xoff = new int[]{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};
			yoff = new int[]{0, 0, nx, -nx, 0, 0, nx, -nx, -nx, nx, 0, 0, 0, 0, nx, -nx, nx, -nx, nx, nx, -nx, -nx, nx, nx, -nx, -nx};
			zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny, 0, 0, 0, 0, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny};
		} else {
			System.out.println("Using default 6-neighbour connectivity");
			connectivity = 6;
			xoff = new int[]{1, -1, 0, 0, 0, 0};
			yoff = new int[]{0, 0, nx, -nx, 0, 0};
			zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		}
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[x+nx*y+nx*ny*z] = false;
		}
		
		
	}
	
		
	public void finalize() {
		data = null;
		mask = null;
		lvlset = null;
		
	}

	
	public final float[][][] exportDilatedData() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z] = dilData[xyz];
		}
		return res;
	}
	
	public final float[] getDilatedData() {
		return dilData;
	}
	
	public void dilateData() {
       
		BitSet boundary = new BitSet(nx*ny*nz);
        BitSet nextboundary = new BitSet(nx*ny*nz);
        BitSet used = new BitSet(nx*ny*nz);
        BitSet nextused = new BitSet(nx*ny*nz);
        
        System.out.println("\n Set initial processed data");
        //set initial used
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;	
			if (FastMath.abs(lvlset[xyz]) <= (float)FastMath.sqrt(3.0)/2.0f && datamask[xyz]) {
				used.set(xyz);
				dilData[xyz] = data[xyz];
			}
        }
       	
        System.out.println("\n Set initial boundary");
        //set initial boundary
        for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;	
			
			if (used.get(xyz)) { 
				for (int k=0; k<connectivity; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if ((mask[xyzn]) && (!used.get(xyzn))&& (!boundary.get(xyzn)) ) {
						boundary.set(xyzn);
					}
				}
			}
        }
        
		//dilate boundary
		boolean dilate = true;
        float sum = 0.0f;
        int count = 0;
		
        nextused = (BitSet)used.clone();
        
    	System.out.println("\n Dilate");
               
        while (dilate) {
        	System.out.println(".");
        	dilate = false;
        	
			for (int xyz = boundary.nextSetBit(0); xyz >= 0; xyz = boundary.nextSetBit(xyz+1)) {
				sum = 0.0f;
				count = 0;
				
				// dilate and build the next boundary
				for (int k=0; k<connectivity; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
										
					if (used.get(xyzn)) {
						sum += dilData[xyzn];
						count++;
					} else {
						if ((mask[xyzn]) && (!boundary.get(xyzn))) { 
							nextboundary.set(xyzn);
							dilate = true;
						}
					}
				}
				
				dilData[xyz] = sum/count;
			
				// add to the processed list
				nextused.set(xyz);
			}
			
			// replace active boundary
			boundary = (BitSet) nextboundary.clone();
			nextboundary.clear();

			// replace active processed list
			used = (BitSet) nextused.clone();		
        }
		
    }
	
	
	
}

