package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

import de.mpg.cbs.utilities.*;

/*
 *
 *  This class computes various morphological operations on 2D and 3D images
 *	
 *	@version    October 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class Morphology {
	
	// no data: used as a library of functions
	
	/*
	* 	erode image with a square kernel 
	*	scalar erosion: eroded = min_kernel (img)
	*/
    public static float[][][] erodeImage(float[][][] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        float[][][] eroded = new float[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			
			eroded[x][y][z] = img[x][y][z];
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i][y+j][z+k] < eroded[x][y][z]) eroded[x][y][z] = img[x+i][y+j][z+k];
				}
			}
		}
        return eroded;
    }
    
    /*
	* 	dilate image with a square kernel 
	*	scalar dilation: dilated = max_kernel (img)
	*/
    public static float[][][] dilateImage(float[][][] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        float[][][] dilated = new float[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			
			dilated[x][y][z] = img[x][y][z];
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i][y+j][z+k] > dilated[x][y][z]) dilated[x][y][z] = img[x+i][y+j][z+k];
				}
			}
		}
        return dilated;
    }

    /*
	* 	erode image with a square kernel 
	*	absolute scalar dilation: eroded = minabs_kernel (img)
	*/
   public static float[][][] erodeAbsImage(float[][][] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        float[][][] eroded = new float[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			
			eroded[x][y][z] = img[x][y][z];
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (Numerics.abs(img[x+i][y+j][z+k])<Numerics.abs(eroded[x][y][z])) { 
						eroded[x][y][z] = img[x+i][y+j][z+k];
					}
				}
			}
		}
        return eroded;
    }
    

	/** erode binary object with a custom kernel */
	public static boolean[][][] erodeObject(boolean[][][] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        boolean[][][] eroded = new boolean[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			eroded[x][y][z] = true;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
					if (img[x+i][y+j][z+k]==false) { 
						eroded[x][y][z] = false;
					}
				}
			}
		}
        return eroded;
    }
    
	/** erode binary object with a custom kernel */
	public static byte[][][] erodeObject(byte[][][] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        byte[][][] eroded = new byte[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			eroded[x][y][z] = 1;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
					if (img[x+i][y+j][z+k]==0) { 
						eroded[x][y][z] = 0;
					}
				}
			}
		}
        return eroded;
    }
    
	/** erode binary object with a custom kernel */
	public static boolean[] erodeObject(boolean[] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        boolean[] eroded = new boolean[nx*ny*nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	int xyz = x+nx*y+nx*ny*z;
			eroded[xyz] = true;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
					if (img[x+i+nx*(y+j)+nx*ny*(z+k)]==false) { 
						eroded[xyz] = false;
						break;
					}
				}
			}
		}
        return eroded;
    }
    
	/** erode binary object with a custom kernel */
	public static byte[] erodeObject(byte[] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        byte[] eroded = new byte[nx*ny*nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	int xyz = x+nx*y+nx*ny*z;
			eroded[xyz] = 1;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
					if (img[x+i+nx*(y+j)+nx*ny*(z+k)]==0) { 
						eroded[xyz] = 0;
						break;
					}
				}
			}
		}
        return eroded;
    }
    
	/** erode binary object with a custom kernel */
	public static int[] erodeObject(int[] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        int[] eroded = new int[nx*ny*nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	int xyz = x+nx*y+nx*ny*z;
			eroded[xyz] = 1;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
					if (img[x+i+nx*(y+j)+nx*ny*(z+k)]==0) { 
						eroded[xyz] = 0;
						break;
					}
				}
			}
		}
        return eroded;
    }
    
    /** dilate binary object with a square kernel */
	public static boolean[][][] dilateObject(boolean[][][] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        boolean[][][] dilated = new boolean[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			
			dilated[x][y][z] = false;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i][y+j][z+k]==true) {
						dilated[x][y][z] = true;
						break;
					}
				}
			}
		}
        return dilated;
    }

 	public static boolean[] dilateObject(boolean[] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        boolean[] dilated = new boolean[nx*ny*nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			dilated[xyz] = false;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i + nx*(y+j) + nx*ny*(z+k)]==true) {
						dilated[xyz] = true;
						break;
					}
				}
			}
		}
        return dilated;
    }

 	public static byte[] dilateObject(byte[] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        byte[] dilated = new byte[nx*ny*nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			dilated[xyz] = 0;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i + nx*(y+j) + nx*ny*(z+k)]>0) {
						dilated[xyz] = 1;
						break;
					}
				}
			}
		}
        return dilated;
    }

 	public static int[] dilateObject(int[] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        int[] dilated = new int[nx*ny*nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			dilated[xyz] = 0;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i + nx*(y+j) + nx*ny*(z+k)]>0) {
						dilated[xyz] = 1;
						break;
					}
				}
			}
		}
        return dilated;
    }

   /** dilate binary object with a square kernel */
	public static byte[][][] dilateObject(byte[][][] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        byte[][][] dilated = new byte[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			
			dilated[x][y][z] = 0;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i][y+j][z+k]>0) {
						dilated[x][y][z] = 1;
						break;
					}
				}
			}
		}
        return dilated;
    }

	/** erode binary object with a custom kernel */
	public static boolean[][][] erodeObject(boolean[][][] img, int nx, int ny, int nz, boolean[][][] mask, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        boolean[][][] eroded = new boolean[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			eroded[x][y][z] = true;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
					if ( (mask[i+dx][j+dy][k+dz]) && (img[x+i][y+j][z+k]==false)) { 
						eroded[x][y][z] = false;
					}
				}
			}
		}
        return eroded;
    }
    
    /** dilate binary object with a custom kernel */
	public static boolean[][][] dilateObject(boolean[][][] img, int nx, int ny, int nz, boolean[][][] mask, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        boolean[][][] dilated = new boolean[nx][ny][nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			dilated[x][y][z] = false;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {			
					if ( (mask[i+dx][j+dy][k+dz]) && (img[x+i][y+j][z+k]==true) ) {
						dilated[x][y][z] = true;
					}
				}
			}
		}
        return dilated;
    }

	/*
	* 	erode binary object with a custom kernel
	*	using the BitSet structure with indexing convention
	*	index = x + nx*y + nx*ny*z 
	*/
	public static BitSet erodeObject(BitSet img, int nx, int ny, int nz, BitSet mask, int dx, int dy, int dz) {
        BitSet eroded = new BitSet(nx*ny*nz);
		int ndx = 2*dx+1;
		int ndxndy = ndx*(2*dy+1);
		int nxny = nx*ny;
		int x,y,z;
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )
		for(int index=img.nextSetBit(0); index>=0; index=img.nextSetBit(index+1)) {
			z = index/nxny;
			y = (index%nxny)/nx;
			x = ((index%nxny)%nx);
			eroded.set( index, true );
			for (int i=-dx;i<=dx;i++) for (int j=-dy;j<=dy;j++) for (int k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
					if ( (mask.get( i+dx + ndx*(j+dy) + ndxndy*(k+dz)) ) 
						&& (!img.get( x+i + nx*(y+j) + nxny*(z+k)) )) { 
							eroded.set( index, false );
					}
				}
			}
		}
        return eroded;
    }
    
    /*
	* 	dilate binary object with a custom kernel
	*	using the BitSet structure with indexing convention
	*	index = x + nx*y + nx*ny*z 
	*/
	public static BitSet dilateObject(BitSet img, int nx, int ny, int nz, BitSet mask, int dx, int dy, int dz) {
        BitSet dilated = new BitSet(nx*ny*nz);
		int ndx = 2*dx+1;
		int ndxndy = ndx*(2*dy+1);
		int nxny = nx*ny;
		int x,y,z;
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )
		for(int index=img.nextSetBit(0); index>=0; index=img.nextSetBit(index+1)) {
			z = index/nxny;
			y = (index%nxny)/nx;
			x = ((index%nxny)%nx);
			for (int i=-dx;i<=dx;i++) for (int j=-dy;j<=dy;j++) for (int k=-dz;k<=dz;k++) {
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {		
					if ( (mask.get( i+dx + ndx*(j+dy) + ndxndy*(k+dz)) )) { 
						dilated.set( x+i + nx*(y+j) + nxny*(z+k) );
					}
				}
			}
		}
        return dilated;
    }
}
