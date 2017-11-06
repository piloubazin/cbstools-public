package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;


/*
 *
 *  This class computes various morphological operations on 2D and 3D images
 *	
 *	@version    October 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ObjectMorphology {
	
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

	/** erode binary object with a square kernel */
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
						break;
					}
				}
			}
		}
        return eroded;
    }
    
	/** erode binary object with a square kernel */
	public static boolean[] erodeObject(boolean[] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        boolean[] eroded = new boolean[nx*ny*nz];
		
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			
			eroded[x+nx*y+nx*ny*z] = true;
			for (i=-dx;i<=dx;i++) for (j=-dy;j<=dy;j++) for (k=-dz;k<=dz;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[(x+i)+nx*(y+j)+nx*ny*(z+k)]==false) { 
						eroded[x+nx*y+nx*ny*z] = false;
						break;
					}
				}
			}
		}
        return eroded;
    }
    
	public static void fastErodeObject(boolean[] img, int nx, int ny, int nz, int d) {
        int xyz, x, y, z, t;
		boolean[] eroded = new boolean[nx*ny*nz];
		
		for (t=0;t<d;t++) {
			for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
				xyz = x+nx*y+nx*ny*z;
				eroded[xyz] = img[xyz];
				if (eroded[xyz]) {
					if (!img[xyz+1]) eroded[xyz] = false;
					else if (!img[xyz-1]) eroded[xyz] = false;
					else if (!img[xyz+nx]) eroded[xyz] = false;
					else if (!img[xyz-nx]) eroded[xyz] = false;
					else if (!img[xyz+nx*ny]) eroded[xyz] = false;
					else if (!img[xyz-nx*ny]) eroded[xyz] = false;
				}
			}
			for (xyz=0;xyz<nx*ny*nz;xyz++) {
				img[xyz] = eroded[xyz];
			}
		}
        return;
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

	public static boolean[] dilateObject(boolean[] img, int nx, int ny, int nz, int dx, int dy, int dz) {
        int x,y,z;
		int i,j,k;
        boolean[] dilated = new boolean[nx*ny*nz];
		boolean search;
        
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) if (dilated[x+nx*y+nx*ny*z]==false) {
			search = true;
			for (i=-dx;i<=dx && search;i++) for (j=-dy;j<=dy && search;j++) for (k=-dz;k<=dz && search;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i+nx*(y+j)+nx*ny*(z+k)]==true) {
						dilated[x+nx*y+nx*ny*z] = true;
						search = false;
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
		boolean search;
        
		// dx,dy,dz describe the structuring element ( x+/-dx, y+/-dy, z+/-dz )

        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) if (dilated[x+nx*y+nx*ny*z]==0) {
			search = true;
			for (i=-dx;i<=dx && search;i++) for (j=-dy;j<=dy && search;j++) for (k=-dz;k<=dz && search;k++) {
				
				if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
			
					if (img[x+i+nx*(y+j)+nx*ny*(z+k)]>0) {
						dilated[x+nx*y+nx*ny*z] = 1;
						search = false;
					}
				}
			}
		}
        return dilated;
    }

    
	public static void fastDilateObject(boolean[] img, int nx, int ny, int nz, int d) {
        int xyz, x, y, z, t;
		boolean[] dilated = new boolean[nx*ny*nz];
		
		for (t=0;t<d;t++) {
			for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
				xyz = x+nx*y+nx*ny*z;
				dilated[xyz] = img[xyz];
				if (!dilated[xyz]) {
					if (img[xyz+1]) dilated[xyz] = true;
					else if (img[xyz-1]) dilated[xyz] = true;
					else if (img[xyz+nx]) dilated[xyz] = true;
					else if (img[xyz-nx]) dilated[xyz] = true;
					else if (img[xyz+nx*ny]) dilated[xyz] = true;
					else if (img[xyz-nx*ny]) dilated[xyz] = true;
				}
			}
			for (xyz=0;xyz<nx*ny*nz;xyz++) {
				img[xyz] = dilated[xyz];
			}
		}
        return;
    }

	public static void fastDilateFunction(float[][][] img, int nx, int ny, int nz, int d, float ratio) {
        int x, y, z, t;
		float[][][] dilated = new float[nx][ny][nz];
		
		for (t=0;t<d;t++) {
			for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
				dilated[x][y][z] = img[x][y][z];
				
				if (ratio*img[x+1][y][z]>dilated[x][y][z]) dilated[x][y][z] = ratio*img[x+1][y][z];
				if (ratio*img[x-1][y][z]>dilated[x][y][z]) dilated[x][y][z] = ratio*img[x-1][y][z];
				if (ratio*img[x][y+1][z]>dilated[x][y][z]) dilated[x][y][z] = ratio*img[x][y+1][z];
				if (ratio*img[x][y-1][z]>dilated[x][y][z]) dilated[x][y][z] = ratio*img[x][y-1][z];
				if (ratio*img[x][y][z+1]>dilated[x][y][z]) dilated[x][y][z] = ratio*img[x][y][z+1];
				if (ratio*img[x][y][z-1]>dilated[x][y][z]) dilated[x][y][z] = ratio*img[x][y][z-1];
				
			}
			for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
				img[x][y][z] = dilated[x][y][z];
			}
		}
        return;
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
