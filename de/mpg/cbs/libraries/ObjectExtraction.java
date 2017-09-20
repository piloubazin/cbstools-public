package de.mpg.cbs.libraries;


/**
 *
 *  This class computes labels and properties for sets of binary objects
 *	
 *  @version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ObjectExtraction {
	
	// no data: used as a library of functions
	private final static int MaxObject = 1000000;
	
	public final static int SUPERIOR =  10;
	public final static int SUPEQUAL =  11;
	public final static int INFERIOR =  12;
	public final static int INFEQUAL =  13;
	public final static int EQUAL =  	14;
	public final static int UNEQUAL =  	15;
	public final static int NONE =      25;
	public final static int AND =       26;
	public final static int OR =        27;
	public final static int XOR =       28;
	
	public final static float SQR2 = (float)Math.sqrt(2.0f);
	public final static float SQR3 = (float)Math.sqrt(3.0f);
   // object manipulation
    
    public static final boolean[][][] objectFromImage(float[][][] img, int nx, int ny, int nz, float level1, float level2, int type1, int type2) {
        int x,y,z;
        boolean[][][] obj = new boolean[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            float   tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[][][] objectFromImage(float[][][] img, int nx, int ny, int nz, float level, int type) {
        int x,y,z;
        boolean[][][] obj = new boolean[nx][ny][nz];
        
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = true;
            if ( (type==INFERIOR) && (img[x][y][z] >=level) ) 	obj[x][y][z] = false;
            if ( (type==INFEQUAL) && (img[x][y][z] > level) ) 	obj[x][y][z] = false;
            if ( (type==SUPERIOR) && (img[x][y][z] <=level) ) 	obj[x][y][z] = false;
            if ( (type==SUPEQUAL) && (img[x][y][z] < level) ) 	obj[x][y][z] = false;
			if ( (type==EQUAL) 	  && (img[x][y][z]!=level) ) 	obj[x][y][z] = false;
            if ( (type==UNEQUAL)  && (img[x][y][z]==level) ) 	obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[] objectFromImage(float[] img, int nx, int ny, int nz, float level, int type) {
        boolean[] obj = new boolean[nx*ny*nz];
        
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            obj[xyz] = true;
            if ( (type==INFERIOR) && (img[xyz] >=level) ) 	obj[xyz] = false;
            if ( (type==INFEQUAL) && (img[xyz] > level) ) 	obj[xyz] = false;
            if ( (type==SUPERIOR) && (img[xyz] <=level) ) 	obj[xyz] = false;
            if ( (type==SUPEQUAL) && (img[xyz] < level) ) 	obj[xyz] = false;
            if ( (type==EQUAL) 	  && (img[xyz]!=level) ) 	obj[xyz] = false;
            if ( (type==UNEQUAL)  && (img[xyz]==level) ) 	obj[xyz] = false;
        }
        return obj;
    }
    public static final float[] floatObjectFromImage(float[] img, int nx, int ny, int nz, float level, int type) {
        float[] obj = new float[nx*ny*nz];
        
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            obj[xyz] = 1.0f;
            if ( (type==INFERIOR) && (img[xyz] >=level) ) 	obj[xyz] = 0.0f;
            if ( (type==INFEQUAL) && (img[xyz] > level) ) 	obj[xyz] = 0.0f;
            if ( (type==SUPERIOR) && (img[xyz] <=level) ) 	obj[xyz] = 0.0f;
            if ( (type==SUPEQUAL) && (img[xyz] < level) ) 	obj[xyz] = 0.0f;
            if ( (type==EQUAL) 	  && (img[xyz]!=level) ) 	obj[xyz] = 0.0f;
            if ( (type==UNEQUAL)  && (img[xyz]==level) ) 	obj[xyz] = 0.0f;
        }
        return obj;
    }
    public static enum Comparator {SUP,SUP_EQ,INF,INF_EQ};
    public static final boolean[][][] objectFromImage(float[][][] img,float thresh,Comparator comp) {
        int x,y,z;
        int nx=img.length;
        int ny=img[0].length;
        int nz=img[0][0].length;
        boolean[][][] obj = new boolean[nx][ny][nz];
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	switch(comp){
        		case SUP:obj[x][y][z]=(img[x][y][z]>thresh);break;
        		case SUP_EQ:obj[x][y][z]=(img[x][y][z]>=thresh);break;
        		case INF:obj[x][y][z]=(img[x][y][z]<thresh);break;
        		case INF_EQ:obj[x][y][z]=(img[x][y][z]<=thresh);break;
        	}
        }
        return obj;
    }
    public static final boolean[][][] objectFromLabelImage(int[][][] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        boolean[][][] obj = new boolean[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = false;
            if ( (type1==EQUAL)    && (img[x][y][z]!= level1) ) obj[x][y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = false;
            if ( (type2==EQUAL)    && (img[x][y][z]!= level2) ) obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[][][] objectFromLabelImage(int[][][][] img, int nx, int ny, int nz, int nc, int level1, int level2, int type1, int type2) {
        int x,y,z,c;
        boolean[][][] obj = new boolean[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (c=0;c<nc;c++){
            obj[x][y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z][c] >=level1) ) obj[x][y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z][c] > level1) ) obj[x][y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z][c] <=level1) ) obj[x][y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z][c] < level1) ) obj[x][y][z] = false;
            if ( (type1==EQUAL)    && (img[x][y][z][c]!= level1) ) obj[x][y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z][c] <=level2) ) obj[x][y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z][c] < level2) ) obj[x][y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z][c] >=level2) ) obj[x][y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z][c] > level2) ) obj[x][y][z] = false;
            if ( (type2==EQUAL)    && (img[x][y][z][c]!= level2) ) obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[][][] objectFromLabelImage(byte[][][] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        boolean[][][] obj = new boolean[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = false;
            if ( (type1==EQUAL)    && (img[x][y][z]!= level1) ) obj[x][y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = false;
            if ( (type2==EQUAL)    && (img[x][y][z]!= level2) ) obj[x][y][z] = false;
        }
        return obj;
    }
    public static final boolean[][][] objectFromLabelImage(byte[][][] img, int nx, int ny, int nz, int level, int type) {
		return objectFromLabelImage(img, nx,ny,nz, level, level, type, type);
	}
    public static final boolean[][][] objectFromLabelImage(int[][][] img, int nx, int ny, int nz, int level, int type) {
		return objectFromLabelImage(img, nx,ny,nz, level, level, type, type);
	}
	public static final boolean[][] objectFromLabelImage(int[][] img, int nx, int ny, int level1, int level2, int type1, int type2) {
        int x,y;
        boolean[][] obj = new boolean[nx][ny];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            obj[x][y] = true;
            if ( (type1==INFERIOR) && (img[x][y] >=level1) ) obj[x][y] = false;
            if ( (type1==INFEQUAL) && (img[x][y] > level1) ) obj[x][y] = false;
            if ( (type1==SUPERIOR) && (img[x][y] <=level1) ) obj[x][y] = false;
            if ( (type1==SUPEQUAL) && (img[x][y] < level1) ) obj[x][y] = false;
            if ( (type1==EQUAL)    && (img[x][y]!= level1) ) obj[x][y] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y] <=level2) ) obj[x][y] = false;
            if ( (type2==SUPEQUAL) && (img[x][y] < level2) ) obj[x][y] = false;
            if ( (type2==INFERIOR) && (img[x][y] >=level2) ) obj[x][y] = false;
            if ( (type2==INFEQUAL) && (img[x][y] > level2) ) obj[x][y] = false;
            if ( (type2==EQUAL)    && (img[x][y]!= level2) ) obj[x][y] = false;
        }
        return obj;
    }
	public static final byte[][][] objectFromLabelImageToByte(byte[][][] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        byte[][][] obj = new byte[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[x][y][z] = 1;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y][z] = 0;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y][z] = 0;
            if ( (type1==EQUAL)    && (img[x][y][z]!= level1) ) obj[x][y][z] = 0;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y][z] = 0;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y][z] = 0;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y][z] = 0;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y][z] = 0;
            if ( (type2==EQUAL)    && (img[x][y][z]!= level2) ) obj[x][y][z] = 0;
        }
        return obj;
    }
    public static final byte[][][] objectFromLabelImageToByte(byte[] img, int nx, int ny, int nz, int level1, int level2, int type1, int type2) {
        int x,y,z;
        byte[][][] obj = new byte[nx][ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            int   	tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        int xyz =0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	
            obj[x][y][z] = 1;
            if ( (type1==INFERIOR) && (img[xyz] >=level1) ) obj[x][y][z] = 0;
            if ( (type1==INFEQUAL) && (img[xyz] > level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPERIOR) && (img[xyz] <=level1) ) obj[x][y][z] = 0;
            if ( (type1==SUPEQUAL) && (img[xyz] < level1) ) obj[x][y][z] = 0;
            if ( (type1==EQUAL)    && (img[xyz]!= level1) ) obj[x][y][z] = 0;
            
            if ( (type2==SUPERIOR) && (img[xyz] <=level2) ) obj[x][y][z] = 0;
            if ( (type2==SUPEQUAL) && (img[xyz] < level2) ) obj[x][y][z] = 0;
            if ( (type2==INFERIOR) && (img[xyz] >=level2) ) obj[x][y][z] = 0;
            if ( (type2==INFEQUAL) && (img[xyz] > level2) ) obj[x][y][z] = 0;
            if ( (type2==EQUAL)    && (img[xyz]!= level2) ) obj[x][y][z] = 0;
            
            xyz++;
        }
        return obj;
    }
    public static final byte[][][] objectFromLabelImageToByte(byte[] img, int nx, int ny, int nz, int level, int type) {
    	return objectFromLabelImageToByte(img,nx,ny,nz,level,level,type,type);
    }
    public static final byte[][][] objectFromLabelImageToByte(byte[][][] img, int nx, int ny, int nz, int level, int type) {
    	return objectFromLabelImageToByte(img,nx,ny,nz,level,level,type,type);
    }
    public static final float[][][] labelFromImage(float[][][] img, int nx, int ny, int nz, float level1, float level2, int type1, int type2) {
		int x,y,z;
		float[][][] label = new float[nx][ny][nz];
		boolean[][][] obj = objectFromImage(img, nx, ny, nz, level1, level2, type1, type2);
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++)
			if (obj[x][y][z]) label[x][y][z] = 1.0f;
			else label[x][y][z] = 0.0f;
		obj = null;
		return label;
	}
    public static final boolean[][] objectFromImageXSlice(float[][][] img, int nx, int ny, int nz, int x, float level1, float level2, int type1, int type2) {
        int y,z;
        boolean[][] obj = new boolean[ny][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            float   tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            obj[y][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[y][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[y][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[y][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[y][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[y][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[y][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[y][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[y][z] = false;
        }
        return obj;
    }
    
    public static final boolean[][] objectFromImageYSlice(float[][][] img, int nx, int ny, int nz, int y, float level1, float level2, int type1, int type2) {
        int x,z;
        boolean[][] obj = new boolean[nx][nz];
        // make sure level1 < level2
        if (level1 > level2) {
            float   tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (z=0;z<nz;z++) {
            obj[x][z] = true;
            if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][z] = false;
            if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][z] = false;
            if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][z] = false;
            if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][z] = false;
            
            if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][z] = false;
            if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][z] = false;
            if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][z] = false;
            if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][z] = false;
        }
        return obj;
    }
    
    public static final boolean[][] objectFromImageZSlice(float[][][] img, int nx, int ny, int nz, int z, float level1, float level2, int type1, int andor, int type2) {
        int x,y;
        boolean[][] obj = new boolean[nx][ny];
        // make sure level1 < level2
        if (level1 > level2) {
            float   tmp = level1; level1 = level2; level2 = tmp;
            int     typ = type1; type1 = type2; type2 = typ; 
        }
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            if (andor==AND) {
                // and: if either of the properties is false, then no object
                obj[x][y] = true;
                if ( (type1==INFERIOR) && (img[x][y][z] >=level1) ) obj[x][y] = false;
                if ( (type1==INFEQUAL) && (img[x][y][z] > level1) ) obj[x][y] = false;
                if ( (type1==SUPERIOR) && (img[x][y][z] <=level1) ) obj[x][y] = false;
                if ( (type1==SUPEQUAL) && (img[x][y][z] < level1) ) obj[x][y] = false;

                if ( (type2==SUPERIOR) && (img[x][y][z] <=level2) ) obj[x][y] = false;
                if ( (type2==SUPEQUAL) && (img[x][y][z] < level2) ) obj[x][y] = false;
                if ( (type2==INFERIOR) && (img[x][y][z] >=level2) ) obj[x][y] = false;
                if ( (type2==INFEQUAL) && (img[x][y][z] > level2) ) obj[x][y] = false;
            } else {
                // or: if either of the properties is true, then object
                obj[x][y] = false;
                if ( (type1==INFERIOR) && (img[x][y][z] < level1) ) obj[x][y] = true;
                if ( (type1==INFEQUAL) && (img[x][y][z] <=level1) ) obj[x][y] = true;
                if ( (type1==SUPERIOR) && (img[x][y][z] > level1) ) obj[x][y] = true;
                if ( (type1==SUPEQUAL) && (img[x][y][z] >=level1) ) obj[x][y] = true;

                if ( (type2==SUPERIOR) && (img[x][y][z] > level2) ) obj[x][y] = true;
                if ( (type2==SUPEQUAL) && (img[x][y][z] >=level2) ) obj[x][y] = true;
                if ( (type2==INFERIOR) && (img[x][y][z] < level2) ) obj[x][y] = true;
                if ( (type2==INFEQUAL) && (img[x][y][z] <=level2) ) obj[x][y] = true;
            }
        }
        return obj;
    }
	public static final boolean[] objectFromLabelImage(byte[] img, int nx, int ny, int nz, byte level, int type) {
        boolean[] obj = new boolean[nx*ny*nz];
        
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            obj[xyz] = true;
            if ( (type==INFERIOR) && (img[xyz] >=level) ) 	obj[xyz] = false;
            if ( (type==INFEQUAL) && (img[xyz] > level) ) 	obj[xyz] = false;
            if ( (type==SUPERIOR) && (img[xyz] <=level) ) 	obj[xyz] = false;
            if ( (type==SUPEQUAL) && (img[xyz] < level) ) 	obj[xyz] = false;
            if ( (type==EQUAL) 	  && (img[xyz]!=level) ) 	obj[xyz] = false;
            if ( (type==UNEQUAL)  && (img[xyz]==level) ) 	obj[xyz] = false;
        }
        return obj;
    }
	public static final boolean[] objectFromLabelImage(int[] img, int nx, int ny, int nz, int level, int type) {
        boolean[] obj = new boolean[nx*ny*nz];
        
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            obj[xyz] = true;
            if ( (type==INFERIOR) && (img[xyz] >=level) ) 	obj[xyz] = false;
            if ( (type==INFEQUAL) && (img[xyz] > level) ) 	obj[xyz] = false;
            if ( (type==SUPERIOR) && (img[xyz] <=level) ) 	obj[xyz] = false;
            if ( (type==SUPEQUAL) && (img[xyz] < level) ) 	obj[xyz] = false;
            if ( (type==EQUAL) 	  && (img[xyz]!=level) ) 	obj[xyz] = false;
            if ( (type==UNEQUAL)  && (img[xyz]==level) ) 	obj[xyz] = false;
        }
        return obj;
    }

}
