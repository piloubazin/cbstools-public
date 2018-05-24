package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This class computes labels and properties for sets of binary objects
 *	
 *  @version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ObjectLabeling {
	
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

    public static final int countLabels(int[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb, label[x][y][z]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final byte countLabels(byte[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        byte Nlb;
        ArrayList<Byte> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb, label[x][y][z]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final byte countLabels(byte[] label, int nx, int ny, int nz) {
        int xyz,n;
        byte Nlb;
        ArrayList<Byte> lb = new ArrayList();
        boolean newLabel;
        
        Nlb = 0;
        for (xyz=0;xyz<nx*ny*nz;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb, label[xyz]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int countLabels(int[] label, int nx, int ny, int nz) {
        int xyz,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        Nlb = 0;
        for (xyz=0;xyz<nx*ny*nz;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb, label[xyz]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int countLabels(float[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int countLabels(float[] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
			int xyz = x+nx*y+nx*ny*z;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int countLabels(int[][] label, int nx, int ny) {
        int x,y,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
		Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y]);
                Nlb++;
            }
        }
        return Nlb;
    }
    
    public static final int[] listLabels(int[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    
    public static final int[] listLabels(int[][] label, int nx, int ny) {
        int x,y,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }

    public static final float[] listLabels(float[][] label, int nx, int ny) {
        int x,y,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y]);
                Nlb++;
            }
        }
		float[] list = new float[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }

    public static final float[] listLabels(float[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		float[] list = new float[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }

    public static final float[] listLabels(float[] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
		float[] list = new float[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    public static final int[] listLabels(int[] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        int ind = -1;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	ind = x + y*nx + z*nx*ny;
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[ind]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[ind]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    public static final int[] listLabels(int[] label, int nx, int ny, int nz, int offset) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        int ind = -1;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	ind = offset + x + y*nx + z*nx*ny;
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[ind]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[ind]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    public static final int[] listLabels(byte[] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        int ind = -1;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        	ind = x + y*nx + z*nx*ny;
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[ind]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,(int)label[ind]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) list[n] = lb.get(n);
		
        return list;
    }
    public static final float[] listOrderedLabels(float[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Float> lb = new ArrayList();
        boolean newLabel;
        
        //lb.add(0,0.0f);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		float[] list = new float[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final int[] listOrderedLabels(int[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }
    
    public static final int[] listOrderedLabels4d(int[][][][] label, int nx, int ny, int nz, int nc) {
        int x,y,z,c,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (c=0;c<nc;c++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z][c]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y][z][c]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }
    
    public static final int[] listOrderedLabels(int[][] label, int nx, int ny) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++)  {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[x][y]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final int[] listOrderedLabels(int[] label, int nx, int ny, int nz) {
        int xyz,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (xyz=0;xyz<nx*ny*nz;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }
    
    public static final byte[] listOrderedLabels(byte[] label, int nx, int ny, int nz) {
        int xyz,n;
        int Nlb;
        ArrayList<Byte> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (xyz=0;xyz<nx*ny*nz;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
		byte[] list = new byte[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }
    
    public static final int[] listOrderedLabels(int[] label, int nx, int ny) {
        int xyz,n;
        int Nlb;
        ArrayList<Integer> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,0);
        Nlb = 0;
        for (xyz=0;xyz<nx*ny;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add(Nlb,label[xyz]);
                Nlb++;
            }
        }
		int[] list = new int[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final short[] listOrderedLabels(short[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Short> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,(short)0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add((short)Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		short[] list = new short[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final short[] listOrderedLabels(short[] label, int nx, int ny, int nz) {
        int xyz,n;
        int Nlb;
        ArrayList<Short> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,(short)0);
        Nlb = 0;
        for (xyz=0;xyz<nx*ny*nz;xyz++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[xyz]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add((short)Nlb,label[xyz]);
                Nlb++;
            }
        }
		short[] list = new short[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }

    public static final byte[] listOrderedLabels(byte[][][] label, int nx, int ny, int nz) {
        int x,y,z,n;
        int Nlb;
        ArrayList<Byte> lb = new ArrayList();
		boolean newLabel;
        
        //lb.add(0,(byte)0);
        Nlb = 0;
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            newLabel=true;
            for (n=0;n<Nlb;n++)
                if (label[x][y][z]==lb.get(n)) { newLabel=false; break; }
            if (newLabel) {
                lb.add((byte)Nlb,label[x][y][z]);
                Nlb++;
            }
        }
		byte[] list = new byte[Nlb];
		for (n=0;n<Nlb;n++) {
			list[n] = lb.get(n);
			for (int m=n+1;m<Nlb;m++) {
				if (lb.get(m)<list[n]) {
					// switch place
					list[n] = lb.get(m);
					lb.set(m, lb.get(n) );
					lb.set(n, list[n]);
				}
			}
		}
		
        return list;
    }
    
    public static final int[] invertList(int[] list) {
    	int maxlb = -1;
    	for (int n=0;n<list.length;n++) {
    		if (list[n]>maxlb) maxlb = list[n];
    	}
    	int[] lb = new int[maxlb+1];
		for (int n=0;n<lb.length;n++) lb[n] = -1;
		
    	for (int n=0;n<list.length;n++) {
    		lb[list[n]] = n;
    	}
    	return lb;
    }

    public static final boolean[][][] largestObjectFromLabel(int[][][] label, int nlb, int nx, int ny, int nz) {
        int x,y,z,n;
        boolean[][][] obj = new boolean[nx][ny][nz];
        int[] Nobj = new int[nlb];
        int best,size;
        
        for (n=0;n<nlb;n++) Nobj[n]=0;
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
            Nobj[ label[x][y][z] ]++;
        }
        size=0;best=0;
        for (n=1;n<nlb;n++) if (Nobj[n]>size) {
            size = Nobj[n];
            best = n;
        }
		if (best>0)
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
				obj[x][y][z] = (label[x][y][z]==best);
			}
		else
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
				obj[x][y][z] = false;
			}
        return obj;
    }
    
   public static final boolean[] largestObjectFromLabel(int[] label, int nlb, int nx, int ny, int nz) {
        boolean[] obj = new boolean[nx*ny*nz];
        
        int[] Nobj = new int[nlb];
        int best,size;
        
        for (int n=0;n<nlb;n++) Nobj[n]=0;
        
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
            Nobj[ label[xyz] ]++;
        }
        size=0;best=0;
        for (int n=1;n<nlb;n++) if (Nobj[n]>size) {
            size = Nobj[n];
            best = n;
        }
		if (best>0)
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				obj[xyz] = (label[xyz]==best);
			}
		else
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				obj[xyz] = false;
			}
        return obj;
    }
    
    public static final boolean[][][] largestObjectFromLabel(int[][][] label, int nx, int ny, int nz) {
		return largestObjectFromLabel(label, countLabels(label,nx,ny,nz), nx,ny,nz);
    }
    
    public static final boolean[] largestObjectFromLabel(int[] label, int nx, int ny, int nz) {
		return largestObjectFromLabel(label, countLabels(label,nx,ny,nz), nx,ny,nz);
    }
    
    public static final boolean[][] largestObjectFromLabel(int[][] label, int nlb, int nx, int ny) {
        int x,y,n;
        boolean[][] obj = new boolean[nx][ny];
        int[] Nobj = new int[nlb];
        int best,size;
        
        for (n=0;n<nlb;n++) Nobj[n]=0;
        
        for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
            Nobj[ label[x][y] ]++;
        }
        size=0;best=0;
        for (n=1;n<nlb;n++) if (Nobj[n]>size) {
            size = Nobj[n];
            best = n;
        }
		if (best>0)
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
				obj[x][y] = (label[x][y]==best);
			}
		else
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
				obj[x][y] = false;
			}
        return obj;
    }
    
    public static final boolean[][] largestObjectFromLabel(int[][] label, int nx, int ny) {
		return largestObjectFromLabel(label, countLabels(label,nx,ny), nx,ny);
    }
    
    
    public static final boolean[][][] connectedObjectAt(boolean[][][] object, int x0, int y0, int z0, int nx, int ny, int nz) {
        boolean[][][] obj = new boolean[nx][ny][nz];
        
    	int[][][] label = connected6Object3D(object, nx, ny, nz);
		 
        int lb = label[x0][y0][z0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			obj[x][y][z] = (label[x][y][z]==lb);
		}
		return obj;
    }
    
    /*
     * @brief Fill holes that are not 'conn' connected to the background.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     * @param conn Connectivity to use.
     */
    public static final boolean[][][] removeHoles(boolean[][][] object, int nx, int ny, int nz, int conn) {		
		for (int x = 0; x < nx; x++)
			for (int y = 0; y < ny; y++)
				for (int z = 0; z < nz; z++) {
					object[x][y][z] = !object[x][y][z];
				}
		
    	int[][][] lb;
		
		if (conn == 6) {
			lb = connected6Object3D(object, nx, ny, nz);
		} else if (conn == 18) {
			lb = connected18Object3D(object, nx, ny, nz);
		} else if (conn == 26) {
			lb = connected26Object3D(object, nx, ny, nz);
		} else {
			System.out.println("Unsupported connectivity: " + conn + " \n");
			return null;
		}
		
		object = largestObjectFromLabel(lb, nx, ny, nz);
		
		for (int x = 0; x < nx; x++)
			for (int y = 0; y < ny; y++)
				for (int z = 0; z < nz; z++) {
					object[x][y][z] = !object[x][y][z];
				}
		
		return object;
    }
    
    public static final boolean[] removeHoles(boolean[] object, int nx, int ny, int nz, int conn) {		
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			object[xyz] = !object[xyz];
		}
		
    	int[] lb;
		
		if (conn == 6) {
			lb = connected6Object3D(object, nx, ny, nz);
		} else if (conn == 18) {
			lb = connected18Object3D(object, nx, ny, nz);
		} else if (conn == 26) {
			lb = connected26Object3D(object, nx, ny, nz);
		} else {
			System.out.println("Unsupported connectivity: " + conn + " \n");
			return null;
		}
		
		object = largestObjectFromLabel(lb, nx, ny, nz);
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			object[xyz] = !object[xyz];
		}
		
		return object;
    }
    
    /*
     * @brief Fill holes that are not 6 connected to the background.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] removeHoles6(boolean[][][] object, int nx, int ny, int nz) {
    	return removeHoles(object, nx, ny, nz, 6);
    }
    
    
    /*
     * @brief Fill holes that are not 18 connected to the background.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] removeHoles18(boolean[][][] object, int nx, int ny, int nz) {
    	return removeHoles(object, nx, ny, nz, 18);
    }
    
    
    /*
     * @brief Returns the largest 'conn' connected object.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     * @param conn Connectivity to use.
     */
    public static final boolean[][][] largestObject(boolean[][][] object, int nx, int ny, int nz, int conn) {		
    	int[][][] lb;
		
		if (conn == 6) {
			lb = connected6Object3D(object, nx, ny, nz);
		} else if (conn == 18) {
			lb = connected18Object3D(object, nx, ny, nz);
		} else if (conn == 26) {
			lb = connected26Object3D(object, nx, ny, nz);
		} else {
			System.out.println("Unsupported connectivity: " + conn + " \n");
			return null;
		}
		
		object = largestObjectFromLabel(lb, nx, ny, nz);
		
		return object;
    }
    
    public static final boolean[] largestObject(boolean[] object, int nx, int ny, int nz, int conn) {		
    	int[] lb;
		
		if (conn == 6) {
			lb = connected6Object3D(object, nx, ny, nz);
		} else if (conn == 18) {
			lb = connected18Object3D(object, nx, ny, nz);
		} else if (conn == 26) {
			lb = connected26Object3D(object, nx, ny, nz);
		} else {
			System.out.println("Unsupported connectivity: " + conn + " \n");
			return null;
		}
		object = largestObjectFromLabel(lb, countLabels(lb, nx, ny, nz), nx, ny, nz);
		
		return object;
    }
    
   
    /*
     * @brief Returns the largest 6 connected object.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] largest6Object(boolean[][][] object, int nx, int ny, int nz) {		
    	return largestObject(object, nx, ny, nz, 6);
    }
    
    
    /*
     * @brief Returns the largest 18 connected object.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] largest18Object(boolean[][][] object, int nx, int ny, int nz) {		
    	return largestObject(object, nx, ny, nz, 18);
    }
   
    
    /*
     * @brief Returns the largest 18 connected object.
     * @param object Boolean volume to work on.
     * @param nx X dimension.
     * @param ny Y dimension.
     * @param nz Z dimension.
     */
    public static final boolean[][][] largest26Object(boolean[][][] object, int nx, int ny, int nz) {		
    	return largestObject(object, nx, ny, nz, 26);
    }
    	
	/** 
	 *	Connected components of an object.
     *  2D images: 4-connectivity
	 */
	public static final int[][] connected4Object2D(boolean img[][], int nx, int ny) {
		int Nlabel;
		int[][] label = new int[nx][ny];
		ArrayList<Integer>   lb = new ArrayList();
		int lbMin;
		int x,y,c,i,j,k;
		int Nlb;
		int[]   connect = new int[4];
		int Nconnect;
		int AddLabel=0;
	
		// the input is a 3x3 binary image (0 out, 1 in)
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				label[x][y] = 0;
			}
		}
		
		lb.add(0,0);
		Nlabel = 1;
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				if (img[x][y]) {
					// object point: neighbors ?
					Nconnect = 0;
                    for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) {
                        if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) ) {
                            if (i*i+j*j < 2) {
                                if (label[x+i][y+j] > 0) {
                                    connect[Nconnect] = lb.get( label[x+i][y+j] );
                                    Nconnect++;
                                }
                            }
                        }
                    }
					// if connected values, find the smallest lb label and attribute it
					// to all others (-> join labels)
					if (Nconnect>0) {
						lbMin = lb.get(connect[0]);
						for (k=1;k<Nconnect;k++) lbMin = Math.min(lbMin,lb.get(connect[k]) );
						for (k=0;k<Nconnect;k++) lb.set(connect[k],lbMin);
						label[x][y] = lbMin;
					} else {
						// new, unconnected region
						label[x][y] = Nlabel;
						lb.add(Nlabel,Nlabel);
						Nlabel++;
						/* not needed with ArrayLists
						// check if the number of labels is above the threshold
						if (Nlabel>=lb.length-1) {
							int[] tmp = new int[2*lb.length];
							for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
							lb = tmp;
						}
						*/
					}
				}
			}
		}
		// only one level of labels
		for (k=1;k<Nlabel;k++) {
			c = k;
			while (lb.get(c)!=c) c = lb.get(c);
			lb.set(k, c);
		}
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (k=1;k<Nlabel;k++) {
			if (lb.get(k)==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// copy on label image
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				label[x][y] = lb2[ lb.get( label[x][y] ) ];
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	
	/** 
	 *	Connected components of an object.
     *  2D images: 8-neighborhood 
	 */
	public static final int[][] connected8Object2D(boolean img[][], int nx, int ny) {
		int Nlabel;
		int[][] label = new int[nx][ny];
		int[]   lb = new int[MaxObject];
		int lbMin;
		int x,y,c,i,j,k;
		int Nlb;
		int[]   connect = new int[4];
		int Nconnect;
		int AddLabel=0;
	
		// the input is a 3x3 binary image (0 out, 1 in)
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				label[x][y] = 0;
			}
		}
		
		lb[0] = 0;
		Nlabel = 1;
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				if (img[x][y]) {
					// object point: neighbors ?
					// object point: neighbors ?
					Nconnect = 0;
                    for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) {
                        if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) ) {
                            if (i*i+j*j < 3) {
                                if (label[x+i][y+j] > 0) {
                                    connect[Nconnect] = lb[ label[x+i][y+j] ];
                                    Nconnect++;
                                }
                            }
                        }
                    }
					// if connected values, find the smallest lb label and attribute it
					// to all others (-> join labels)
					if (Nconnect>0) {
						lbMin = lb[connect[0]];
						for (k=1;k<Nconnect;k++) lbMin = Math.min(lbMin,lb[connect[k]]);
						for (k=0;k<Nconnect;k++) lb[connect[k]] = lbMin;
						label[x][y] = lbMin;
					} else {
						// new, unconnected region
						label[x][y] = Nlabel;
						lb[Nlabel] = Nlabel;
						Nlabel++;
						// check if the number of labels is above the threshold
						if (Nlabel>=lb.length-1) {
							int[] tmp = new int[2*lb.length];
							for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
							lb = tmp;
						}
					}
				}
			}
		}
		// only one level of labels
		for (k=1;k<Nlabel;k++) {
			c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// copy on label image
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				label[x][y] = lb2[ lb[ label[x][y] ] ];
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}

	/** 
	 *	Connected components of an object.
     *  3D images: 6-neighborhood
	 *  caution: there seems to be a problem with some configurations
	 *	(two connected regions are labelled differently
	 */
	/* 
	public static final int[][][] connected6Object3Dwrong(boolean img[][][], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int x,y,z,i,j,k,l,c,n;
		int Nlb;
		int[]   connect = new int[7];
		int Nconnect;
		int AddLabel;
		
		for (x=0;x<nx;x++)
			for (y=0;y<ny;y++)
				for (z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					if (img[x][y][z]) {
						// object point: neighbors ?
						Nconnect = 0;
                        for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 2) {
                                    if (label[x+i][y+j][z+k] > 0) {
                                        connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
                                        Nconnect++;
                                    }
                                }
                            }
                        }
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[x][y][z] = lbMin;
						} else {
							// new, unconnected region
							label[x][y][z] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (k=1;k<Nlabel;k++) {
			c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// copy on label image
        for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
                for (z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	*/
	
	/** 
	 *	Connected components of an object.
     *  3D images: 6-neighborhood
	 *  slower but exact method (hopefully)
	 */
	public static final int[][][] connected6Object3D(boolean img[][][], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[6];
		int Nconnect;
		int AddLabel;
		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if (img[x][y][z]) {
						// object point: neighbors ?
						Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 2) {
                                    if (label[x+i][y+j][z+k] > 0) {
                                        connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
                                        Nconnect++;
                                    }
                                }
                            }
                        }
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[x][y][z] = lbMin;
						} else {
							// new, unconnected region
							label[x][y][z] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (int k=1;k<Nlabel;k++) {
			int c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					if (lb2[ lb[ label[x][y][z] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 2) {
									if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>0) {
										if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>lb2[ lb[ label[x][y][z] ] ]) {
											//System.out.print("labelling problem!");
											int badLb = lb2[ lb[ label[x+i][y+j][z+k] ] ];
											lb2[ lb[ label[x+i][y+j][z+k] ] ] = lb2[ lb[ label[x][y][z] ] ];
											// if there are higher labels, decrease them
											for (int n=1;n<Nlabel;n++) {
												if (lb2[n] > badLb) lb2[n]--;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	/** 
	 *	Connected components of an object.
     *  3D images: 6-neighborhood
	 *  slower but exact method (hopefully)
	 */
	public static final int[] connected6Object3D(boolean img[], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[]   label = new int[nx*ny*nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[6];
		int Nconnect;
		int AddLabel;
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++)
			label[xyz] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (img[xyz]) {
						// object point: neighbors ?
						Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 2) {
                                    if (label[xyz+i+nx*j+nx*ny*k] > 0) {
                                        connect[Nconnect] = lb[ label[xyz+i+nx*j+nx*ny*k] ];
                                        Nconnect++;
                                    }
                                }
                            }
                        }
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[xyz] = lbMin;
						} else {
							// new, unconnected region
							label[xyz] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (int k=1;k<Nlabel;k++) {
			int c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					if (lb2[ lb[ label[xyz] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 2) {
									if (lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ]>0) {
										if (lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ]>lb2[ lb[ label[xyz] ] ]) {
											//System.out.print("labelling problem!");
											int badLb = lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ];
											lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ] = lb2[ lb[ label[xyz] ] ];
											// if there are higher labels, decrease them
											for (int n=1;n<Nlabel;n++) {
												if (lb2[n] > badLb) lb2[n]--;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			label[xyz] = lb2[ lb[ label[xyz] ] ];
        }    
		
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	
    /**
	 *	Connected components of an object.
     *  3D images: 18-neighborhood
     */
    public static final int[][][] connected18Object3D(boolean img[][][], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[18];
		int Nconnect;
		int AddLabel;
		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if (img[x][y][z]) {
                        // object point: neighbors ?
                        Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 3) {
                                    if (label[x+i][y+j][z+k] > 0) {
                                        connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
                                        Nconnect++;
                                    }
                                }
                            }
                        }
                        // if connected values, find the smallest lb label and attribute it
                        // to all others (-> join labels)
                        if (Nconnect>0) {
                            //printf("c:%d",Nconnect);
                            lbMin = lb[connect[0]];
                            for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                            for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                            label[x][y][z] = lbMin;
                        } else {
                            // new, unconnected region
                            label[x][y][z] = Nlabel;
                            lb[Nlabel] = Nlabel;
                            //printf("l:%d", Nlabel);
                            Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
                        }
                    }
                }
            }
        }
        // only one level of labels
        for (int k=1;k<Nlabel;k++) {
            int c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					if (lb2[ lb[ label[x][y][z] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 3) {
									if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>0) {
										if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>lb2[ lb[ label[x][y][z] ] ]) {
											//System.out.print("labelling problem!");
											int badLb = lb2[ lb[ label[x+i][y+j][z+k] ] ];
											lb2[ lb[ label[x+i][y+j][z+k] ] ] = lb2[ lb[ label[x][y][z] ] ];
											// if there are higher labels, decrease them
											for (int n=1;n<Nlabel;n++) {
												if (lb2[n] > badLb) lb2[n]--;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
    }
    /**
	 *	Connected components of an object.
     *  3D images: 18-neighborhood
     */
    public static final int[] connected18Object3D(boolean img[], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[]   label = new int[nx*ny*nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[18];
		int Nconnect;
		int AddLabel;
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++)
			label[xyz] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					int xyz = x + nx*y + nx*ny*z;
					if (img[xyz]) {
                        // object point: neighbors ?
                        Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 3) {
                                    if (label[xyz+i+nx*j+nx*ny*k] > 0) {
                                        connect[Nconnect] = lb[ label[xyz+i+nx*j+nx*ny*k] ];
                                        Nconnect++;
                                    }
                                }
                            }
                        }
                        // if connected values, find the smallest lb label and attribute it
                        // to all others (-> join labels)
                        if (Nconnect>0) {
                            //printf("c:%d",Nconnect);
                            lbMin = lb[connect[0]];
                            for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                            for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                            label[xyz] = lbMin;
                        } else {
                            // new, unconnected region
                            label[xyz] = Nlabel;
                            lb[Nlabel] = Nlabel;
                            //printf("l:%d", Nlabel);
                            Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
                        }
                    }
                }
            }
        }
        // only one level of labels
        for (int k=1;k<Nlabel;k++) {
            int c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					int xyz = x + nx*y + nx*ny*z;
					if (lb2[ lb[ label[xyz] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 3) {
									if (lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ]>0) {
										if (lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ]>lb2[ lb[ label[xyz] ] ]) {
											//System.out.print("labelling problem!");
											int badLb = lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ];
											lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ] = lb2[ lb[ label[xyz] ] ];
											// if there are higher labels, decrease them
											for (int n=1;n<Nlabel;n++) {
												if (lb2[n] > badLb) lb2[n]--;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			label[xyz] = lb2[ lb[ label[xyz] ] ];
        }
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
    }

    /**
     *	Connected components of an object.
     *  3D images: 26-neighborhood
     */
    public static final int[][][] connected26Object3D(boolean img[][][], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[26];
		int Nconnect;
		int AddLabel;
		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if (img[x][y][z]) {
                        // object point: neighbors ?
                        Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (label[x+i][y+j][z+k] > 0) {
                                    connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
                                    Nconnect++;
                                }
                            }
                        }
                        // if connected values, find the smallest lb label and attribute it
                        // to all others (-> join labels)
                        if (Nconnect>0) {
                            //printf("c:%d",Nconnect);
                            lbMin = lb[connect[0]];
                            for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                            for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                            label[x][y][z] = lbMin;
                        } else {
                            // new, unconnected region
                            label[x][y][z] = Nlabel;
                            lb[Nlabel] = Nlabel;
                            //printf("l:%d", Nlabel);
                            Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
                        }
                    }
                }
            }
        }
        // only one level of labels
        for (int k=1;k<Nlabel;k++) {
            int c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					if (lb2[ lb[ label[x][y][z] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>0) {
									if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>lb2[ lb[ label[x][y][z] ] ]) {
										//System.out.print("labelling problem!");
										int badLb = lb2[ lb[ label[x+i][y+j][z+k] ] ];
										lb2[ lb[ label[x+i][y+j][z+k] ] ] = lb2[ lb[ label[x][y][z] ] ];
										// if there are higher labels, decrease them
										for (int n=1;n<Nlabel;n++) {
											if (lb2[n] > badLb) lb2[n]--;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
    }   
    /**
	 *	Connected components of an object.
     *  3D images: 26-neighborhood
     */
    public static final int[] connected26Object3D(boolean img[], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[]   label = new int[nx*ny*nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[26];
		int Nconnect;
		int AddLabel;
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++)
			label[xyz] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					int xyz = x + nx*y + nx*ny*z;
					if (img[xyz]) {
                        // object point: neighbors ?
                        Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (label[xyz+i+nx*j+nx*ny*k] > 0) {
									connect[Nconnect] = lb[ label[xyz+i+nx*j+nx*ny*k] ];
									Nconnect++;
								}
							}
                        }
                        // if connected values, find the smallest lb label and attribute it
                        // to all others (-> join labels)
                        if (Nconnect>0) {
                            //printf("c:%d",Nconnect);
                            lbMin = lb[connect[0]];
                            for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                            for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                            label[xyz] = lbMin;
                        } else {
                            // new, unconnected region
                            label[xyz] = Nlabel;
                            lb[Nlabel] = Nlabel;
                            //printf("l:%d", Nlabel);
                            Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
                        }
                    }
                }
            }
        }
        // only one level of labels
        for (int k=1;k<Nlabel;k++) {
            int c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					int xyz = x + nx*y + nx*ny*z;
					if (lb2[ lb[ label[xyz] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ]>0) {
									if (lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ]>lb2[ lb[ label[xyz] ] ]) {
										//System.out.print("labelling problem!");
										int badLb = lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ];
										lb2[ lb[ label[xyz+i+nx*j+nx*ny*k] ] ] = lb2[ lb[ label[xyz] ] ];
										// if there are higher labels, decrease them
										for (int n=1;n<Nlabel;n++) {
											if (lb2[n] > badLb) lb2[n]--;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			label[xyz] = lb2[ lb[ label[xyz] ] ];
        }
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
    }
	/** 
	 *	Connected components in images separated by edges
     *  3D images: 6-neighborhood
	 *  slower but exact method (hopefully)
	 */
	public static final int[][][] connected6Object3D(float img[][][], float threshold, boolean[][][] mask, int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[6];
		int Nconnect;
		int AddLabel;
		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if (mask[x][y][z]) {
						// object point: neighbors ?
						Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 2) {
                                	if (Numerics.abs(img[x+i][y+j][z+k]-img[x][y][z])<threshold) {
										if (label[x+i][y+j][z+k] > 0) {
											connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
											Nconnect++;
										}
									}
                                }
                            }
                        }
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[x][y][z] = lbMin;
						} else {
							// new, unconnected region
							label[x][y][z] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (int k=1;k<Nlabel;k++) {
			int c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					if (lb2[ lb[ label[x][y][z] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 2) {
									if (Numerics.abs(img[x+i][y+j][z+k]-img[x][y][z])<threshold) {
										if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>0) {
											if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>lb2[ lb[ label[x][y][z] ] ]) {
												//System.out.print("labelling problem!");
												int badLb = lb2[ lb[ label[x+i][y+j][z+k] ] ];
												lb2[ lb[ label[x+i][y+j][z+k] ] ] = lb2[ lb[ label[x][y][z] ] ];
												// if there are higher labels, decrease them
												for (int n=1;n<Nlabel;n++) {
													if (lb2[n] > badLb) lb2[n]--;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	
	/** 
	 *	Connected components in images separated by edges
     *  3D images: 6-neighborhood
	 *  slower but exact method (hopefully)
	 */
	public static final int[][][] connected6Object3D(int img[][][], int nx, int ny, int nz) {
		int Nlabel = 0;
		int[][][]   label = new int[nx][ny][nz];
		int[]       lb = new int[MaxObject];
		int lbMin;
		int Nlb;
		int[]   connect = new int[6];
		int Nconnect;
		int AddLabel;
		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if (img[x][y][z]>0) {
						// object point: neighbors ?
						Nconnect = 0;
                        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
                                if (i*i+j*j+k*k < 2) {
                                	if (img[x+i][y+j][z+k]==img[x][y][z]) {
										if (label[x+i][y+j][z+k] > 0) {
											connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
											Nconnect++;
										}
									}
                                }
                            }
                        }
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (int l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (int l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[x][y][z] = lbMin;
						} else {
							// new, unconnected region
							label[x][y][z] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
							// check if the number of labels is above the threshold
							if (Nlabel>=lb.length-1) {
								int[] tmp = new int[2*lb.length];
								for (int n=0;n<Nlabel;n++) tmp[n] = lb[n];
								lb = tmp;
							}
						}
					}
				}
			}
		}
		// only one level of labels
		for (int k=1;k<Nlabel;k++) {
			int c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		
		// count the valid labels and rearrange labels to have regular increment
		Nlb = 0;
		int[] lb2 = new int[Nlabel];
		lb2[0] = 0;
		for (int k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				Nlb++;
				lb2[k] = Nlb;
			}
		}
		// check for problems ?
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
					if (lb2[ lb[ label[x][y][z] ] ]>0) {
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) {
							if ( (x+i>=0) && (x+i<nx) && (y+j>=0) && (y+j<ny) && (z+k>=0) && (z+k<nz) ) {
								if (i*i+j*j+k*k < 2) {
									if (img[x+i][y+j][z+k]==img[x][y][z]) {
										if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>0) {
											if (lb2[ lb[ label[x+i][y+j][z+k] ] ]>lb2[ lb[ label[x][y][z] ] ]) {
												//System.out.print("labelling problem!");
												int badLb = lb2[ lb[ label[x+i][y+j][z+k] ] ];
												lb2[ lb[ label[x+i][y+j][z+k] ] ] = lb2[ lb[ label[x][y][z] ] ];
												// if there are higher labels, decrease them
												for (int n=1;n<Nlabel;n++) {
													if (lb2[n] > badLb) lb2[n]--;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		// copy on label image
        for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
                for (int z=0;z<nz;z++) {
                    label[x][y][z] = lb2[ lb[ label[x][y][z] ] ];
                }
			}
		}
        // clean up
        lb = null;
		lb2 = null;
        connect = null;
	   
		return label;
	}
	
}
