package de.mpg.cbs.structures;

import java.io.*;
import java.util.*;
import java.util.zip.*;
import java.lang.*;
import java.net.URL;

/**
 *
 *  This class handles a LUT for critical points in 3D topology
 *	with the connectivity from a given filename
 *	
 *	@version    June 2005
 *	@author     Pierre-Louis Bazin
 *
 */

public class CriticalPointLUT {
	
	private static final int B0 = 1;
	private static final int B1 = 2;
	private static final int B2 = 4;
	private static final int B3 = 8;
	private static final int B4 = 16;
	private static final int B5 = 32;
	private static final int B6 = 64;
	private static final int B7 = 128;
	private static final int B8 = 256;
	private static final int B9 = 512;
	private static final int B10 = 1024;
	private static final int B11 = 2048;
	private static final int B12 = 4096;
	private static final int B13 = 8192;
	private static final int B14 = 16384;
	private static final int B15 = 32768;
	private static final int B16 = 65536;
	private static final int B17 = 131072;
	private static final int B18 = 262144;
	private static final int B19 = 524288;
	private static final int B20 = 1048576;
	private static final int B21 = 2097152;
	private static final int B22 = 4194304;
	private static final int B23 = 8388608;
	private static final int B24 = 16777216;
	private static final int B25 = 33554432;
	private static final int B26 = 67108864;
	
	private static final int MAX = 200;

	private BitSet 		isRegular;
	
	private	String		filename;			// file to open
	private	String		filepath = null;	// root directory for the file
	
	private static final boolean	debug = true;
	
	public CriticalPointLUT(String filename_, int size) {
		isRegular = new BitSet(size);
		filename = filename_;
	}
	
	public CriticalPointLUT(String filepath_, String filename_, int size) {
		isRegular = new BitSet(size);
		filepath = filepath_;
		filename = filename_;
	}
	
	public void finalize() {
		isRegular = null;
	}
	
	/**
	 *  add a new value
	 */
	public final void set(int ind, boolean val) { isRegular.set(ind, val); }
	
	/**
	 *  get a value
	 */
	 public final boolean get(int ind) { return isRegular.get(ind); }
	
	 /** translate a pattern into the key number */
	 public final int keyFromPattern(byte[][][] img, int x, int y, int z) {
		 int key=0;
		 
		 if (img[x-1][y-1][z-1]>img[x][y][z]) key += B0;
		 if (img[x-1][y-1][z  ]>img[x][y][z]) key += B1;
		 if (img[x-1][y-1][z+1]>img[x][y][z]) key += B2;
		 
		 if (img[x-1][y  ][z-1]>img[x][y][z]) key += B3;
		 if (img[x-1][y  ][z  ]>img[x][y][z]) key += B4;
		 if (img[x-1][y  ][z+1]>img[x][y][z]) key += B5;
		 
		 if (img[x-1][y+1][z-1]>img[x][y][z]) key += B6;
		 if (img[x-1][y+1][z  ]>img[x][y][z]) key += B7;
		 if (img[x-1][y+1][z+1]>img[x][y][z]) key += B8;
		 
		 if (img[x  ][y-1][z-1]>img[x][y][z]) key += B9;
		 if (img[x  ][y-1][z  ]>img[x][y][z]) key += B10;
		 if (img[x  ][y-1][z+1]>img[x][y][z]) key += B11;
		 
		 if (img[x  ][y  ][z-1]>img[x][y][z]) key += B12;
		 if (img[x  ][y  ][z+1]>img[x][y][z]) key += B13;
		 
		 if (img[x  ][y+1][z-1]>img[x][y][z]) key += B14;
		 if (img[x  ][y+1][z  ]>img[x][y][z]) key += B15;
		 if (img[x  ][y+1][z+1]>img[x][y][z]) key += B16;
		 
		 if (img[x+1][y-1][z-1]>img[x][y][z]) key += B17;
		 if (img[x+1][y-1][z  ]>img[x][y][z]) key += B18;
		 if (img[x+1][y-1][z+1]>img[x][y][z]) key += B19;
		 
		 if (img[x+1][y  ][z-1]>img[x][y][z]) key += B20;
		 if (img[x+1][y  ][z  ]>img[x][y][z]) key += B21;
		 if (img[x+1][y  ][z+1]>img[x][y][z]) key += B22;
		 
		 if (img[x+1][y+1][z-1]>img[x][y][z]) key += B23;
		 if (img[x+1][y+1][z  ]>img[x][y][z]) key += B24;
		 if (img[x+1][y+1][z+1]>img[x][y][z]) key += B25;
		 
		 return key;
	 }
	 
	 /** translate a pattern into the key number */
	 public final int keyFromPattern(byte[][][][] img, int x, int y, int z, int k) {
		 int key=0;
		 
		 if (img[x-1][y-1][z-1][k]>img[x][y][z][k]) key += B0;
		 if (img[x-1][y-1][z  ][k]>img[x][y][z][k]) key += B1;
		 if (img[x-1][y-1][z+1][k]>img[x][y][z][k]) key += B2;
		 
		 if (img[x-1][y  ][z-1][k]>img[x][y][z][k]) key += B3;
		 if (img[x-1][y  ][z  ][k]>img[x][y][z][k]) key += B4;
		 if (img[x-1][y  ][z+1][k]>img[x][y][z][k]) key += B5;
		 
		 if (img[x-1][y+1][z-1][k]>img[x][y][z][k]) key += B6;
		 if (img[x-1][y+1][z  ][k]>img[x][y][z][k]) key += B7;
		 if (img[x-1][y+1][z+1][k]>img[x][y][z][k]) key += B8;
		 
		 if (img[x  ][y-1][z-1][k]>img[x][y][z][k]) key += B9;
		 if (img[x  ][y-1][z  ][k]>img[x][y][z][k]) key += B10;
		 if (img[x  ][y-1][z+1][k]>img[x][y][z][k]) key += B11;
		 
		 if (img[x  ][y  ][z-1][k]>img[x][y][z][k]) key += B12;
		 if (img[x  ][y  ][z+1][k]>img[x][y][z][k]) key += B13;
		 
		 if (img[x  ][y+1][z-1][k]>img[x][y][z][k]) key += B14;
		 if (img[x  ][y+1][z  ][k]>img[x][y][z][k]) key += B15;
		 if (img[x  ][y+1][z+1][k]>img[x][y][z][k]) key += B16;
		 
		 if (img[x+1][y-1][z-1][k]>img[x][y][z][k]) key += B17;
		 if (img[x+1][y-1][z  ][k]>img[x][y][z][k]) key += B18;
		 if (img[x+1][y-1][z+1][k]>img[x][y][z][k]) key += B19;
		 
		 if (img[x+1][y  ][z-1][k]>img[x][y][z][k]) key += B20;
		 if (img[x+1][y  ][z  ][k]>img[x][y][z][k]) key += B21;
		 if (img[x+1][y  ][z+1][k]>img[x][y][z][k]) key += B22;
		 
		 if (img[x+1][y+1][z-1][k]>img[x][y][z][k]) key += B23;
		 if (img[x+1][y+1][z  ][k]>img[x][y][z][k]) key += B24;
		 if (img[x+1][y+1][z+1][k]>img[x][y][z][k]) key += B25;
		 
		 return key;
	 }
	 
	 /** translate a pattern into the key number */
	 public final int keyFromPattern(boolean[][][] img, int x, int y, int z) {
		 int key=0;
		 
		 if (img[x-1][y-1][z-1]) key += B0;
		 if (img[x-1][y-1][z  ]) key += B1;
		 if (img[x-1][y-1][z+1]) key += B2;
		 
		 if (img[x-1][y  ][z-1]) key += B3;
		 if (img[x-1][y  ][z  ]) key += B4;
		 if (img[x-1][y  ][z+1]) key += B5;
		 
		 if (img[x-1][y+1][z-1]) key += B6;
		 if (img[x-1][y+1][z  ]) key += B7;
		 if (img[x-1][y+1][z+1]) key += B8;
		 
		 if (img[x  ][y-1][z-1]) key += B9;
		 if (img[x  ][y-1][z  ]) key += B10;
		 if (img[x  ][y-1][z+1]) key += B11;
		 
		 if (img[x  ][y  ][z-1]) key += B12;
		 if (img[x  ][y  ][z+1]) key += B13;
		 
		 if (img[x  ][y+1][z-1]) key += B14;
		 if (img[x  ][y+1][z  ]) key += B15;
		 if (img[x  ][y+1][z+1]) key += B16;
		 
		 if (img[x+1][y-1][z-1]) key += B17;
		 if (img[x+1][y-1][z  ]) key += B18;
		 if (img[x+1][y-1][z+1]) key += B19;
		 
		 if (img[x+1][y  ][z-1]) key += B20;
		 if (img[x+1][y  ][z  ]) key += B21;
		 if (img[x+1][y  ][z+1]) key += B22;
		 
		 if (img[x+1][y+1][z-1]) key += B23;
		 if (img[x+1][y+1][z  ]) key += B24;
		 if (img[x+1][y+1][z+1]) key += B25;
		 
		 return key;
	 }
	 
	 /** translate a pattern into the key number */
	 public final int keyFromPattern(byte[] img, int idx, int dx, int dy, int dz) {
		 int key=0;
		 
		 if (img[idx-dx-dy-dz]>img[idx]) key += B0;
		 if (img[idx-dx-dy   ]>img[idx]) key += B1;
		 if (img[idx-dx-dy+dz]>img[idx]) key += B2;
		 
		 if (img[idx-dx   -dz]>img[idx]) key += B3;
		 if (img[idx-dx      ]>img[idx]) key += B4;
		 if (img[idx-dx   +dz]>img[idx]) key += B5;
		 
		 if (img[idx-dx+dy-dz]>img[idx]) key += B6;
		 if (img[idx-dx+dy   ]>img[idx]) key += B7;
		 if (img[idx-dx+dy+dz]>img[idx]) key += B8;
		 
		 if (img[idx   -dy-dz]>img[idx]) key += B9;
		 if (img[idx   -dy   ]>img[idx]) key += B10;
		 if (img[idx   -dy+dz]>img[idx]) key += B11;
		 
		 if (img[idx      -dz]>img[idx]) key += B12;
		 if (img[idx      +dz]>img[idx]) key += B13;
		 
		 if (img[idx   +dy-dz]>img[idx]) key += B14;
		 if (img[idx   +dy   ]>img[idx]) key += B15;
		 if (img[idx   +dy+dz]>img[idx]) key += B16;
		 
		 if (img[idx+dx-dy-dz]>img[idx]) key += B17;
		 if (img[idx+dx-dy   ]>img[idx]) key += B18;
		 if (img[idx+dx-dy+dz]>img[idx]) key += B19;
		 
		 if (img[idx+dx   -dz]>img[idx]) key += B20;
		 if (img[idx+dx      ]>img[idx]) key += B21;
		 if (img[idx+dx   +dz]>img[idx]) key += B22;
		 
		 if (img[idx+dx+dy-dz]>img[idx]) key += B23;
		 if (img[idx+dx+dy   ]>img[idx]) key += B24;
		 if (img[idx+dx+dy+dz]>img[idx]) key += B25;
		 
		 return key;
	 }
	 
	 /** translate a pattern into the key number */
	 public final int keyFromPattern(boolean[] img, int idx, int dx, int dy, int dz) {
		 int key=0;
		 
		 if (img[idx-dx-dy-dz]) key += B0;
		 if (img[idx-dx-dy   ]) key += B1;
		 if (img[idx-dx-dy+dz]) key += B2;
		 
		 if (img[idx-dx   -dz]) key += B3;
		 if (img[idx-dx      ]) key += B4;
		 if (img[idx-dx   +dz]) key += B5;
		 
		 if (img[idx-dx+dy-dz]) key += B6;
		 if (img[idx-dx+dy   ]) key += B7;
		 if (img[idx-dx+dy+dz]) key += B8;
		 
		 if (img[idx   -dy-dz]) key += B9;
		 if (img[idx   -dy   ]) key += B10;
		 if (img[idx   -dy+dz]) key += B11;
		 
		 if (img[idx      -dz]) key += B12;
		 if (img[idx      +dz]) key += B13;
		 
		 if (img[idx   +dy-dz]) key += B14;
		 if (img[idx   +dy   ]) key += B15;
		 if (img[idx   +dy+dz]) key += B16;
		 
		 if (img[idx+dx-dy-dz]) key += B17;
		 if (img[idx+dx-dy   ]) key += B18;
		 if (img[idx+dx-dy+dz]) key += B19;
		 
		 if (img[idx+dx   -dz]) key += B20;
		 if (img[idx+dx      ]) key += B21;
		 if (img[idx+dx   +dz]) key += B22;
		 
		 if (img[idx+dx+dy-dz]) key += B23;
		 if (img[idx+dx+dy   ]) key += B24;
		 if (img[idx+dx+dy+dz]) key += B25;
		 
		 return key;
	 }
	 
	 /** translate a key into the corresponding pattern */
	 public final byte[][][] patternFromKey(int key) {
		 byte[][][] img = new byte[3][3][3];
		 
		 img[2][2][2] = (byte)(key/B25); key -= img[2][2][2]*B25;
		 img[2][2][1] = (byte)(key/B24); key -= img[2][2][1]*B24;
		 img[2][2][0] = (byte)(key/B23); key -= img[2][2][0]*B23;
		 
		 img[2][1][2] = (byte)(key/B22); key -= img[2][1][2]*B22;
		 img[2][1][1] = (byte)(key/B21); key -= img[2][1][1]*B21;
		 img[2][1][0] = (byte)(key/B20); key -= img[2][1][0]*B20;
		 
		 img[2][0][2] = (byte)(key/B19); key -= img[2][0][2]*B19;
		 img[2][0][1] = (byte)(key/B18); key -= img[2][0][1]*B18;
		 img[2][0][0] = (byte)(key/B17); key -= img[2][0][0]*B17;
		 
		 img[1][2][2] = (byte)(key/B16); key -= img[1][2][2]*B16;
		 img[1][2][1] = (byte)(key/B15); key -= img[1][2][1]*B15;
		 img[1][2][0] = (byte)(key/B14); key -= img[1][2][0]*B14;
		 
		 img[1][1][2] = (byte)(key/B13); key -= img[1][1][2]*B13;
		 img[1][1][0] = (byte)(key/B12); key -= img[1][1][0]*B12;
		 
		 img[1][0][2] = (byte)(key/B11); key -= img[1][0][2]*B11;
		 img[1][0][1] = (byte)(key/B10); key -= img[1][0][1]*B10;
		 img[1][0][0] = (byte)(key/B9 ); key -= img[1][0][0]*B9 ;
		 
		 img[0][2][2] = (byte)(key/B8 ); key -= img[0][2][2]*B8 ;
		 img[0][2][1] = (byte)(key/B7 ); key -= img[0][2][1]*B7 ;
		 img[0][2][0] = (byte)(key/B6 ); key -= img[0][2][0]*B6 ;
		 
		 img[0][1][2] = (byte)(key/B5 ); key -= img[0][1][2]*B5 ;
		 img[0][1][1] = (byte)(key/B4 ); key -= img[0][1][1]*B4 ;
		 img[0][1][0] = (byte)(key/B3 ); key -= img[0][1][0]*B3 ;
		 
		 img[0][0][2] = (byte)(key/B2 ); key -= img[0][0][2]*B2 ;
		 img[0][0][1] = (byte)(key/B1 ); key -= img[0][0][1]*B1 ;
		 img[0][0][0] = (byte)(key/B0 ); key -= img[0][0][0]*B0 ;
		 
		 return img;
	 }

	 public final void loadPattern() {
		 //URL res = CriticalPointLUT.class.getResource( filename );
		 //String current = System.getProperty("user.dir");
		 //System.out.println("path: "+current);
		 //System.setProperty("user.dir", res.getPath());
		 //System.out.println("pattern path: "+System.getProperty("user.dir"));
		 try {
			 //File f = new File( res.getFile() );
			 //File f = new File( filename );
			 //FileInputStream fis = new FileInputStream( f );
			 //if (debug) System.out.println("Opening LUT: "+res.toString());
			 //InputStream fis = res.openStream();
			 InputStream fis = CriticalPointLUT.class.getResourceAsStream( filename );
			 byte[] list = new byte[B26];
			
			 fis.read(list);
			 fis.close();
			 int N = 0;
			 for (int n=0;n<B26;n++) {
				 isRegular.set( n, (list[n]==1) );
			 	 if (list[n]==1) N++;
			 }
			 if (debug) System.out.println("Simple points: "+N);
			 
			 list = null;
		 } catch (FileNotFoundException e) {
			 System.out.println("File not found:");
			 System.out.println(e.getMessage());
		 } catch (IOException e ) {
			 System.out.println("i/o exception:");
			 System.out.println(e.getMessage());
		 }
	 }
	 
	 public final boolean loadCompressedPattern() {
		 /*
		 URL res = null;
		 try {
			res = CriticalPointLUT.class.getResource( filename );
		 } catch (Exception e) {
			System.err.println("error:");
            System.err.println(e.getMessage()+ "\n");
		 }
		 */
		 //String current = System.getProperty("user.dir");
		 //System.out.println("path: "+current);
		 //System.setProperty("user.dir", res.getPath());
		 //System.out.println("pattern path: "+System.getProperty("user.dir"));
		 try {
			 //File f = new File( res.getFile() );
			 //File f = new File( filename );
			 //FileInputStream fis = new FileInputStream( f );
			 //if (debug) System.out.println("Opening LUT: "+res.toString());
			 //FileInputStream fis = new FileInputStream( res.getPath() );
			 //FileInputStream fis = res.openStream();
			 InputStream fis;
			 if (filepath==null) fis = CriticalPointLUT.class.getResourceAsStream( filename );
			 else fis = new FileInputStream( filepath+filename );
			 System.out.println("Opening LUT: "+fis.toString());
			 System.out.flush();
			 GZIPInputStream gzfis = new GZIPInputStream(fis, B26);
			 //GZIPInputStream gzfis = new GZIPInputStream(res.openStream(), B26);
			 byte[] list = new byte[B26];
			
			 gzfis.read(list, 0, B26);
			 gzfis.close();
			 int N = 0;
			 for (int n=0;n<B26;n++) {
				 isRegular.set( n, (list[n]==1) );
				 if (list[n]==1) N++;
			 }
			 if (debug) System.out.println("Simple points: "+N);
			 
			 list = null;
			 return true;
		 } catch (FileNotFoundException e) {
			 System.out.println("File not found:");
			 System.out.println(e.getMessage());
			 isRegular = null;
			 return false;
		 } catch (IOException e ) {
			 System.out.println("i/o exception:");
			 System.out.println(e.getMessage());
			 isRegular = null;
			 return false;
		 }
	 }
	 
	 public final String getFilename() {
		 URL res = CriticalPointLUT.class.getResource( filename );
		 return res.toString();
	 }
}
