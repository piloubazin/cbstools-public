package de.mpg.cbs.utilities;


/**
 *
 *  This class computes various basic numerical functions.
 *	<p> 
 *	Includes min, max, gaussian random numbers, special functions, 
 *	mathematical constants. This duplicates some of Java's Math functions
 *	where these have limited type support (e.g. double, but not float) or
 *	perform multiple casting operations (slow).
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 */

public class Numerics {
	
	/** mathematical constants */
	public static final float ZERO = 1e-30f;
	/** mathematical constants */
	public static final float INF = 1e+30f;
	/** mathematical constants */
	public static final float PI = 3.1416f;
	
	/**
	 *	rounding function (faster than regular Java)
	 */
    public static final int round(float num) {
		if (num==(int)num) return (int)num;
		else if (num>0) return (int)(num+0.5f);
		else return (int)(num-0.5f);
    }
    public static final int floor(float num) {
		if (num==(int)num) return (int)num;
		else if (num>0) return (int)(num);
		else return (int)(num)-1;
	}
    public static final int ceil(float num) {		
		if (num==(int)num) return (int)num;
		else if (num>0) return (int)(num)+1;
		else return (int)(num);
    }
	
    public static final int round(double num) {		
		if (num==(int)num) return (int)num;
		else if (num>0) return (int)(num+0.5);
		else return (int)(num-0.5);
    }
    public static final int floor(double num) {		
		if (num==(int)num) return (int)num;
		else if (num>0) return (int)(num);
		else return (int)(num)-1;
    }
    public static final int ceil(double num) {		
		if (num==(int)num) return (int)num;
		else if (num>0) return (int)(num)+1;
		else return (int)(num);
    }
	/* min, max */
	public static final float min( float a, float b) {
		if (a < b) return a;
		else return b;
	}
	public static final double min( double a, double b) {
		if (a < b) return a;
		else return b;
	}
	public static final byte min( byte a, byte b) {
		if (a < b) return a;
		else return b;
	}
	public static final short min( short a, short b) {
		if (a < b) return a;
		else return b;
	}
	public static final int min( int a, int b) {
		if (a < b) return a;
		else return b;
	}
	public static final long min( long a, long b) {
		if (a < b) return a;
		else return b;
	}
	public static final float min( float a, float b, float c) {
		return min(a,min(b,c));
	}
	public static final double min( double a, double b, double c) {
		return min(a,min(b,c));
	}
	public static final double min( double[] a) {
		double m = a[0];
		for (int n=1;n<a.length;n++) if (a[n]<m) m = a[n];
		return m;
	}
	public static final float min( float[] a) {
		float m = a[0];
		for (int n=1;n<a.length;n++) if (a[n]<m) m = a[n];
		return m;
	}
	public static final float min( float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8,
									float x9, float x10, float x11, float x12, float x13, float x14, float x15, float x16,
									float x17, float x18, float x19, float x20, float x21, float x22, float x23, float x24 ) {
		return min(x1,min(x2,min(x3,min(x4,min(x5,min(x6,min(x7,min(x8,
				min(x9,min(x10,min(x11,min(x12,min(x13,min(x14,min(x15,
				 min(x16,min(x17,min(x18,min(x19,min(x20,min(x21,min(x22,min(x23,x24)))))))))))))))))))))));
	}
	public static final float min( float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8,
									float x9, float x10, float x11, float x12, float x13, float x14, float x15, float x16,
									float x17, float x18 ) {
		return min(x1,min(x2,min(x3,min(x4,min(x5,min(x6,min(x7,min(x8,
				min(x9,min(x10,min(x11,min(x12,min(x13,min(x14,min(x15,
				 min(x16,min(x17,x18)))))))))))))))));
	}
	public static final float min( float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8 ) {
		return min(x1,min(x2,min(x3,min(x4,min(x5,min(x6,min(x7,x8)))))));
	}
	public static final float min( float x1, float x2, float x3, float x4, float x5, float x6 ) {
		return min(x1,min(x2,min(x3,min(x4,min(x5,x6)))));
	}
	public static final float min( float x1, float x2, float x3, float x4 ) {
		return min(x1,min(x2,min(x3,x4)));
	}

	public static final float max( float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8,
									float x9, float x10, float x11, float x12, float x13, float x14, float x15, float x16,
									float x17, float x18, float x19, float x20, float x21, float x22, float x23, float x24 ) {
		return max(x1,max(x2,max(x3,max(x4,max(x5,max(x6,max(x7,max(x8,
				max(x9,max(x10,max(x11,max(x12,max(x13,max(x14,max(x15,
				 max(x16,max(x17,max(x18,max(x19,max(x20,max(x21,max(x22,max(x23,x24)))))))))))))))))))))));
	}
	public static final float max( float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8,
									float x9, float x10, float x11, float x12, float x13, float x14, float x15, float x16,
									float x17, float x18 ) {
		return max(x1,max(x2,max(x3,max(x4,max(x5,max(x6,max(x7,max(x8,
				max(x9,max(x10,max(x11,max(x12,max(x13,max(x14,max(x15,
				 max(x16,max(x17,x18)))))))))))))))));
	}
	public static final float max( float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
		return max(x1,max(x2,max(x3,max(x4,max(x5,max(x6,max(x7,x8)))))));
	}

	public static final float max( float[][][] img, int x, int y, int z, int d) {
		float max = -INF;
		for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int k=-d;k<=d;k++) {
			if (img[x+i][y+j][z+k]>max) max = img[x+i][y+j][z+k];
		}
		return max;
	}
	public static final float min( float[][][] img, int x, int y, int z, int d) {
		float min = INF;
		for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int k=-d;k<=d;k++) {
			if (img[x+i][y+j][z+k]<min) min = img[x+i][y+j][z+k];
		}
		return min;
	}

	public static final float max( float a, float b) {
		if (a > b) return a;
		else return b;
	}
	public static final double max( double a, double b) {
		if (a > b) return a;
		else return b;
	}
	public static final int max( int a, int b) {
		if (a > b) return a;
		else return b;
	}
	public static final short max( short a, short b) {
		if (a > b) return a;
		else return b;
	}
	public static final int max( int a, int b, int c) {
		return max(a,max(b,c));
	}
	public static final double max( double a, double b, double c) {
		return max(a,max(b,c));
	}
	public static final double max( double a, double b, double c, double d) {
		return max(a,max(b,max(c,d)));
	}
	public static final double max( double a, double b, double c, double d, double e, double f) {
		return max(a,max(b,max(c,max(d,max(e,f)))));
	}
	public static final float max( float a, float b, float c) {
		return max(a,max(b,c));
	}
	public static final float max( float a, float b, float c, float d, float e) {
		return max(a,max(b,max(c,max(d,e))));
	}
	public static final double max( double[] val) {
		double max = val[0];
		for (int n=1;n<val.length;n++) if (val[n]>max) max = val[n];
		return max;
	}
	public static final float max( float[] val) {
		float max = val[0];
		for (int n=1;n<val.length;n++) if (val[n]>max) max = val[n];
		return max;
	}
	public static final int max( int[] val) {
		int max = val[0];
		for (int n=1;n<val.length;n++) if (val[n]>max) max = val[n];
		return max;
	}
	public static final byte max( byte[] val) {
		byte max = val[0];
		for (int n=1;n<val.length;n++) if (val[n]>max) max = val[n];
		return max;
	}
	public static final float max( float[] val, int first, int last) {
		float max = val[first];
		for (int n=first+1;n<last;n++) if (val[n]>max) max = val[n];
		return max;
	}
	
	public static final short bounded( short x, short a, short b) {
		return max(a,min(x,b));
	}
	
	public static final int bounded( int x, int a, int b) {
		return max(a,min(x,b));
	}
	
	public static final float bounded( float x, float a, float b) {
		return max(a,min(x,b));
	}
	
	public static final double bounded( double x, double a, double b) {
		return max(a,min(x,b));
	}

	public static final float minmag( float a, float b) {
		if (a*a < b*b) return a;
		else return b;
	}
	
	public static final float minmag( float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8 ) {
		float min = x1;
		if (x2*x2<min*min) min = x2;
		if (x3*x3<min*min) min = x3;
		if (x4*x4<min*min) min = x4;
		if (x5*x5<min*min) min = x5;
		if (x6*x6<min*min) min = x6;
		if (x7*x7<min*min) min = x7;
		if (x8*x8<min*min) min = x8;
		return min;
	}
	public static final float minmag( float x1, float x2, float x3, float x4, float x5, float x6 ) {
		float min = x1;
		if (x2*x2<min*min) min = x2;
		if (x3*x3<min*min) min = x3;
		if (x4*x4<min*min) min = x4;
		if (x5*x5<min*min) min = x5;
		if (x6*x6<min*min) min = x6;
		return min;
	}
	public static final float minmag( float x1, float x2, float x3, float x4 ) {
		float min = x1;
		if (x2*x2<min*min) min = x2;
		if (x3*x3<min*min) min = x3;
		if (x4*x4<min*min) min = x4;
		return min;
	}
	public static final float minmag( float x1, float x2, float x3 ) {
		float min = x1;
		if (x2*x2<min*min) min = x2;
		if (x3*x3<min*min) min = x3;
		return min;
	}

	public static final float minmag( float[] x, int last) {
		float min = x[0];
		for (int n=1;n<last;n++) if (x[n]*x[n]<min*min) min = x[n];
		return min;
	}
	public static final float maxmag( float a, float b) {
		if (a*a > b*b) return a;
		else return b;
	}
	
	public static final double maxmag( double a, double b) {
		if (a*a > b*b) return a;
		else return b;
	}
	
	/** returns x in [a,b] linearly mapped to [0,1] */
	public static final float mapped( float x, float a, float b) {
		return (max(a,min(x,b))-a)/(b-a);
	}
	
	/* absolute value, sign */
	public static final float abs( float a) {
		if (a > 0) return a;
		else return -a;
	}
	public static final double abs( double a) {
		if (a > 0) return a;
		else return -a;
	}
	public static final int abs( int a) {
		if (a > 0) return a;
		else return -a;
	}
	public static final int sign( float a) {
		if (a > 0) return 1;
		else if (a < 0) return -1;
		else return 0;
	}
	public static final int sign( double a) {
		if (a > 0) return 1;
		else if (a < 0) return -1;
		else return 0;
	}
	
	public static final float pow( float a, float b) {
			 if (b==0.0f) 		return 1.0f;
		else if (b==1.0f)	 	return a;
		else if (b==2.0f)		return a*a;
		else if (b==3.0f)	 	return a*a*a;
		else if (b==0.5f)		return (float)Math.sqrt(a);
		else if (b==1.0f/3.0f)	return (float)Math.cbrt(a);
		else if (b==2.0f/3.0f)	return (float)Math.cbrt(a*a);
		else					return (float)Math.pow(a,b);
	}
	
	public static final float inverse( float a) {
		if (a>ZERO || a<-ZERO) return 1.0f/a;
		else return INF;
	}
	
	/** find the index of the smallest values in val */
	public static final byte argmin(float[] val) {
		byte nmin=0;
		for (byte m=1;m<val.length;m++) if (val[m]<val[nmin]) {
			nmin = m;
		}
		return nmin;
	}	
			
	/** gives the position of the lowest value*/
	public static final int argmin(float a, float b) {
		if (a==b) return -1;
		else if (a<b) return 0;
		else return 1;
	}	
			
	/** gives the position of the lowest value*/
	public static final int argmin(float a, float b, float c) {
		if (a==b && b==c) return -1;
		if (a<b) {
			if (a<c) return 0;
			else return 2;
		} else { 
			if (b<c) return 1;
			else return 2;
		}
	}	
			
	/** find the index of the largest values in val */
	public static final int argmax(float[] val) {
		int nmax=0;
		for (int m=1;m<val.length;m++) if (val[m]>val[nmax]) {
			nmax = m;
		}
		return nmax;
	}	
			
	/** find the indices of the num largest values in val (>=0) */
	public static final byte[] argmax(float[] val, int num) {
		byte[] id = new byte[num];
		for (int n=0;n<num;n++) {
			byte nmax=0;
			for (byte m=1;m<val.length;m++) if (val[m]>val[nmax]) {
				nmax = m;
			}
			id[n] = nmax;
			val[nmax] *= -1;
		}
		// rewrite the values
		for (int n=0;n<num;n++) {
			val[id[n]] *= -1;
		}
		return id;
	}	
			
	/** find the indices of the num largest values in val (>=0) */
	public static final byte[] argmax(float[] val, int size, int num) {
		byte[] id = new byte[num];
		for (int n=0;n<num;n++) {
			byte nmax=0;
			for (byte m=1;m<size;m++) if (val[m]>val[nmax]) {
				nmax = m;
			}
			id[n] = nmax;
			val[nmax] *= -1;
		}
		// rewrite the values
		for (int n=0;n<num;n++) {
			val[id[n]] *= -1;
		}
		return id;
	}	
			
	/** find the indices of the num smallest values in val */
	public static final byte[] argmin(float[] val, int size, int num) {
		float max = val[0];
		for (byte m=1;m<size;m++) if (val[m]>max) max = val[m];
		byte[] id = new byte[num];
		for (int n=0;n<num;n++) {
			byte nmin=0;
			for (byte m=1;m<size;m++) if (val[m]<val[nmin]) {
				nmin = m;
			}
			id[n] = nmin;
			val[nmin] += max+1.0f;
		}
		// rewrite the values
		for (int n=0;n<num;n++) {
			val[id[n]] -= max+1.0f;
		}
		return id;
	}	
			
	/** find the indices of the num largest values in val (>=0) */
	public static final void argmax(byte[] id, float[] val, int num) {
		for (int n=0;n<num;n++) {
			byte nmax=0;
			for (byte m=1;m<val.length;m++) if (val[m]>val[nmax]) {
				nmax = m;
			}
			id[n] = nmax;
			val[nmax] *= -1;
		}
		// rewrite the values
		for (int n=0;n<num;n++) {
			val[id[n]] *= -1;
		}
		return;
	}	
			
	/** find the indices of the num largest values in val (any sign) */
	public static final void argmax(byte[] id, float[] bestval, float[] val, int num) {
		for (int n=0;n<num;n++) {
			byte nmax=0;
			for (byte m=1;m<val.length;m++) if (val[m]>val[nmax]) {
				nmax = m;
			}
			id[n] = nmax;
			bestval[n] = val[nmax];
			val[nmax] = -INF;
		}
		return;
	}	
			
	/** find the indices of the num largest values in val (any sign) */
	public static final void argmax(short[] id, float[] bestval, float[] val, int num) {
		for (int n=0;n<num;n++) {
			short nmax=0;
			for (short m=1;m<val.length;m++) if (val[m]>val[nmax]) {
				nmax = m;
			}
			id[n] = nmax;
			bestval[n] = val[nmax];
			val[nmax] = -INF;
		}
		return;
	}	
			
	/** find the indices of the num largest values in val (any sign) */
	public static final void argmax(byte[] id, double[] bestval, double[] val, int num) {
		for (int n=0;n<num;n++) {
			byte nmax=0;
			for (byte m=1;m<val.length;m++) if (val[m]>val[nmax]) {
				nmax = m;
			}
			id[n] = nmax;
			bestval[n] = val[nmax];
			val[nmax] = -INF;
		}
		return;
	}	
	
	/** find the indices of the num largest values in val (>=0) */
	public static final void argmax(byte[] id, float[] bestval, float[] val, int num, byte first) {
		for (int n=0;n<num;n++) {
			byte nmax=first;
			for (byte m=0;m<val.length;m++) if (val[m]>val[nmax]) {
				nmax = m;
			}
			id[n] = nmax;
			bestval[n] = val[nmax];
			val[nmax] = -INF;
		}
		return;
	}	
			
	/** sort the array by decreasing values (destructive!) */
	public static final void sortDown(float[] val) {
		float[] sorted = new float[val.length];
		for (int n=0;n<val.length;n++) {
			byte nmax=0;
			for (byte m=1;m<val.length;m++) if (val[m]>val[nmax]) {
				nmax = m;
			}
			sorted[n] = val[nmax];
			val[nmax] = -INF;
		}
		for (int n=0;n<val.length;n++) {
			val[n] = sorted[n];
		}
		return;
	}	
			
	/** sort the array by increasing values (destructive!) */
	public static final void sortUp(float[] val) {
		float[] sorted = new float[val.length];
		for (int n=0;n<val.length;n++) {
			byte nmax=0;
			for (byte m=1;m<val.length;m++) if (val[m]<val[nmax]) {
				nmax = m;
			}
			sorted[n] = val[nmax];
			val[nmax] = INF;
		}
		for (int n=0;n<val.length;n++) {
			val[n] = sorted[n];
		}
		return;
	}	
	
	public static final void sort(float[] val) {
		for (int n=0;n<val.length;n++) {
			for (int m=n+1;m<val.length;m++) {
				if (val[m]<val[n]) {
					// switch place
					float tmp = val[n];
					val[n] = val[m];
					val[m] = tmp;
				}
			}
		}
	}
			
	public static final void sortmag(float[] val, int size) {
		for (int n=0;n<size;n++) {
			for (int m=n+1;m<size;m++) {
				if (val[m]*val[m]<val[n]*val[n]) {
					// switch place
					float tmp = val[n];
					val[n] = val[m];
					val[m] = tmp;
				}
			}
		}
	}
			
	public static final float avg(float[] val) {
		float avg=0.0f;
		for (int n=0;n<val.length;n++) avg += val[n];
		return avg/val.length;
	}
	
	public static final float sum(float[] val) {
		float sum=0.0f;
		for (int n=0;n<val.length;n++) sum += val[n];
		return sum;
	}
	
	public static final float stdev(float[] val) {
		float avg=0.0f;
		float stdev=0.0f;
		for (int n=0;n<val.length;n++) avg += val[n];
		avg = avg/val.length;
		for (int n=0;n<val.length;n++) stdev += square(val[n]-avg);
		stdev = sqrt(stdev/val.length);
		return stdev;
	}
	
	/** gives the position of the highest value*/
	public static final int argmax(float a, float b, float c) {
		if (a==b && b==c) return -1;
		if (a>b) {
			if (a>c) return 0;
			else return 2;
		} else { 
			if (b>c) return 1;
			else return 2;
		}
	}	
	/** gives the position of the highest value*/
	public static final int argmax(float a, float b) {
		if (a==b) return -1;
		else if (a>b) return 0;
		else return 1;
	}	
			
	/** gives the position of the highest absolute value*/
	public static final int argmaxmag(double a, double b, double c) {
		//if (a==b && b==c) return -1;
		if (a*a>b*b) {
			if (a*a>c*c) return 0;
			else return 2;
		} else { 
			if (b*b>c*c) return 1;
			else return 2;
		}
	}	
			
	/** gives the position of the second highest value*/
	public static final int argsec(double a, double b, double c) {
		if (a>b) {
			if (b>c) return 1;
			else if (a>c) return 2;
			else return 0;
		} else { 
			if (a>c) return 0;
			else if (b>c) return 2;
			else return 1;
		}
	}	
			
	/** gives the position of the second highest absolute value*/
	public static final int argsecmag(double a, double b, double c) {
		//if (a==b && b==c) return -1;
		if (a*a>b*b) {
			if (b*b>c*c) return 1;
			else if (a*a>c*c) return 2;
			else return 0;
		} else { 
			if (a*a>c*c) return 0;
			else if (b*b>c*c) return 2;
			else return 1;
		}
	}	
			
	/** gives the position of the lowest absolute value*/
	public static final int argminmag(double a, double b, double c) {
		if (a*a<b*b) {
			if (a*a<c*c) return 0;
			else return 2;
		} else { 
			if (b*b<c*c) return 1;
			else return 2;
		}
	}	
			
	/** gives the position of the highest value*/
	public static final int argmax(double a, double b, double c) {
		if (a==b && b==c) return -1;
		if (a>b) {
			if (a>c) return 0;
			else return 2;
		} else { 
			if (b>c) return 1;
			else return 2;
		}
	}	
			
	/** gives the position of the highest value*/
	public static final int argmax(float a, float b, float c, float d) {
		if (a==b && b==c & c==d) return -1;
		if (a>b) {
			if (a>c) {
				if (a>d) return 0;
				else return 3;
			} else {
				if (c>d) return 2;
				else return 3;
			}
		} else { 
			if (b>c) {
				if (b>d) return 1;
				else return 3;
			} else {
				if (c>d) return 2;
				else return 3;
			}
		}
	}	
			
	/** boolean to number value */
	public static final byte bin(boolean val) {
		if (val) return 1;
		else return 0;
	}
	
	/** squaring */
	public static final float square(float val) {
		return val*val;
	}
	public static final double square(double val) {
		return val*val;
	}
	
	/** cubing */
	public static final float cube(float val) {
		return val*val*val;
	}
	public static final double cube(double val) {
		return val*val*val;
	}
	public static final int cube(int val) {
		return val*val*val;
	}
	
	/** quading */
	public static final float quad(float val) {
		return val*val*val*val;
	}
	public static final double quad(double val) {
		return val*val*val*val;
	}
	
	/** boolean into 0/1 switch */
	public static final float kronecker(boolean val) {
		if (val) return 1.0f;
		else return 0.0f;
	}
	
	/** encapsulated functions for floats */
	public static final float exp(float val) {
		return (float)Math.exp(val);
	}
	public static final float log(float val) {
		return (float)Math.log(val);
	}
	public static final float sqrt(float val) {
		return (float)Math.sqrt(val);
	}
	
	public static final double modulo(double val, double mod) {
	    return val-Numerics.round(val/mod)*mod;
	    /*
	    if (val>mod/2.0) return val-mod;  
	    else if (val<-mod/2.0) return val+mod;
	    else return val;
	    */
	}
	public static final float modulo(float val, float mod) {
	    return val-Numerics.round(val/mod)*mod;
	    /*
	    if (val>mod/2.0) return val-mod;  
	    else if (val<-mod/2.0) return val+mod;
	    else return val;
	    */
	}
	public static final int wrap(float val, float mod) {
	    return Numerics.round(val/mod);
	}
}
