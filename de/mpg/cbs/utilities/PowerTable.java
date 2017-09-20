package de.mpg.cbs.utilities;

/**
 * PowerTable.java
 *
 *	@version    February 2007
 *	@author     Joseph Hennessey
 */

/**
 * convenience class to compute p^q quickly with p and q within certain bounds.
 * the values are pre-computed and then accessed in the computations.
 * non-covered cases use the regular power function instead.
 */
public class PowerTable{
    private double pmin;
    private double pmax;
    private double pstep;
    private double prange;
    private double prangesteps;
    private int intPrangesteps;

    private double qCached;

    private double[] cache;

	private final double	ZERO = 1e-30;
    private final double	INF = 1e+30;
	
	private final boolean 	debug = false;
    
     public double lookup(double p, double q)
    {
        double result = 0.0;
        if((p<pmin)||(p>pmax)){
        	if (debug) System.err.println(p+" out of range: valid range is " + pmax + " >= p >= " + pmin + " with a precision of " + pstep);
		result = Math.pow(p, q);

	  } else if(q != qCached){
        	if (debug) System.err.println(q+" value not cached: cached q value is " + qCached );
		result = Math.pow(p, q);
	  } else {
		result = intLookup(p);
	  }
	  return(result);
    }
	
	/** look-up of value with no check */
    public double lookup(double p) {
		return intLookup(p);
	}
	
     /** Creates a new instance of powertable*/
   public PowerTable(double inPmin, double inPmax, double inPstep, double inQcached){
		pmin = inPmin;
  		pmax = inPmax;
		pstep = inPstep;
		prange = pmax - pmin;
		qCached = inQcached;
		prangesteps = prange / pstep;
		intPrangesteps = (int) Math.ceil(prangesteps);
		intPrangesteps++;
		int count;
		cache = new double[intPrangesteps];
		for (count = 0; count <intPrangesteps; count ++){
			double tempP = pmin + (count * pstep);
			if (qCached<0) {
				if (tempP<ZERO) cache[count] = INF;
				else cache[count] = Math.pow(tempP, qCached);
			} else cache[count] = Math.pow(tempP, qCached);
		}	
	}    

   private double intLookup(double p){
         double tempP = p - pmin;
	 tempP = tempP / pstep;
	 //int tempIntValue = (int) Math.round(tempP);
	 int tempIntValue = (int) (tempP + 0.5);
	 double result = cache[tempIntValue];
         return(result);
   }
      
   /* for testing
    public static void main(String[] args) {
	
	 double tempPmin = 0.0;
	 double tempPmax = 65535.0;
	 double tempPstep = 0.01;
       double tempQ = 1.625;

	 long time0 = 0;
	 long time1 = 0;
	 long time2 = 0;
	 long time3 = 0;

	 long dtime0 = 0;
	 long dtime1 = 0;
	 long dtime2 = 0;

	 time0 = java.lang.System.nanoTime();

	 powertable temppowertable  = new powertable( tempPmin , tempPmax , tempPstep , tempQ );
        
	 time1 = java.lang.System.nanoTime();

	 dtime0 = time1 - time0;
	 System.out.println("dtime0=" + dtime0 );

	 int count;
       double tempResult = 0.0;
       double tempValue = 0.0;
       int reps = 256;
       int repCount;
       
       for (repCount =0;repCount<reps;repCount++){     
	 tempValue = tempPmin;
	 for (count = 0; count < 65536; count ++){
            tempValue += tempPstep ;
		tempResult = temppowertable.lookup(tempValue,  tempQ );
	 }
       }
          System.out.println("tempValue = " + tempValue + ", tempResult = " + tempResult);
	 time2 = java.lang.System.nanoTime();

	 dtime1 = time2 - time1;
	 System.out.println("dtime1=" + dtime1 );

         
       for (repCount =0;repCount<reps;repCount++){     
	 tempValue = tempPmin;
	 for (count = 0; count < 65536; count ++){
            tempValue += tempPstep ;
		tempResult = Math.pow(tempValue, tempQ);
	 }
       }
         System.out.println("tempValue = " + tempValue + ", tempResult = " + tempResult);
 	 time3 = java.lang.System.nanoTime();

	 dtime2 = time3 - time2;
	 System.out.println("dtime2=" + dtime2 );


   }
    */
}
