package de.mpg.cbs.utilities;


/**
 *
 *  This class computes various basic numerical functions.
 *	<p> 
 *
 *	@version    May 17, 2005
 *	@author     Pierre-Louis Bazin
 */

public class ProcessMatrix {
	
	/** mathematical constants */
	public static final float ZERO = 1e-30f;
	/** mathematical constants */
	public static final float INF = 1e+30f;
	/** mathematical constants */
	public static final float PI = 3.1416f;
	
	/**
	 * This method compute the inverse of a matrix 
	 * @param m : the given matrix which should be a square array of floats 
	 * @param inverse : is the result of inversion
	 * @param nColumn : number of column of squared array
	 * @return : determinant of the matrix
	 */
	
	public static final float invertMatrix(float m[][], float[][] inverse, int nColumn){
		
		float det = 0.0f;
		//float[][] inverse = new float[nColumn][nColumn];
		
		if (nColumn == 1) {
			det = m[0][0];
			if (det != 0.0f)
				inverse[0][0] = 1.0f/m[0][0];
			else
				System.out.println("Matrix is singular");
			
			return det;
			//return inverse;
		} else if (nColumn == 2) {
			det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
			if (det != 0.0f) {
				inverse [0][0] =m[1][1]/det;
				inverse [0][1] = -m[1][0]/det;
				inverse [1][0] = -m[0][1]/det;
				inverse [1][1] = m[1][1]/det;
				
			} else
				System.out.println("Matrix is singular");
			return det;
			//return inverse;
		} else if (nColumn == 3) {
			det = m[0][0]*m[1][1]*m[2][2] + m[1][0]*m[2][1]*m[0][2] + m[2][0]*m[0][1]*m[1][2]
			    - m[2][0]*m[1][1]*m[0][2] - m[2][1]*m[1][2]*m[0][0] - m[2][2]*m[1][0]*m[0][1];
			if (det != 0.0f) {
				inverse[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1])/det;
				inverse[0][1]= (m[2][1]*m[0][2]-m[2][2]*m[0][1])/det;
				inverse[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1])/det;
				inverse[1][0] = (m[1][2]*m[2][0]-m[1][0]*m[2][2])/det;
				inverse[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0])/det;
				inverse[1][2] = (m[1][0]*m[0][2]-m[1][2]*m[0][0])/det;
				inverse[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0])/det;
				inverse[2][1]= (m[0][1]*m[2][0]-m[0][0]*m[2][1])/det;
				inverse[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0])/det;
			} else
				System.out.println("Matrix is singular");
			return det;
			//return inverse;
		} else {
			int[] pivlst = new int[2*nColumn+1];
			boolean[] pivchk = new boolean[nColumn];
			int leRow,leCol;
			leRow = leCol =0;
			float piv,t,leval,temp;
			int i,j,k; //loop counters
			int zeroMat = -1;
			det = 1.0f;
			for (i =0; i<nColumn ;i++ ){
				pivchk[i] = false;
				for (j =0; j<nColumn; j++)
					inverse[i][j] = m[i][j];
			}
			
			for (i=0; i<nColumn; i++) {
				leval = 0.0f;
				for ( j = 0; j < nColumn; j++ ) 
		            if ( ! (pivchk[j]) ) 
		               for ( k = 0; k < nColumn; k++ ) 
		                  if ( ! (pivchk[k]) ) {
							 temp = Math.abs(inverse[j][k]);
		                     if ( temp > leval ) {
		                        leRow = j;
		                        leCol = k;
		                        leval = temp;
								zeroMat=0;
		                     }
		                  }
				//System.out.println("leCol= "+ leCol + "leRow= "+leRow);
				if(zeroMat != 0) {
					System.out.println("Cannot invert the zero matrix\n");
					det=0.0f;
					inverse = null;
				 }
		         pivchk[leCol] = true;
		         pivlst[i*2] = leRow;
		         pivlst[i*2+1] = leCol;
		         if ( leRow != leCol ) {
		            det = -det;
		            for ( int l = 0; l < nColumn; l++ ) {
		               float swap = inverse[leRow][l];
		               inverse[leRow][l] = inverse[leCol][l];
		               inverse[leCol][l] = swap;
		            }
		         }
		         piv = inverse[leCol][leCol];
		         det = det * piv;
		         //System.out.println(det);
		         if ( det > 1.0e+30 ) {
		            det = 1.0f;
		         }
		         inverse[leCol][leCol] = 1.0f;
		         for ( int l = 0; l < nColumn; l++ ) {
		            inverse[leCol][l] = inverse[leCol][l] / piv;
		            //System.out.println("inverse["+leCol+"]["+l+"]= "+ inverse[leCol][l]);
		         }
		         for ( int l1 = 0; l1 < nColumn; l1++ ) {
		            if ( l1 != leCol ) {
		               t = inverse[l1][leCol];
		               inverse[l1][leCol] = 0.0f;
		               for ( int l = 0; l < nColumn; l++ ) {
		                  inverse[l1][l] = inverse[l1][l] - inverse[leCol][l] * t;
		               }
		            }
		         }
		      }
		      
			  for ( i = 0; i < nColumn; i++ ) {
		         int l = nColumn - i - 1;
		         if ( pivlst[l*2] != pivlst[l*2+1] ) {
		            leRow = pivlst[l*2];
		            leCol = pivlst[l*2+1];
		            for ( k = 0; k < nColumn; k++ ) {
		            	float swap = inverse[k][leRow];
			            inverse[k][leRow] = inverse[k][leCol];
			            inverse[k][leCol] = swap;
		            }
		         }
		      }
		      
		   }
		return det;
		
	}
	
}
