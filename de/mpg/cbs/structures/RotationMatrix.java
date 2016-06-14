package de.mpg.cbs.structures;

import java.io.*;
import java.util.*;

/**
 *
 *  This class handles 3x3 rotation matrices with
 *	inverse Rodrigues parametrization.
 *  <p>
 *	No matrix class encapsulation is needed, 
 *  as most matrix computations are simplified with rotations.
 *	
 *	@version    November 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class RotationMatrix {
	
	private float[][]	matrix;
	private float[]		parameters;
	
	public RotationMatrix() {
		matrix = new float[3][3];
		parameters = new float[3];
        
        for (int i=0;i<3;i++) {
            for (int j=0;j<3;j++) matrix[i][j] = 0.0f;
            matrix[i][i] = 1.0f;
            parameters[i] = 0.0f;
        }       
	}
    
    public float getParameter(int i) {
        return parameters[i];
    }
    
    public float[] getParameters() {
        return parameters;
    }
    
    public float getMatrix(int i,int j) {
        return matrix[i][j];
    }
        
    public float[][] getMatrix() {
        return matrix;
    }
        
    public void setParameters(float a, float b, float c) {
        double x, a2, b2, c2;
        double an, bn, cn, n, m;

		// reproject in appropriate space
        n = Math.sqrt( a*a + b*b + c*c );
        if (n > 1.0) {
			// reproject onto [-1,1] instead
            m = n % 2;
			if (m>1.0) {
				an = m * a / n - 2 * a / n;
				bn = m * b / n - 2 * b / n;
				cn = m * c / n - 2 * c / n;
			} else {
				an = a * m / n;
				bn = b * m / n;
				cn = c * m / n;
			}
        } else {
            an = a;
            bn = b;
            cn = c;
        }

        // set the parameters
        parameters[0] = (float)an; parameters[1] = (float)bn; parameters[2] = (float)cn;
        
        // compute the new matrix
        a2 = an*an; b2 = bn*bn; c2 = cn*cn;

		x = 1.0 - a2 - b2 - c2;
		// for roundoff errors ?
		if ( x<0.0 && x>-1e-10) x = 0.0;
        x = Math.sqrt(x);

		if (Double.isNaN(x)) System.out.print("! x^2 < 0: "+a2+", "+b2+", "+c2+"\n");
		
        matrix[0][0] = (float)(1.0 - 2.0 * b2 - 2.0 * c2);
        matrix[1][0] = (float)(2.0 * an * bn + 2.0 * x * cn);
        matrix[2][0] = (float)(2.0 * an * cn - 2.0 * x * bn);
        matrix[0][1] = (float)(2.0 * an * bn - 2.0 * x * cn);
        matrix[1][1] = (float)(1.0 - 2.0 * a2 - 2.0 * c2);
        matrix[2][1] = (float)(2.0 * bn * cn + 2.0 * x * an);
        matrix[0][2] = (float)(2.0 * an * cn + 2.0 * x * bn);
        matrix[1][2] = (float)(2.0 * bn * cn - 2.0 * x * an);
        matrix[2][2] = (float)(1.0 - 2.0 * a2 - 2.0 * b2);

        return;
    }
    public void setParameters(float[] param) {
        setParameters(param[0],param[1],param[2]);
    }
    public void setParameters(float[] param, int offset) {
        setParameters(param[offset+0],param[offset+1],param[offset+2]);
    }
    
    public void setMatrix(float[][] mat) {
        float x, a2, b2, c2;

        // set the matrix
        for (int i=0;i<3;i++) for (int j=0;j<3;j++)
            matrix[i][j] = mat[i][j];
        
        // compute the new parameters
        x = (float)Math.sqrt(mat[0][0] + mat[1][1] + mat[2][2] + 1.0f)/ 2.0f;

        if (x > 0.0) {
            parameters[0] = (mat[2][1] / x - mat[1][2] / x) / 4.0f;
            parameters[1] = (mat[0][2] / x - mat[2][0] / x) / 4.0f;
            parameters[2] = (mat[1][0] / x - mat[0][1] / x) / 4.0f;
        } else {
            if ((mat[0][0] < 1.0) && (mat[1][1] < 1.0) && (mat[2][2] < 1.0)) {
                parameters[0] = (float)Math.sqrt( (mat[0][1]*mat[0][1] + mat[0][2]*mat[0][2]) 
                                                    / (2.0f*(1.0 - mat[0][0])) );
                parameters[1] = (float)Math.sqrt( (mat[1][0]*mat[1][0] + mat[1][2]*mat[1][2]) 
                                                    / (2.0f*(1.0 - mat[1][1])) );
                parameters[2] = (float)Math.sqrt( (mat[2][0]*mat[2][0] + mat[2][1]*mat[2][1]) 
                                                    / (2.0f*(1.0 - mat[2][2])) );
                if (mat[0][1] < 0.0) parameters[1] = -parameters[1];
                if (mat[0][2] < 0.0) parameters[2] = -parameters[2];
            } else if (mat[0][0] == 1.0f) {
                parameters[0] = 1.0f;
                parameters[1] = 0.0f;
                parameters[2] = 0.0f;
            } else if (mat[1][1] == 1.0f) {
                parameters[0] = 0.0f;
                parameters[1] = 1.0f;
                parameters[2] = 0.0f;
            } else if (mat[2][2] == 1.0f) {
                parameters[0] = 0.0f;
                parameters[1] = 0.0f;
                parameters[2] = 1.0f;
            } 
        }
    }

    public void setMatrix(double[][] mat) {
        float x, a2, b2, c2;

        // set the matrix
        for (int i=0;i<3;i++) for (int j=0;j<3;j++)
            matrix[i][j] = (float)mat[i][j];
        
        // compute the new parameters
        x = (float)Math.sqrt(mat[0][0] + mat[1][1] + mat[2][2] + 1.0f)/ 2.0f;

        if (x > 0.0) {
            parameters[0] = (float) ( (mat[2][1] / x - mat[1][2] / x) / 4.0 );
            parameters[1] = (float) ( (mat[0][2] / x - mat[2][0] / x) / 4.0 );
            parameters[2] = (float) ( (mat[1][0] / x - mat[0][1] / x) / 4.0 );
        } else {
            if ((mat[0][0] < 1.0) && (mat[1][1] < 1.0) && (mat[2][2] < 1.0)) {
                parameters[0] = (float)Math.sqrt( (mat[0][1]*mat[0][1] + mat[0][2]*mat[0][2]) 
                                                    / (2.0f*(1.0 - mat[0][0])) );
                parameters[1] = (float)Math.sqrt( (mat[1][0]*mat[1][0] + mat[1][2]*mat[1][2]) 
                                                    / (2.0f*(1.0 - mat[1][1])) );
                parameters[2] = (float)Math.sqrt( (mat[2][0]*mat[2][0] + mat[2][1]*mat[2][1]) 
                                                    / (2.0f*(1.0 - mat[2][2])) );
                if (mat[0][1] < 0.0) parameters[1] = -parameters[1];
                if (mat[0][2] < 0.0) parameters[2] = -parameters[2];
            } else if (mat[0][0] == 1.0f) {
                parameters[0] = 1.0f;
                parameters[1] = 0.0f;
                parameters[2] = 0.0f;
            } else if (mat[1][1] == 1.0f) {
                parameters[0] = 0.0f;
                parameters[1] = 1.0f;
                parameters[2] = 0.0f;
            } else if (mat[2][2] == 1.0f) {
                parameters[0] = 0.0f;
                parameters[1] = 0.0f;
                parameters[2] = 1.0f;
            } 
        }
    }

    public void derivatives(float[][] dR, float da, float db, float dc) {
        float x, x2, a2, b2, c2, ab, ac, bc, ax, bx, cx;
        float an, bn, cn, n, m;
        float ZERO = 1e-6f;

        n = (float)Math.sqrt(parameters[0]*parameters[0] + parameters[1]*parameters[1] + parameters[2]*parameters[2]);
        if (n > 1.0) {
			// reproject onto [-1,1] instead ?
			m = n % 2;
			if (m>1.0) {
				an = m * parameters[0] / n - 2 * parameters[0] / n;
				bn = m * parameters[1] / n - 2 * parameters[1] / n;
				cn = m * parameters[2] / n - 2 * parameters[2] / n;
			} else {
				an = parameters[0] * m / n;
				bn = parameters[1] * m / n;
				cn = parameters[2] * m / n;
			}
         } else {
            an = parameters[0];
            bn = parameters[1];
            cn = parameters[2];
        }

        a2 = an*an; b2 = bn*bn; c2 = cn*cn;
        ab = an*bn; ac = an*cn; bc = bn*cn;
        x2 = 1.0f - a2 - b2 - c2;
        x = (float)Math.sqrt(x2);
        ax = x*an; bx = x*bn; cx = x*cn;

        /* should not occur: if x=0 then a=b=c=0 then x=1... */
        if (x < ZERO) System.out.println("x<0 in rotation");

        dR[0][0] = -4.0f * (bn * db + cn * dc);
        dR[1][0] = 2.0f / x * ((bx - ac) * da + (ax - bc) * db + (x2 - c2) * dc);
        dR[2][0] = 2.0f / x * ((cx + ab) * da + (b2 - x2) * db + (ax + bc) * dc);
        dR[0][1] = 2.0f / x * ((bx + ac) * da + (ax + bc) * db + (c2 - x2) * dc);
        dR[1][1] = -4.0f * (an * da + cn * dc);
        dR[2][1] = 2.0f / x * ((x2 - a2) * da + (cx - ab) * db + (bx - ac) * dc);
        dR[0][2] = 2.0f / x * ((cx - ab) * da + (x2 - b2) * db + (ax - bc) * dc);
        dR[1][2] = 2.0f / x * ((a2 - x2) * da + (cx + ab) * db + (bx + ac) * dc);
        dR[2][2] = -4.0f * (an * da + bn * db);

        return;
    }
    public float[][] derivatives(float da, float db, float dc) {
        float[][] dR = new float[3][3];
        derivatives(dR, da, db, dc);
        return dR;
    }
    
	
	public float[] computeEulerAngles() {
		float[] euler = new float[3];
		double cphi,sphi,ctheta,stheta,cpsi,spsi;
		
		ctheta = matrix[2][2];
		stheta = Math.sqrt(1.0f-ctheta*ctheta);
		
		if (stheta>0) {
			cphi = -matrix[2][1]/stheta;
			sphi =  matrix[2][0]/stheta;
			
			cpsi =  matrix[1][2]/stheta;
			spsi =  matrix[0][2]/stheta;
			
			euler[0] = (float)Math.atan2(sphi,cphi);
			euler[1] = (float)Math.atan2(stheta,ctheta);
			euler[2] = (float)Math.atan2(spsi,cpsi);
			
		} else {
			cphi = matrix[0][0];
			sphi = matrix[0][1];
			
			euler[0] = (float)Math.atan2(sphi,cphi);
			euler[1] = 0.0f;
			euler[2] = 0.0f;
		}
		
        return euler;
    }
        

}
