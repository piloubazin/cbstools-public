package edu.jhu.ece.iacl.algorithms.graphics.smooth;

import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import edu.jhu.ece.iacl.algorithms.VersionUtil;
import edu.jhu.ece.iacl.algorithms.graphics.XuGraphicsWrapper;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.DistanceFieldOnSurface;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface.Edge;

import Jama.*;

/**
 * Inflate surface to specified mean curvature using Lorentzian norm as
 * termination criterion. A weighted umbrella operator is used to regularize the
 * surface. This method is useful for generating surfaces that are at multiple
 * levels of geometric complexity, although cusp formation is common with this
 * inflation procedure.
 * 
 * @author Blake Lucas
 * 
 */
public class SurfaceInflate extends XuGraphicsWrapper {
	public static String getVersion() {
		return VersionUtil.parseRevisionNumber("$Revision: 1.13 $");
	}

	public SurfaceInflate(AbstractCalculation parent, EmbeddedSurface surf) {
		super(parent, surf.clone(), true);
		setLabel("Surface Inflation");
		this.surf.setName(surf.getName() + "_inf");
	}

	public SurfaceInflate(EmbeddedSurface surf) {
		super(surf.clone(), true);
		setLabel("Surface Inflation");
		this.surf.setName(surf.getName() + "_inf");
	}

	public void inflate() {
		inflate(0.75, 3, 10, 2000, true);
	}

	protected class Vertex {
		Vector3f v;
		int N;
		int[] nbd;
	}
	/**
	 * Inflate surface
	 * @param BETA step size
	 * @param step number of iterations at which to rescale surface to preserve volume
	 * @param EPS curvature threshold for termination
	 * @param imax maximum number of iterations
	 * @param useLorentzian true if lorentzian scaling should be used
	 * @return inflated surface
	 */
	public EmbeddedSurface inflate(double BETA, int step, float EPS, int imax,
			boolean useLorentzian) {
		float scaleL = 0;
		int k, j;
		int N = 1, MR_N = 1; // number of multi-resolution steps
		float MR_step = 2; // difference on the L2 norm of mean curvature
							// between different multiresolution levels
		setTotalUnits(100);
		int Vnum = XUGetVertexTableVN();
		System.out.printf("mean curvature measure threshold = %.2f\n", EPS);
		if (MR_N > 1)
			System.out.printf(
					"delta_measure between multiresolution levels = %g\n",
					MR_step);

		double areaO = 0, value;
		int PN = XUGetPolyTablePN();
		int i;
		int[] Mask = new int[Vnum];
		for (i = 0; i < PN; ++i)
			areaO += XUPolyArea(i);
		if (useLorentzian) {
			if (scaleL == 0)
				scaleL = (float) updateScaleL(Mask, areaO, Vnum);
			System.out.printf("Lorentzian scale factor = %g\n", scaleL);
			value = Math.sqrt(Math.log(1 + 0.5 / (scaleL * scaleL)));
			System.out
					.printf(
							"... the mean curvature measure threshold will be scaled accordingly\n... i.e., %g*%g=%g\n",
							EPS, value, EPS * value);
			EPS *= value;
			if (MR_N > 1) {
				System.out
						.printf(
								"... the delta_measure between multires. levels will be scaled accordingly\n... i.e., %g*%g=%g\n",
								MR_step, value, MR_step * value);
				MR_step *= value;
			}
		}

		System.out.printf("\n");
		float max;
		Vector3f C;
		double hCurv, hCurv_0, hCurvOld;
		double origVolume = surf.volume();
		// create nbd list for each vertex --- won't be changed
		Vertex[] VList = new Vertex[Vnum];

		for (i = 0; i < Vnum; ++i) {
			VList[i] = new Vertex();
			VList[i].N = XUGetNeighborEN(i);
			VList[i].nbd = new int[VList[i].N];
			for (k = 0; k < VList[i].N; k++) {
				VList[i].nbd[k] = XUGetNeighborPolyId(i, k);
			}
		}
		int pid;
		Vector3f pnt = new Vector3f(), V;
		double area, totalA, areaN, scale;
		hCurv = meanCurv(Mask, areaO, Vnum, useLorentzian, scaleL);
		hCurv_0 = hCurv;
		hCurvOld = hCurv;
		j = 1;
		double hCurvFirst = -1;
		do {
			for (i = 0; i < Vnum; ++i) {
				totalA = 0;
				pnt.x = 0;
				pnt.y = 0;
				pnt.z = 0;
				// Current vertex
				V = XUGetVertex(i);

				// Regularize mesh by small amount BETA by moving point to
				// center
				for (k = 0; k < VList[i].N; ++k) {
					pid = VList[i].nbd[k];
					// Calculate center of surrounding region
					C = XUPolyCenter(pid);
					area = XUPolyArea(pid);
					totalA += area;

					pnt.x += area * C.x;
					pnt.y += area * C.y;
					pnt.z += area * C.z;
				}
				if (totalA > 0) {
					pnt.x /= totalA;
					pnt.y /= totalA;
					pnt.z /= totalA;
				}
				V.x = (float) (V.x * (1 - BETA) + BETA * pnt.x);
				V.y = (float) (V.y * (1 - BETA) + BETA * pnt.y);
				V.z = (float) (V.z * (1 - BETA) + BETA * pnt.z);
				// System.out.println("POINT "+i+" "+V);
				surf.setVertex(i, V);
			}

			if (j % step == 0) {
				// Rescale vertices to preserve volume
				scale = Math.pow((origVolume / surf.volume()), 0.3333333333);
				Point3f centroid = surf.getCenterOfMass();
				for (i = 0; i < Vnum; ++i) {
					V = XUGetVertex(i);
					V.x = (float) ((V.x - centroid.x) * scale + centroid.x);
					V.y = (float) ((V.y - centroid.y) * scale + centroid.y);
					V.z = (float) ((V.z - centroid.z) * scale + centroid.z);
					surf.setVertex(i, V);
				}
				areaN = 0;
				for (i = 0; i < PN; ++i)
					areaN += XUPolyArea(i);
				/*** curvature calculation ***/
				// Calculate Total Mean Curvature
				hCurv = meanCurv(Mask, areaN, Vnum, useLorentzian, scaleL);
				max = (float) (Math.abs(hCurv - hCurvOld) / hCurv_0);
				hCurvOld = hCurv;

				if ((hCurv < (EPS + (N - 1) * MR_step)) || (j == (imax))) {
					// Compute Shape Not Used
					N -= 1;
				}
				if (hCurvFirst == -1)
					hCurvFirst = hCurv;
				setCompletedUnits(1 - (hCurv - EPS) / (hCurvFirst - EPS));
				System.out.println("Mean Curvature " + hCurv);
			}
			++j;
		} while ((hCurv > EPS) && (j <= imax));
		markCompleted();
		return surf;
	}

	/********************************************************/
	protected double meanCurv(int[] mask, double area, int VN,
			boolean useLorentzian, float scaleL) {
		int i, j;
		Vector3f v0, v1, v2;
		Vector3f v, w, u;
		int vid, pid, EN;
		double curv;
		Vector3f Curv = new Vector3f();
		double alpha, temp = 0, Narea;
		double kcurv;

		curv = 0;
		for (i = 0; i < VN; ++i) {
			if (mask[i] == 0) {
				v0 = XUGetVertex(i);
				Curv.x = 0;
				Curv.y = 0;
				Curv.z = 0;
				EN = XUGetNeighborEN(i);
				vid = XUGetVertexNeighborVId(i, EN - 1);
				v1 = XUGetVertex(vid);
				Narea = 0.0;
				for (j = 0; j < EN; ++j) {
					pid = XUGetNeighborPolyId(i, j);
					Narea += XUPolyArea(pid);
					vid = XUGetVertexNeighborVId(i, j);
					v2 = XUGetVertex(vid);

					v = XUVecSub(v0, v2);
					w = XUVecSub(v1, v2);
					alpha = Math.acos(Math.max(-1, Math.min(1, XUVecDot(v, w)
							/ (XUVecMag(v) * XUVecMag(w)))));
					if (alpha != 0) {
						temp = 1.0f / Math.tan(alpha);
					} else {
						temp = 0;
						System.err.println("Alpha is " + alpha + " " + v + " "
								+ w);
					}
					vid = XUGetVertexNeighborVId(i, ((EN + j - 2) % EN));
					v2 = XUGetVertex(vid);

					v = XUVecSub(v0, v2);
					w = XUVecSub(v1, v2);
					alpha = Math.acos(Math.max(-1, Math.min(1, XUVecDot(v, w)
							/ (XUVecMag(v) * XUVecMag(w)))));
					if (alpha != 0) {
						temp += 1 / Math.tan(alpha);
					} else {
						System.err.println("Alpha is " + alpha + " " + v + " "
								+ w);
					}
					if (Double.isNaN(temp) || Double.isInfinite(temp)) {
						temp = 0;
					}
					v = XUVecSub(v0, v1);

					Curv.x += v.x * temp * 0.5f;
					Curv.y += v.y * temp * 0.5f;
					Curv.z += v.z * temp * 0.5f;
					vid = XUGetVertexNeighborVId(i, j);
					v1 = XUGetVertex(vid);
				}
				Narea /= 3.0;
				kcurv = (Narea > 0) ? Curv.length() / Narea : 0;
				if (useLorentzian)
					curv += Math.log(1 + 0.5 * kcurv * kcurv
							/ (scaleL * scaleL))
							* Narea;
				else
					curv += kcurv * kcurv * Narea;
			}
		}
		// curv /= area;
		curv *= curvNorm;
		curv = Math.sqrt(curv);
		return (curv);
	}

	private static final float curvNorm = (float) (1.0f / (4 * Math.acos(-1)));

	protected double updateScaleL(int[] mask, double area, int VN) {
		int i, j;
		Vector3f v0, v1, v2;
		Vector3f v, w;
		int vid, pid, EN;
		Vector3f Curv = new Vector3f();
		double alpha, temp, Narea;
		float scaleL = 0;
		float median;

		float[] kcurv = new float[VN];
		;
		for (i = 0; i < VN; ++i) {
			if (mask[i] == 0) {
				v0 = XUGetVertex(i);
				Curv.x = 0;
				Curv.y = 0;
				Curv.z = 0;
				EN = XUGetNeighborEN(i);
				vid = XUGetVertexNeighborVId(i, EN - 1);
				v1 = XUGetVertex(vid);
				Narea = 0.0;
				for (j = 0; j < EN; ++j) {
					pid = XUGetNeighborPolyId(i, j);
					Narea += XUPolyArea(pid);
					vid = XUGetVertexNeighborVId(i, j);
					v2 = XUGetVertex(vid);
					v = XUVecSub(v0, v2);
					w = XUVecSub(v1, v2);
					alpha = Math.acos(XUVecDot(v, w)
							/ (XUVecMag(v) * XUVecMag(w)));
					temp = 1 / Math.tan(alpha);
					vid = XUGetVertexNeighborVId(i, ((EN + j - 2) % EN));
					v2 = XUGetVertex(vid);
					v = XUVecSub(v0, v2);
					w = XUVecSub(v1, v2);
					alpha = Math.acos(XUVecDot(v, w)
							/ (XUVecMag(v) * XUVecMag(w)));
					if (alpha != 0) {
						temp += 1 / Math.tan(alpha);
					}
					if (Double.isNaN(temp) || Double.isInfinite(temp)) {
						temp = 0;
					}
					v = XUVecSub(v0, v1);
					Curv.x += v.x * temp / 2;
					Curv.y += v.y * temp / 2;
					Curv.z += v.z * temp / 2;
					vid = XUGetVertexNeighborVId(i, j);
					v1 = XUGetVertex(vid);
				}
				Narea /= 3.0f;
				kcurv[i] = (Narea > 0) ? (float) (Curv.length() / Narea) : 0;
			}
		}
		hpsort(VN, kcurv);
		median = (kcurv[(int) Math.floor(VN / 2)] + kcurv[(int) Math
				.ceil(VN / 2)]) / 2;
		// printf("median = %f ---> ",median);
		for (i = 0; i < VN; ++i)
			kcurv[i] = Math.abs(kcurv[i] - median);
		hpsort(VN, kcurv);
		median = (kcurv[(int) Math.floor(VN / 2)] + kcurv[(int) Math
				.ceil(VN / 2)]) / 2;
		// printf("%f\n",median);
		// if (median > 1)
		// Magic number for Median Absolute Deviation MAD
		scaleL = MedianAbsoluteDeviation * median;
		// else
		// scaleL = 1.4826;
		// printf("scaleL = %f\n",scaleL);
		return (scaleL);
	}

	protected static final float MedianAbsoluteDeviation = 1.4826f;

	void hpsort(int n, float[] r) {
		int i, ir, j, l;
		float rra;
		int k;

		float[] ra;
		int[] x;
		int xx;

		x = new int[n + 1];
		ra = new float[n + 1];
		for (k = 1; k <= n; ++k) {
			x[k] = k - 1;
			ra[k] = r[k - 1];
		}

		if (n < 2) {
			System.out.printf("size of the array is < 2! \n");
			System.exit(1);
		}
		l = (n >> 1) + 1;
		ir = n;

		/*
		 * * The index l will be decremented from its initial value down to 1
		 * during the "hiring"* (heap creation) phase. Once it reaches 1, the
		 * index ir will be decremented from its* initial value down to 1 during
		 * the "retirement-and-promotion" (heap selection) phase.
		 */

		for (;;) {
			if (l > 1) { /* still in hiring phase */
				k = --l;
				xx = x[k];
				rra = ra[k];
			} else { /* In retirment-and-promotion phase. */
				rra = ra[ir]; /* Clear a space at the end of array. */
				xx = x[ir];
				ra[ir] = ra[1]; /* Retire the top of the array. */
				x[ir] = x[1];
				if (--ir == 0) { /* Done with the last promotion. */
					ra[1] = rra; /* The last competent worker of all! */
					break;
				}
			}
			i = l; /* Whether in hiring phase or promotion phase, */
			j = l + l; /*
						 * here set up to shift down element rra to its proper
						 * level.
						 */
			while (j <= ir) {
				if (j < ir && ra[j] < ra[j + 1])
					j++; /* Compare to the better underling. */
				if (rra < ra[j]) { /* Demote rra. */
					ra[i] = ra[j];
					x[i] = x[j];
					i = j;
					j <<= 1;
				} else
					break; /* Found rra's level. Terminate the sift-down. */
			}
			ra[i] = rra; /* Put rra into slot. */
			x[i] = xx;
		}

		for (k = 0; k < n; ++k) {
			r[k] = ra[k + 1];
		}
	}

	int NBR_MAX_NUM = 5000;

}
