package edu.jhu.ece.iacl.algorithms.graphics.utilities;

import java.util.Hashtable;

import javax.vecmath.Vector3f;

import edu.jhu.ece.iacl.algorithms.VersionUtil;
import edu.jhu.ece.iacl.algorithms.graphics.GeometricUtilities;
import edu.jhu.ece.iacl.algorithms.graphics.XuGraphicsWrapper;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.data.BinaryMinHeap;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.geom.VertexIndexed;

/**
 * Compute fast-marching distance from a set of vertices that are considered to
 * be the zero level set. There is no inside or outside for surfaces, so all
 * distances should be positive.
 * 
 * @author Blake
 * 
 */
public class DistanceFieldOnSurface extends XuGraphicsWrapper {
	public static String getVersion() {
		return VersionUtil.parseRevisionNumber("$Revision: 1.14 $");
	}

	public DistanceFieldOnSurface(EmbeddedSurface surf) {
		super(surf, true);
	}
	
	// added by Clara 2011-09-13
	// specify table to init for the mesh
	public DistanceFieldOnSurface(EmbeddedSurface surf, int initTables) {
		super(surf, false);
		/*
		edges = EmbeddedSurface.buildEdgeTable(surf);
		Hashtable<Long, EmbeddedSurface.Edge> hash = EmbeddedSurface
				.buildEdgeHash(edges);
		faces = EmbeddedSurface.buildFaceTable(surf, edges, hash);

		neighborEdgeFaceTable = EmbeddedSurface.buildNeighborEdgeFaceTable(
				surf, faces, hash, edges);
		neighborVertexEdgeTable = EmbeddedSurface.buildNeighborVertexEdgeTable(
				neighborVertexVertexTable, edges, hash);
		neighborVertexFaceTable = EmbeddedSurface.buildNeighborVertexFaceTable(
				neighborVertexEdgeTable, neighborEdgeFaceTable);
        */
	}
	

	public DistanceFieldOnSurface(AbstractCalculation parent,
			EmbeddedSurface surf) {
		super(parent, surf, true);
	}

	public DistanceFieldOnSurface(XuGraphicsWrapper wrappedSurf) {
		super(wrappedSurf);
	}

	private static final byte ALIVE = 1;

	private static final byte NBAND = 2;

	private static final byte FAWAY = 3;
	private static double INFINITY = 1000000000;
	public double[] solve(int[] contour, double[] F) {
		return solve(contour, contour.length, F);
	}
	double[] constantSpeed=null;
	public double[] solve(int[] contour) {
		if(constantSpeed==null){
			constantSpeed=new double[surf.getVertexCount()];
			for(int i=0;i<constantSpeed.length;i++)constantSpeed[i]=1;
		}
		return solve(contour, contour.length,constantSpeed);
	}
	public double[] solve(int[] contour, int numinitvert, double[] F) {
		double[] T; /* final distances of each vertex to the contour */
		double tempvalue, newvalue;
		byte[] label; /* "Alive" => 1; "Narrow Band" => 2; "Far Away" => 3; */
		// int[] BackPointer; /* backpointer to the narrowband heap */
		// Xheap H; /* Narrrow band heap */
		VertexIndexed he;
		int heapindex;
		double currentF;

		int i, j, k, NN, NN1, NN2, VN, VId, neighVId, n0, n1;

		VN = XUGetVertexTableVN();
		T = new double[VN];
		label = new byte[VN];

		/* Initialize heap structure */
		// H = xhInitEmpty();
		/* Initialization for marching inwards and outwards simultaneously */
		for (i = 0; i < VN; i++) {
			label[i] = (byte) FAWAY; /* All points are labelled as FAR AWAY */
			T[i] = INFINITY;
		}

		for (i = 0; i < numinitvert; i++) {
			label[contour[i]] = (byte) ALIVE; /* all vertices on contour are alive */
			T[contour[i]] = 0; /* all vertices on contour have distance zero */
		}
		BinaryMinHeap heap = new BinaryMinHeap(VN, VN, 1, 1);
		/* Initialize NarrowBand Heap */
		for (i = 0; i < VN; i++) {
			if (label[i] == (byte) ALIVE) {

				/* put its neighbors into the NarrowBand */
				NN1 = XUGetNeighborEN(i);
				for (j = 0; j < NN1; j++) {
					neighVId = XUGetVertexNeighborVId(i, j);

					if (label[neighVId] == (byte) FAWAY) {
						/*
						 * compute the distance to this point and add it to the
						 * heap
						 */
						label[neighVId] = (byte) NBAND;
						/*
						 * Check all possible triangles surrounding vertex
						 * neighVId1
						 */
						/*
						 * Note: Only ALIVE points contribute to the distance
						 * computation
						 */
						NN2 = XUGetNeighborEN(neighVId);
						newvalue = INFINITY;
						currentF = (F == null) ? 1 : F[neighVId];

						for (k = 0; k < NN2; k++) {
							n0 = XUGetVertexNeighborVId(neighVId, k);
							n1 = XUGetVertexNeighborVId(neighVId, (k + 1) % NN2);
							tempvalue = ReComputev2(neighVId, n0, n1, T, label,
									currentF, heap);
							if (tempvalue < newvalue)
								newvalue = tempvalue;
						} /* end for */

						T[neighVId] = newvalue;
						VertexIndexed vox = new VertexIndexed(newvalue);
						vox.setPosition(neighVId);
						heap.add(vox);

					} /* end if */
				} /* end for */
			} /* end if */
		} /* end for */

		/* End of Initialization */

		/*
		 * Begin Fast Marching to get the unsigned distance function inwords and
		 * outwards simultaneously since points inside and outside the contour
		 * won't interfere with each other
		 */
		while (!heap.isEmpty()) { /* There are still points not yet accepted */
			he = (VertexIndexed) heap.remove(); /*
												 * Label the point with smallest
												 * value among all NarrowBand
												 * points as ALIVE
												 */
			VId = he.getRow();

			T[VId] = he.getValue();
			label[VId] = (byte) ALIVE;

			/* Update its neighbor */
			/*
			 * Put FARAWAY neighbors into NarrowBand, Recompute values at
			 * NarrowBand neighbors, Keep ALIVE (Accepted) neighbor unchanged
			 */
			NN1 = XUGetNeighborEN(VId);
			for (i = 0; i < NN1; i++) {
				neighVId = XUGetVertexNeighborVId(VId, i);

				/* Don't change ALIVE neighbors */
				if (label[neighVId] != (byte) ALIVE) {
					NN2 = XUGetNeighborEN(neighVId);
					newvalue = INFINITY;
					currentF = (F == null) ? 1 : F[neighVId];
					for (j = 0; j < NN2; j++) {
						n0 = XUGetVertexNeighborVId(neighVId, j);
						n1 = XUGetVertexNeighborVId(neighVId, (j + 1) % NN2);
						tempvalue = ReComputev2(neighVId, n0, n1, T, label,
								currentF, heap);
						if (tempvalue < newvalue)
							newvalue = tempvalue;
					} /* end for */

					/*
					 * If it was a FARAWAY point, add it to the NarrowBand Heap;
					 * otherwise, just update its value using the backpointer
					 */
					VertexIndexed vox = new VertexIndexed((float) newvalue);
					vox.setPosition(neighVId);
					if (label[neighVId] == (byte) NBAND)
						heap.change(neighVId, 0, 0, vox);
					else {
						heap.add(vox);
						label[neighVId] = (byte) NBAND;
					}

				} /* end if ALIVE */
			} /* end updating neighbors */

		} /* end of marching loop */

		// Bug fix to handle case where nodes with no downwind neighbors are not
		// assigned distances.

		for (VId = 0; VId < VN; VId++) {
			if (T[VId] == INFINITY) {
				NN1 = XUGetNeighborEN(VId);
				newvalue = INFINITY;
				for (i = 0; i < NN1; i++) {
					neighVId = XUGetVertexNeighborVId(VId, i);
					NN2 = XUGetNeighborEN(neighVId);
					currentF = (F == null) ? 1 : F[neighVId];
					for (j = 0; j < NN2; j++) {
						n0 = XUGetVertexNeighborVId(neighVId, j);
						n1 = XUGetVertexNeighborVId(neighVId, (j + 1) % NN2);
						tempvalue = ReComputev2(neighVId, n0, n1, T, label,
								currentF, heap);
						if (tempvalue < newvalue)
							newvalue = tempvalue;
					}
				}
				T[VId] = newvalue;
			}
		}
		return T;
	}

	protected Vector3f vectsubptr(Vector3f va, Vector3f vb)
	/* Compute va - vb */
	{
		Vector3f out = new Vector3f();

		out.x = va.x - vb.x;
		out.y = va.y - vb.y;
		out.z = va.z - vb.z;

		return out;
	}

	protected double compute_dist(Vector3f v0, Vector3f v1) {
		double x, y, z;

		x = v1.x - v0.x;
		y = v1.y - v0.y;
		z = v1.z - v0.z;

		// printf("in the fnct, dist is %f \n", sqrt(x*x+y*y+z*z) );

		return (Math.sqrt(x * x + y * y + z * z));
	}

	/* vIDa, vIDb & vIDc are the three vertices making up the triangle */
	/* mesh is the mesh to compute the distances on */
	/* T is the current transit time values */
	/* label is the label map */
	/* compute the new T for vIDc based on the current T values of vIDa & vIDb */
	/* Unfold the mesh whenever C is an obtuse angle */
	/* see Kimmel and Sethian for derivation of these formulae */
	protected double ReComputev2(int vIDc, int vIDa, int vIDb, double[] T,
			byte[] label, double Fc, BinaryMinHeap heap) {
		double u, a, b, c;
		Vector3f Va, Vb, Vc;
		byte la, lb;
		double Ta, Tb;
		double t1, t2, t;
		int tmpint;
		double tmpdouble;

		int tmpvIDa;
		double tmpa;
		double tmpTa;
		byte tmpla;

		/* Define some auxillary variables */
		int NN, neighVId, i, UNotID, UID = 0, P1ID, P2ID;

		/* SQUARES of lengthes (Captital Letters) */
		double P1A, P1B, P1C, P2A, P2B, P2C, UA, UB, UC, P1P2, P1U, P2U;
		double cos12A, sin12A, cos12B, sin12B, cos12C, sin12C, cos12U, sin12U;
		double AC, BC, AB;

		Vector3f P1, P2, U;
		Vector3f vP1U, vP2U;

		int iters;

		la = label[vIDa];
		lb = label[vIDb];

		if ((la != (byte) ALIVE) && (lb != (byte) ALIVE))
			return INFINITY;

		if (Fc == 0)
			return INFINITY;

		Ta = T[vIDa];
		Tb = T[vIDb];

		if (la == (byte) ALIVE) {
			if (Ta == INFINITY)
				return INFINITY;
		} else
			Ta = INFINITY;

		if (lb == (byte) ALIVE) {
			if (Tb == INFINITY)
				return INFINITY;
		} else
			Tb = INFINITY;

		Va = XUGetVertex(vIDa);
		Vb = XUGetVertex(vIDb);
		Vc = XUGetVertex(vIDc);

		a = compute_dist(Vb, Vc);
		b = compute_dist(Va, Vc);

		if (Ta > Tb) {
			tmpTa = Ta;
			Ta = Tb;
			Tb = tmpTa;

			tmpla = la;
			la = lb;
			lb = tmpla;

			tmpa = a;
			a = b;
			b = tmpa;

			tmpvIDa = vIDa;
			vIDa = vIDb;
			vIDb = tmpvIDa;

			Va = XUGetVertex(vIDa);
			Vb = XUGetVertex(vIDb);
		}

		c = compute_dist(Va, Vb);
		AC = b * b;
		BC = a * a;
		AB = c * c;
		if (AB < AC + BC) {/* Acute Triangle */
			if (la == ALIVE && lb == ALIVE) {
				t = ComputeTInAcutev2(Ta, Tb, a, b, c, Fc);
				return t;
			} else
				return (Math.min(Ta + b / Fc, Tb + a / Fc)); /*
															 * Infinity + sth =
															 * Infinity
															 */
		}

		/* Otherwise, perform unfolding */

		/* Initialization */
		P1ID = vIDa;
		P2ID = vIDb;
		UNotID = vIDc;

		P1A = 0;
		P1B = AB;
		P1C = AC;
		P2B = 0;
		P2A = P1B;
		P2C = BC;
		P1P2 = P1B;

		cos12A = 1;
		sin12A = 0;
		cos12B = 1;
		sin12B = 0;
		cos12C = (P1P2 + P2C - P1C) / (2 * Math.sqrt(P1P2 * P2C));
		sin12C = (1 - cos12C * cos12C); /* Notice: Square of sine */

		/* Now iteratively unfolding */
		iters = 0;
		while (iters < 10) {
			/* Find the newly unfolded vertex ID */
			NN = XUGetNeighborEN(P1ID);
			for (i = 0; i < NN; i++) {
				neighVId = XUGetVertexNeighborVId(P1ID, i);
				if (neighVId == P2ID) {

					if (XUGetVertexNeighborVId(P1ID, (i + NN - 1) % NN) == UNotID) {
						UID = XUGetVertexNeighborVId(P1ID, (i + 1) % NN);
						break;
					}
					if (XUGetVertexNeighborVId(P1ID, (i + 1) % NN) == UNotID) {
						UID = XUGetVertexNeighborVId(P1ID, (i + NN - 1) % NN);
						break;
					}
				}
			}

			P1 = XUGetVertex(P1ID);
			P2 = XUGetVertex(P2ID);
			U = XUGetVertex(UID);

			vP1U = vectsubptr(P1, U);
			vP2U = vectsubptr(P2, U);

			P1U = (vP1U.x * vP1U.x + vP1U.y * vP1U.y + vP1U.z * vP1U.z);
			P2U = (vP2U.x * vP2U.x + vP2U.y * vP2U.y + vP2U.z * vP2U.z);

			cos12U = (P1P2 + P2U - P1U) / (2 * Math.sqrt(P1P2 * P2U));
			sin12U = (1 - cos12U * cos12U); /* Notice: Square of sine */

			/* Now compute three lengthes (squared) */
			UA = P2U + P2A - 2 * Math.sqrt(P2U * P2A)
					* (cos12A * cos12U - Math.sqrt(sin12A * sin12U));
			UB = P2U + P2B - 2 * Math.sqrt(P2U * P2B)
					* (cos12B * cos12U - Math.sqrt(sin12B * sin12U));
			UC = P2U + P2C - 2 * Math.sqrt(P2U * P2C)
					* (cos12C * cos12U - Math.sqrt(sin12C * sin12U));

			/* Now Judge Which Side to continue unfolding */
			if (UA > (UC + AC)) {/* Unfold along P1U */
				UNotID = P2ID;
				P2ID = UID;
				P1P2 = P1U;
				P2A = UA;
				P2B = UB;
				P2C = UC;
			} else if (UB > (UC + BC)) { /* Unfold along P2U */
				UNotID = P1ID;
				P1ID = UID;
				P1P2 = P2U;
				P1A = UA;
				P1B = UB;
				P1C = UC;
			} else { /* Stop Unfolding and compute T */
				/* Compute the actual lengthes */
				UC = Math.sqrt(UC);
				UA = Math.sqrt(UA);
				UB = Math.sqrt(UB);
				if (label[UID] == (byte) ALIVE) {
					if (la == ALIVE) {
						t1 = ComputeTInAcutev2(Ta, T[UID], UC, b, UA, Fc);
						if (t1 < Ta) { /* Reset A to NBAND and put into Heap */

							VertexIndexed vox = new VertexIndexed((float) Ta);
							vox.setPosition(vIDa);
							heap.add(vox);
							label[vIDa] = (byte) NBAND;
						}
					} else
						t1 = INFINITY;
					if (lb == ALIVE) {
						t2 = ComputeTInAcutev2(Tb, T[UID], UC, a, UB, Fc);
						if (t2 < Tb) { /* Reset B to NBAND and put into Heap */
							VertexIndexed vox = new VertexIndexed(Tb);
							vox.setPosition(vIDb);
							heap.add(vox);
							label[vIDb] = (byte) NBAND;
						}
					} else
						t2 = INFINITY;
					return Math.min(t1, t2);
				} else
					return (Math.min(Ta + b / Fc, Tb + a / Fc));
			}

			/* Update angles */
			cos12A = (P1P2 + P2A - P1A) / (2 * Math.sqrt(P1P2 * P2A));
			if (P2B != 0)
				cos12B = (P1P2 + P2B - P1B) / (2 * Math.sqrt(P1P2 * P2B));
			cos12C = (P1P2 + P2C - P1C) / (2 * Math.sqrt(P1P2 * P2C));

			sin12A = 1 - cos12A * cos12A;
			sin12B = 1 - cos12B * cos12B;
			sin12C = 1 - cos12C * cos12C;

			iters++;
		}/* End of while loop */

		return (Math.min(Ta + b / Fc, Tb + a / Fc));
	}

	protected double ComputeTInAcutev2(double Ta, double Tb, double a,
			double b, double c, double Fc) {
		double t1, t2, t, CD, costheta;
		double aa, bb, cc, u, tmp;

		costheta = (a * a + b * b - c * c) / (2 * a * b);

		u = Tb - Ta;

		Fc = 1 / Fc; /* Inverted here ! */

		aa = a * a + b * b - 2 * a * b * costheta;
		bb = 2 * b * u * (a * costheta - b);
		cc = b * b * (u * u - a * a * Fc * Fc * (1 - costheta * costheta));

		tmp = bb * bb - 4 * aa * cc;
		if (tmp < 0)
			return (Math.min(b * Fc + Ta, a * Fc + Tb));

		tmp = Math.sqrt(tmp);

		t1 = (-bb + tmp) / (2 * aa);
		t2 = (-bb - tmp) / (2 * aa);
		t = Math.max(t1, t2);
		CD = (b * (t - u)) / t;

		if ((u < t) && (a * costheta < CD) && (CD < (a / costheta)))
			return (t + Ta);
		else
			return (Math.min(b * Fc + Ta, a * Fc + Tb));

	}

}
