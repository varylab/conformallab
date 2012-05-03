package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.util.HashMap;
import java.util.Map;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.HalfEdgeDataStructure;

public class IsothermicUtility {

	/**
	 * Calculate the angle between the edges that belong to alpha1 and alpha2.
	 * @param alpha1
	 * @param alpha2
	 * @param alpha3
	 * @return
	 */
	public static double calculateTriangleAngle(double alpha1, double alpha2, double alpha3) {
		alpha1 = IsothermicUtility.normalizeAngle(alpha1);
		alpha2 = IsothermicUtility.normalizeAngle(alpha2);
		alpha3 = IsothermicUtility.normalizeAngle(alpha3);
		double beta = abs(alpha2 - alpha1);
		if ((alpha3 > alpha2 && alpha3 > alpha1) || (alpha3 < alpha2 && alpha3 < alpha1)) {
			return beta;
		} else {
			return PI - beta;
		}
	}

	public static double normalizeAngle(double a) {
		a %= 2*PI;
		if (a > PI/2) {
			return a - PI;
		} else if (a < -PI/2) {
			return PI + a;
		} else {
			return a;
		}
	}

	/**
	 * Returns the angle between v1 and v2 in the range ]-pi/2, pi/2]. 
	 * Where the sign is the sign of the determinant |N v1 v2|. 
	 * @param v1
	 * @param v2
	 * @param N
	 * @return
	 */
	public static double getSignedAngle(double[] N, double[] v1, double[] v2) {
		double[][] T = {N, v1, v2};
		double sign = Math.signum(Rn.determinant(T));
		double alpha = Rn.euclideanAngle(v1, v2);
		if (alpha > PI/2) {
			alpha = -(PI - alpha);
		}
		return sign * alpha;
	}

	
	public static Map<Integer, Integer> createUndirectedEdgeMap(HalfEdgeDataStructure<?, ?, ?> hds) {
		Map<Integer, Integer> result = new HashMap<Integer, Integer>();
		Integer i = 0;
		for (Edge<?,?,?> e : hds.getPositiveEdges()) {
			result.put(e.getIndex(), i);
			result.put(e.getOppositeEdge().getIndex(), i);
			i++;
		}
		return result;
	}
	
}
