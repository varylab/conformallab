package de.varylab.discreteconformal.functional;

import static java.lang.Math.PI;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.math.ComplexUtility;

public class HyperIdealHyperellipticUtility {

	public static void calculateCircleIntersections(CoHDS hds) {
		for (CoEdge e : hds.getPositiveEdges()) {
			if (HalfEdgeUtils.isBoundaryEdge(e)) continue;
			double[] a = Pn.dehomogenize(null, e.getStartVertex().P);
			double[] b = Pn.dehomogenize(null, e.getTargetVertex().P);
			double[] c = Pn.dehomogenize(null, e.getNextEdge().getTargetVertex().P);
			double[] d = Pn.dehomogenize(null, e.getOppositeEdge().getNextEdge().getTargetVertex().P);
			double dist1 = Pn.distanceBetween(a, b, Pn.EUCLIDEAN);
			// rotate such that no vertex is north
			double[] center = {0,0,0};
			Rn.add(center, center, a);
			Rn.add(center, center, b);
			Rn.add(center, center, c);
			Rn.normalize(center, center);
			MatrixBuilder m = MatrixBuilder.euclidean().rotateFromTo(center, new double[]{0,0,1});
			m.getMatrix().transformVector(a);
			m.getMatrix().transformVector(b);
			m.getMatrix().transformVector(c);
			m.getMatrix().transformVector(d);
			
			// move all vertices south
//			m = MatrixBuilder.hyperbolic().translateFromTo(new double[]{0,0,0,1}, new double[]{0,0,-0.9,1.0});
//			m.getMatrix().transformVector(a);
//			m.getMatrix().transformVector(b);
//			m.getMatrix().transformVector(c);
//			m.getMatrix().transformVector(d);
//			Pn.dehomogenize(a, a);
//			Pn.dehomogenize(b, b);
//			Pn.dehomogenize(c, c);
//			Pn.dehomogenize(d, d);
			
			double thresh = 1.0 - 1e-3;
			if (a[2] > thresh || b[2] > thresh || c[2] > thresh || d[2] > thresh) {
				assert false : "there sould be no vertex at the north pole";
			}
			
			// calculate angle
			Complex ac = ComplexUtility.stereographic(a);
			Complex bc = ComplexUtility.stereographic(b);
			Complex cc = ComplexUtility.stereographic(c);
			Complex dc = ComplexUtility.stereographic(d);
			double alpha = Math.abs(cc.minus(bc).divide(cc.minus(ac)).arg());
			double beta = Math.abs(dc.minus(bc).divide(dc.minus(ac)).arg());
			double theta = PI - alpha - beta;
			e.setTheta(theta);
			e.getOppositeEdge().setTheta(theta);
		}
	}
	
}
