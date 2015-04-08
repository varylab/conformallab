package de.varylab.discreteconformal.functional;

import static java.lang.Math.PI;
import de.jreality.math.Pn;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.math.ComplexUtility;

public class HyperIdealHyperellipticUtility {

	public static void calculateCircleIntersections(CoHDS hds) {
		for (CoEdge e : hds.getPositiveEdges()) {
			double[] a = Pn.dehomogenize(null, e.getStartVertex().P);
			double[] b = Pn.dehomogenize(null, e.getTargetVertex().P);
			double[] c = Pn.dehomogenize(null, e.getNextEdge().getTargetVertex().P);
			double[] d = Pn.dehomogenize(null, e.getOppositeEdge().getNextEdge().getTargetVertex().P);
			double thresh = 1.0 - 1e-3;
			if (a[2] > thresh || b[2] > thresh || c[2] > thresh || d[2] > thresh) {
				a[2] *= -1; b[2] *= -1;	c[2] *= -1;	d[2] *= -1; // interchange north and south
			}
			Complex ac = ComplexUtility.stereographic(a);
			Complex bc = ComplexUtility.stereographic(b);
			Complex cc = ComplexUtility.stereographic(c);
			Complex dc = ComplexUtility.stereographic(d);
			double alpha = Math.abs(cc.minus(bc).divide(cc.minus(ac)).arg());
			double beta = Math.abs(dc.minus(bc).divide(dc.minus(ac)).arg());
			e.setAlpha(PI - alpha - beta);
			e.getOppositeEdge().setAlpha(e.getAlpha());
		}
	}
	
}
