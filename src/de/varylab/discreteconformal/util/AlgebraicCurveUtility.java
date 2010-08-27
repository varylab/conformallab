package de.varylab.discreteconformal.util;

import static java.lang.Math.abs;
import static java.lang.Math.signum;
import geom3d.Point;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;


public class AlgebraicCurveUtility {

	public static Complex calculateCutModulus(CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		CoVertex v0 = cutInfo.cutRoot;
		CoVertex v1 = cutInfo.vertexCopyMap.get(v0);
		if (v1 == null) throw new RuntimeException("Connot calculate modulus. No cut-root copies found");
		CoVertex v2 = cutInfo.vertexCopyMap.get(v1);
		if (v2 == null) throw new RuntimeException("Connot calculate modulus. No second cut-root copy found");
		Point t0 = v0.getTextureCoord();
		Point t1 = v1.getTextureCoord();
		Point t2 = v2.getTextureCoord();
		Complex z0 = new Complex(t0.x() / t0.z(), t0.y() / t0.z());
		Complex z1 = new Complex(t1.x() / t1.z(), t1.y() / t1.z());
		Complex z2 = new Complex(t2.x() / t2.z(), t2.y() / t2.z());
		Complex w1 = z1.minus(z0);
		Complex w2 = z2.minus(z0);
		Complex tau = w2.divide(w1);
		
		int maxIter = 100;
		// move tau into its fundamental domain
		while ((abs(tau.re) > 0.5 || tau.im < 0 || tau.abs() < 1) && --maxIter > 0) {
			if (abs(tau.re) > 0.5) {
				tau.re -= signum(tau.re);
			}
			if (tau.im < 0) {
				tau = tau.times(-1);
			}
			if (tau.abs() < 1) {
				tau = tau.invert();
			}
		}
		return tau;
	}

	
}
