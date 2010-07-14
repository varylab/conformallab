package de.varylab.discreteconformal.util;

import static java.lang.Math.abs;
import static java.lang.Math.signum;
import geom3d.Point;

import java.util.Iterator;
import java.util.Set;

import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;


public class AlgebraicCurveUtility {

	public static Complex calculateCutModulus(CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		Iterator<Set<CoEdge>> pathIt = cutInfo.paths.iterator();
		Set<CoEdge> path1 = pathIt.next();
		Set<CoEdge> path2 = pathIt.next();
		
		Set<CoVertex> path1V = PathUtility.getVerticesOnPath(path1);
		Set<CoVertex> path2V = PathUtility.getVerticesOnPath(path2);
	
		Complex tau = null;
		double maxIm = 0;
		for (CoVertex v : path1V) {
			for (CoVertex vv : path2V) {
				CoVertex vc = cutInfo.vertexCopyMap.get(v);
				CoVertex vvc = cutInfo.vertexCopyMap.get(vv);
				if (vc == null || vvc == null) {
					continue;
				}
				Point tv = v.getTextureCoord();
				Point tvv = vv.getTextureCoord();
				Point tvc = vc.getTextureCoord();
				Point tvvc = vvc.getTextureCoord();
				Complex z0 = new Complex(tv.x() / tv.z(), tv.y() / tv.z());
				Complex z1 = new Complex(tvv.x() / tvv.z(), tvv.y() / tvv.z());
				Complex z2 = new Complex(tvc.x() / tvc.z(), tvc.y() / tvc.z());
				Complex z3 = new Complex(tvvc.x() / tvvc.z(), tvvc.y() / tvvc.z());
				Complex w1 = z0.minus(z2);
				Complex w2 = z1.minus(z3);
				Complex tauTmp = w2.divide(w1);
				double absIm = Math.abs(tauTmp.im);
				if (absIm > maxIm) {
					tau = tauTmp;
					maxIm = absIm;
				}
			}
		}
		assert tau != null;
		int maxIter = 100;
		while ((abs(tau.re) > 0.5 || tau.abs() < 1) && --maxIter > 0) {
			if (abs(tau.re) > 0.5) {
				tau.re -= signum(tau.re);
			}
			if (tau.abs() < 1) {
				tau = tau.invert();
				tau = tau.times(-1);
			}
		}
		return tau;
	}

	
}
