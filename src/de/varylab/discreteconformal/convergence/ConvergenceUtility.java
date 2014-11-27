package de.varylab.discreteconformal.convergence;

import static java.lang.Math.sqrt;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class ConvergenceUtility {

	public static double electrostaticEnergy(CoHDS hds) {
		double E = 0.0; 
		for (CoVertex v : hds.getVertices()) {
			double[] vPos = v.P;
			Pn.dehomogenize(vPos, vPos);
			for (CoVertex w : hds.getVertices()) {
				if (v == w) continue;
				double[] wPos = w.P;
				Pn.dehomogenize(wPos, wPos);
				double[] dir = Rn.subtract(null, vPos, wPos);
				double dsq = Rn.innerProduct(dir, dir);
				if (dsq == 0) continue;
				E += 1 / dsq;
			}
		}
		return E;
	}
	
	public static double[] getMaxMeanSumCrossRatio(CoHDS hds, double exp) {
		double rMax = 0.0;
		double rSum = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double q = getLengthCrossRatio(e);
			double qfun = (q + 1/q)/2 - 1;
			qfun = Math.pow(qfun, exp);
			rMax = Math.max(rMax, qfun);
			rSum += qfun;
		}
		return new double[] {rMax, rSum / hds.numFaces(), rSum};
	}
	
	public static double[] getMaxMeanSumMultiRatio(CoHDS hds, double exp) {
		double rMax = 0.0;
		double rSum = 0.0;
		for (CoFace f : hds.getFaces()) {
			double q = getLengthMultiRatio(f);
			double qfun = (q + 1/q)/2 - 1;
			qfun = Math.pow(qfun, exp);
			rMax = Math.max(rMax, qfun);
			rSum += qfun;
		}
		return new double[] {rMax, rSum / hds.numFaces(), rSum};
	}

	public static double[] getMaxMeanSumScaleInvariantCircumRadius(CoHDS hds) {
		double sqrtA = Math.sqrt(getTextureTriangleAreaSum(hds));
		double rMax = 0.0;
		double rSum = 0.0;
		for (CoFace f : hds.getFaces()) {
			double rad = getTextureCircumCircleRadius(f) / sqrtA;
			rMax = Math.max(rMax, rad);
			rSum += rad;
		}
		return new double[] {rMax, rSum / hds.numFaces(), rSum};
	}
	
	
	static double getTextureCircumCircleRadius(CoFace f) {
		CoEdge e = f.getBoundaryEdge();
		double a = e.getTexLength();
		double b = e.getNextEdge().getTexLength();
		double c = e.getPreviousEdge().getTexLength();
		double A = getTextureTriangleArea(f);
		return a*b*c / A / 4;
	}
	
	static double getTextureTriangleArea(CoFace f) {
		CoEdge e = f.getBoundaryEdge();
		double a = e.getTexLength();
		double b = e.getNextEdge().getTexLength();
		double c = e.getPreviousEdge().getTexLength();
		return sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) / 4;	
	}
	
	static double getTextureTriangleAreaSum(CoHDS hds) {
		double A = 0.0;
		for (CoFace f : hds.getFaces()) {
			A += getTextureTriangleArea(f);
		}
		return A;
	}
	
	
	static double getLengthCrossRatio(CoEdge e) {
		double a = e.getNextEdge().getLength();
		double b = e.getPreviousEdge().getLength();
		double c = e.getOppositeEdge().getNextEdge().getLength();
		double d = e.getOppositeEdge().getPreviousEdge().getLength();
		return (a * c) / (b * d);
	}
	
	static double getLengthMultiRatio(CoFace f) {
		double q = 1.0;
		for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
			q *= getLengthCrossRatio(e);
		}
		return q;
	}

}
