package de.varylab.discreteconformal.holomorphicformsexperiments;

import java.util.Set;

import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.util.AngleUtilities;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.math.ComplexUtility;

public class Utility {

	public static double calculateLargestEdgeLength(
		CoHDS S,
		AdapterSet a,
		Set<CoVertex> branch,
		boolean branchClustering
	) {
		double maxLength = 0.0;
		double NEIGHBORHOOD_THRESHOLD = 0.05;
		for (CoEdge e : S.getPositiveEdges()) {
			double length = a.get(Length.class, e, Double.class);
			if (length > maxLength) { maxLength = length; }
		}
		if (branchClustering) {
			for (CoVertex b : branch) {
				double[] branchPos = a.getD(Position4d.class, b).clone();
				Matrix T = MatrixBuilder.euclidean().rotateFromTo(branchPos, new double[] {0, 0, 1}).getMatrix();
				for (CoEdge e : S.getPositiveEdges()) {
					double[] es = a.getD(Position4d.class, e.getStartVertex()).clone();
					double[] et = a.getD(Position4d.class, e.getTargetVertex()).clone();
					// rotate branch point up north
					T.transformVector(es);
					T.transformVector(et);
					// project down
					Complex esc = ComplexUtility.stereographic(es);
					Complex etc = ComplexUtility.stereographic(et);
					// only consider a neighborhood of the branch point
					if (esc.abs() > NEIGHBORHOOD_THRESHOLD || etc.abs() > NEIGHBORHOOD_THRESHOLD) { continue; }
					esc = esc.sqr();
					etc = etc.sqr();
					// project back
					double[] esci = ComplexUtility.inverseStereographic(esc);
					double[] etci = ComplexUtility.inverseStereographic(etc);
					// measure distance
					double d = Pn.distanceBetween(esci, etci, Pn.EUCLIDEAN);
					if (d > maxLength) { maxLength = d; }
				}
			}
		}
		return maxLength;
	}
	
	public static double calculateSmallestAngle(
		CoHDS S,
		AdapterSet a
	) {
		double minAngle = Double.MAX_VALUE;
		for (CoEdge e : S.getEdges()) {
			CoEdge e1 = e.getNextEdge();
			CoEdge e2 = e.getPreviousEdge();
			double angle = AngleUtilities.angle(e1, e2, a);
			if (angle < minAngle) { minAngle = angle; }
		}
		return minAngle;
	}
	
}
