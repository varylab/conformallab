package de.varylab.discreteconformal.uniformization;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

import de.jreality.math.P2;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class SurfaceCurveUtility {

	private static Logger
		log = Logger.getLogger(SurfaceCurveUtility.class.getName());
	
	public static CoHDS createSurfaceCurves(
		FundamentalPolygon poly, 
		CoHDS surface,
		AdapterSet aSet,
		int maxElements,
		double maxDrawDistance,
		boolean includePoygon,
		boolean includeAxes
	) {
		List<double[][]> axesSegments = new ArrayList<double[][]>();
		List<double[][]> polySegments = new ArrayList<double[][]>();
		VisualizationUtility.getUniversalCoverSegments(poly, maxElements, maxDrawDistance, includePoygon, includeAxes, Color.BLACK, Color.BLACK, axesSegments, polySegments);
		
		List<double[][][]> allCurves = new ArrayList<double[][][]>();
		if (includeAxes) {
			for (double[][] ds : axesSegments) {
				List<double[][][]> I = intersectTriangulation(surface, ds);
				allCurves.addAll(I);
			}
		}
		if (includePoygon) {
			for (double[][] ds : polySegments) {
				List<double[][][]> I = intersectTriangulation(surface, ds);
				allCurves.addAll(I);
			}
		}
		
		CoHDS result = new CoHDS();
		for (double[][][] s : allCurves) {
			CoEdge e = result.addNewEdge();
			CoEdge eOpp = result.addNewEdge();
			e.linkOppositeEdge(eOpp);
			e.linkNextEdge(eOpp);
			eOpp.linkNextEdge(e);
			CoVertex v0 = result.addNewVertex();
			CoVertex v1 = result.addNewVertex();
			e.setTargetVertex(v0);
			eOpp.setTargetVertex(v1);
			v0.T[0] = s[0][0][0];
			v0.T[1] = s[0][0][1];
			v1.T[0] = s[0][1][0];
			v1.T[1] = s[0][1][1];
			v0.P = s[1][0];
			v1.P = s[1][1];
			Pn.normalize(v0.T, v0.T, Pn.HYPERBOLIC);
			Pn.normalize(v1.T, v1.T, Pn.HYPERBOLIC);
		}
		return result;
	}

	
	public static void createIntersectingEdges(FundamentalPolygon poly, CoHDS surface, CoHDS unwrapped, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo, AdapterSet aSet) {
		List<double[][]> axesSegments = new ArrayList<double[][]>();
		List<double[][]> polySegments = new ArrayList<double[][]>();
		VisualizationUtility.getUniversalCoverSegments(poly, 200, 10, true, false, Color.BLACK, Color.BLACK, axesSegments, polySegments);
		for (double[][] segment : polySegments) {
			double[] sLine = P2.lineFromPoints(null, segment[0], segment[1]);
			for (CoEdge e : unwrapped.getPositiveEdges()) {
				CoVertex s = e.getStartVertex();
				CoVertex t = e.getTargetVertex();
				double[] st = {s.T[0] / s.T[3], s.T[1] / s.T[3], 1};
				double[] tt = {t.T[0] / t.T[3], t.T[1] / t.T[3], 1};
				double[][] edgeSegment = {st, tt};
				double[][] edgePointSegment = {s.P, t.P};
				double[] eLine = P2.lineFromPoints(null, st, tt);
				double[] newPos = P2.pointFromLines(null, sLine, eLine);
				// split only if edge is really intersected
				if (isBetween(newPos, edgeSegment, 1E-6) && isBetween(newPos, segment, 5E-6)) {
					double[] newPoint = getPointOnSegment(newPos, edgeSegment, edgePointSegment);
					List<CoEdge> surfaceEdges = findCorrespondingSurfaceEdges(e, cutInfo, surface);
					if (surfaceEdges.isEmpty()) continue;
					CoEdge splitEdge = null;
					double[] splitPoint = null;
					// select intersected part
					for (CoEdge se : surfaceEdges) {
						CoVertex vs = se.getStartVertex();
						CoVertex vt = se.getTargetVertex();
						double[][] pSegment = {vs.P, vt.P};
						if (isBetween(newPoint, pSegment, 1E-20)) {
							splitEdge = se;
							splitPoint = newPoint;
						}
					}
					if (splitEdge == null) {
						log.warning("no corresponding intersected edge found on surface");
						continue;
					}
					double[] splitStart = splitEdge.getStartVertex().P;
					double[] splitTarget = splitEdge.getTargetVertex().P;
					// do not split if split point snaps to an edge end
					if (Pn.distanceBetween(splitPoint, splitStart, Pn.EUCLIDEAN) < 1.E-6) {
						continue;
					}
					if (Pn.distanceBetween(splitPoint, splitTarget, Pn.EUCLIDEAN) < 1.E-6) {
						continue;
					}
					CoVertex newVertex = TopologyAlgorithms.splitEdge(splitEdge);
					newVertex.P = splitPoint;
					newVertex.T = P2.imbedP2InP3(null, newPos);
				}
			}
		}
	}
	
	
	private static List<CoEdge> findCorrespondingSurfaceEdges(
		CoEdge sourceEdge, 
		CuttingInfo<CoVertex, CoEdge, CoFace> sourceCutInfo,
		CoHDS target
	) {
		CoFace flSource = sourceEdge.getLeftFace();
		CoFace frSource = sourceEdge.getRightFace();
		if (flSource == null) {
			flSource = sourceCutInfo.edgeCutMap.get(sourceEdge).getRightFace();
		}
		if (frSource == null) {
			frSource = sourceCutInfo.edgeCutMap.get(sourceEdge).getLeftFace();
		}
		assert flSource != null;
		assert frSource != null;
		CoFace flTarget = target.getFace(flSource.getIndex());
		CoFace frTarget = target.getFace(frSource.getIndex());
		List<CoEdge> result = HalfEdgeUtils.findEdgesBetweenFaces(flTarget, frTarget);
		return result;
	}
	
	

	private static List<double[][][]> intersectTriangulation(CoHDS T, double[][] segment) {
		List<double[][][]> result = new LinkedList<double[][][]>();
		for (CoFace f : T.getFaces()) {
			CoEdge be = f.getBoundaryEdge();
			CoVertex v0 = be.getStartVertex();
			CoVertex v1 = be.getNextEdge().getStartVertex();
			CoVertex v2 = be.getPreviousEdge().getStartVertex();
			double[] p0 = {v0.T[0] / v0.T[3], v0.T[1] / v0.T[3], 1};
			double[] p1 = {v1.T[0] / v1.T[3], v1.T[1] / v1.T[3], 1};
			double[] p2 = {v2.T[0] / v2.T[3], v2.T[1] / v2.T[3], 1};
			double[][] ts0 = {p0, p1};
			double[][] ts1 = {p1, p2};
			double[][] ts2 = {p2, p0};
			double[][] ps0 = {v0.P, v1.P};
			double[][] ps1 = {v1.P, v2.P};
			double[][] ps2 = {v2.P, v0.P};
			double[] line = P2.lineFromPoints(null, segment[0], segment[1]);
			double[][] tri = {p0, p1, p2};
			double[][] chopped = P2.chopConvexPolygonWithLine(tri, line);
			if (chopped == null || chopped == tri) {
				continue;
			} else {
				double[] c0 = null;
				double[] c1 = null;
				for (double[] p : chopped) {
					if (Arrays.equals(p, p0) || Arrays.equals(p, p1) || Arrays.equals(p, p2)) {
						continue;
					}
					if (c0 == null) {
						c0 = p;
					} else {
						c1 = p;
					}
				}
				if (c0 == null || c1 == null) {
					continue;
				}
				if (isOnSegment(c0, segment) && isOnSegment(c1, segment)) {
					double[] sp0 = v0.P;
					double[] sp1 = v1.P;
					if (isOnSegment(c0, ts0)) {
						sp0 = getPointOnSegment(c0, ts0, ps0);
					} else if (isOnSegment(c0, ts1)) {
						sp0 = getPointOnSegment(c0, ts1, ps1);
					} else if (isOnSegment(c0, ts2)) {
						sp0 = getPointOnSegment(c0, ts2, ps2);
					}
					if (isOnSegment(c1, ts0)) {
						sp1 = getPointOnSegment(c1, ts0, ps0);
					} else if (isOnSegment(c1, ts1)) {
						sp1 = getPointOnSegment(c1, ts1, ps1);
					} else if (isOnSegment(c1, ts2)) {
						sp1 = getPointOnSegment(c1, ts2, ps2);
					}
					double[][][] curveSegment = {{c0, c1}, {sp0, sp1}};
					result.add(curveSegment);
				}
			}
		}
		return result;
	}
	
	static boolean isBetween(double[] p, double[][] s, double snapThresh) {
		Pn.dehomogenize(p, p);
		Pn.dehomogenize(s, s);
		double[] ps0 = Rn.subtract(null, s[0], p);
		double[] ps1 = Rn.subtract(null, s[1], p);
		double d1 = Pn.norm(ps0, Pn.EUCLIDEAN);
		double d2 = Pn.norm(ps0, Pn.EUCLIDEAN);
		if (d1 < snapThresh || d2 < snapThresh) return false;
		double[] s0s1 = Rn.subtract(null, s[0], s[1]);
		double[] cross = Rn.crossProduct(null, ps0, ps1);
		double dot = Rn.innerProduct(ps0, ps1);
		if (Rn.euclideanNorm(cross) > 1E-7) return false;
		if (dot > 0) return false;
		if (dot > Rn.euclideanNormSquared(s0s1)) return false;
	    return true;
	}
	
	
	static boolean isOnSegment(double[] p, double[][] s) {
		double l = Pn.distanceBetween(s[0], s[1], Pn.EUCLIDEAN);
		double l1 = Pn.distanceBetween(s[0], p, Pn.EUCLIDEAN);
		double l2 = Pn.distanceBetween(s[1], p, Pn.EUCLIDEAN);
		return Math.abs(l1 + l2 - l) < 1E-7;
	}
	
	
	private static double[] getPointOnSegment(double[] p, double[][] source, double[][] target) {
//		Pn.dehomogenize(target, target);
		double l = Pn.distanceBetween(source[0], source[1], Pn.HYPERBOLIC);
		double l1 = Pn.distanceBetween(source[0], p, Pn.HYPERBOLIC) / l;
		double l2 = Pn.distanceBetween(source[1], p, Pn.HYPERBOLIC) / l;
		if (l1 < 1E-8 && l2 < 1E-8) {
			log.warning("point is close to both segment ends");
		}
		return Rn.linearCombination(null, l1, target[1], l2, target[0]);
	}
	
	
	
}
