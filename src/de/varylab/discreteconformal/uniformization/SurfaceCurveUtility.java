package de.varylab.discreteconformal.uniformization;

import static de.varylab.discreteconformal.adapter.HyperbolicModel.Klein;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import de.jreality.geometry.IndexedLineSetFactory;
import de.jreality.math.P2;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.scene.IndexedLineSet;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.NodeIndexComparator;

public class SurfaceCurveUtility {

	private static Logger
		log = Logger.getLogger(SurfaceCurveUtility.class.getName());
	
	public static IndexedLineSet createSurfaceCurves(
		FundamentalPolygon poly, 
		CoHDS surface,
		AdapterSet aSet,
		int maxElements,
		double maxDrawDistance,
		boolean includePoygon,
		boolean includeAxes,
		int signature
	) {
		List<double[][]> axesSegments = new ArrayList<double[][]>();
		List<double[][]> polySegments = new ArrayList<double[][]>();
		VisualizationUtility.createUniversalCover(
			poly, 
			Klein,
			maxElements, maxDrawDistance, 
			includePoygon, includeAxes, 
			axesSegments, polySegments, 
			null, null, null
		);
		
		List<double[][][]> allCurves = new ArrayList<double[][][]>();
		if (includeAxes) {
			for (double[][] ds : axesSegments) {
				List<double[][][]> I = intersectTriangulation(surface, ds, signature);
				allCurves.addAll(I);
			}
		}
		if (includePoygon) {
			for (double[][] ds : polySegments) {
				List<double[][][]> I = intersectTriangulation(surface, ds, signature);
				allCurves.addAll(I);
			}
		}
		
//		CoHDS result = new CoHDS();
		double[][] vData = new double[allCurves.size() * 2][];
		double[][] tData = new double[allCurves.size() * 2][2];
		int[][] eData = new int[allCurves.size()][2];
		int index = 0;
		for (double[][][] s : allCurves) {
//			CoEdge e = result.addNewEdge();
//			CoEdge eOpp = result.addNewEdge();
//			e.linkOppositeEdge(eOpp);
//			e.linkNextEdge(eOpp);
//			eOpp.linkNextEdge(e);
//			CoVertex v0 = result.addNewVertex();
//			CoVertex v1 = result.addNewVertex();
//			e.setTargetVertex(v0);
//			eOpp.setTargetVertex(v1);
//			v0.T[0] = s[0][0][0];
//			v0.T[1] = s[0][0][1];
//			v1.T[0] = s[0][1][0];
//			v1.T[1] = s[0][1][1];
//			v0.P = s[1][0];
//			v1.P = s[1][1];
//			Pn.normalize(v0.T, v0.T, Pn.HYPERBOLIC);
//			Pn.normalize(v1.T, v1.T, Pn.HYPERBOLIC);
			vData[index + 0] = s[1][0];
			vData[index + 1] = s[1][1];
			tData[index + 0][0] = s[0][0][0];
			tData[index + 0][1] = s[0][0][1];
			tData[index + 1][0] = s[0][1][0];
			tData[index + 1][1] = s[0][1][1];
			eData[index/2][0] = index;
			eData[index/2][1] = index + 1;
			index += 2;
		}
		IndexedLineSetFactory ilsf = new IndexedLineSetFactory();
		ilsf.setVertexCount(vData.length);
		ilsf.setEdgeCount(eData.length);
		ilsf.setVertexCoordinates(vData);
		ilsf.setVertexTextureCoordinates(tData);
		ilsf.setEdgeIndices(eData);
		ilsf.update();
		return ilsf.getIndexedLineSet();
	}

	
	public static Set<Set<CoVertex>> createIntersectionVertices(
		FundamentalPolygon poly, 
		CoHDS surface, 
		CoHDS domain, 
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo, AdapterSet aSet,
		double snapTolerance,
		int signature
	) {
		Set<Set<CoVertex>> result = new HashSet<Set<CoVertex>>();
		List<double[][]> axesSegments = new ArrayList<double[][]>();
		List<double[][]> polySegments = new ArrayList<double[][]>();
		VisualizationUtility.createUniversalCover(
			poly, 
			Klein, 
			200, 10, 
			true, false, 
			axesSegments, polySegments, 
			null, null, null
		);
		for (double[][] segment : polySegments) {
			Set<CoVertex> segmentVertices = new TreeSet<>(new NodeIndexComparator<CoVertex>());
			createIntersectionVertices(segment, surface, domain, cutInfo, snapTolerance, signature, segmentVertices);
			result.add(segmentVertices);
		}
		return result;
	}


	public static void createIntersectionVertices(
		double[][] segment,
		boolean segmentOnly,
		CoHDS surface, 
		CoHDS domain,
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo,
		double snapTolerance, 
		int signature, 
		Set<CoVertex> result
	) {
		double[] polygonLine = P2.lineFromPoints(null, segment[0], segment[1]);
		Rn.normalize(polygonLine, polygonLine);
		for (CoEdge domainEdge : domain.getPositiveEdges()) {
			CoVertex s = domainEdge.getStartVertex();
			CoVertex t = domainEdge.getTargetVertex();
			double[] st = {s.T[0] / s.T[3], s.T[1] / s.T[3], 1};
			double[] tt = {t.T[0] / t.T[3], t.T[1] / t.T[3], 1};
			double[][] domainSegment = {st, tt};
			double[][] surfaceSegment = {s.P, t.P};
			double[] domainEdgeLine = P2.lineFromPoints(null, st, tt);
			Rn.normalize(domainEdgeLine, domainEdgeLine);
		 	double edgeOnLine = Rn.euclideanNormSquared(Rn.crossProduct(null, polygonLine, domainEdgeLine));
		 	if (edgeOnLine < 1E-8) {
		 		List<CoEdge> surfaceEdges = findCorrespondingSurfaceEdges(domainEdge, cutInfo, surface);
		 		for (CoEdge e : surfaceEdges) {
		 			result.add(e.getStartVertex());
		 			result.add(e.getTargetVertex());
		 		}
		 		continue;
		 	}
			double[] newDomainPoint = P2.pointFromLines(null, polygonLine, domainEdgeLine);
			// split only if edge is really intersected
			if (isOnSegment(newDomainPoint, domainSegment) && (!segmentOnly || isOnSegment(newDomainPoint, segment))) {
				double[] newPoint = getPointOnCorrespondingSegment(newDomainPoint, domainSegment, surfaceSegment, signature);
				List<CoEdge> surfaceEdges = findCorrespondingSurfaceEdges(domainEdge, cutInfo, surface);
				if (surfaceEdges.isEmpty()) continue;
				CoEdge splitEdge = null;
				// select intersected part
				boolean snap = false;
				for (CoEdge se : surfaceEdges) {
					CoVertex vs = se.getStartVertex();
					CoVertex vt = se.getTargetVertex();
					if (result.contains(vs) || result.contains(vt)) {
						log.info("snapped to cycle vertex");
						snap = true;
						break;
					}
					double[][] pSegment = {vs.P, vt.P};
					if (isOnSegment(newPoint, pSegment)) {
						splitEdge = se;
						double d1 = getDistanceBetween(newPoint, vs.P, Pn.EUCLIDEAN);
						double d2 = getDistanceBetween(newPoint, vt.P, Pn.EUCLIDEAN);
						if (d1 < snapTolerance) {
							result.add(vs);
							snap = true;
						}
						if (d2 < snapTolerance) {
							result.add(vt);
							snap = true;
						}
						break;
					}
				}
				if (snap) continue;
				if (splitEdge == null) {
					log.warning("no corresponding intersected edge found on surface");
					continue;
				}
				CoVertex newVertex = TopologyAlgorithms.splitEdge(splitEdge);
				result.add(newVertex);
				newVertex.P = newPoint;
				newVertex.T = P2.imbedP2InP3(null, newDomainPoint);
			}
		}
	}
	
	public static void createIntersectionVertices(
			double[][] segment,
			CoHDS surface, 
			CoHDS domain,
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo,
			double snapTolerance, 
			int signature, 
			Set<CoVertex> result
		) {
		createIntersectionVertices(segment, true, surface, domain, cutInfo, snapTolerance, signature, result);
	}
	
	private static List<CoEdge> findCorrespondingSurfaceEdges(
		CoEdge sourceEdge, 
		CuttingInfo<CoVertex, CoEdge, CoFace> sourceCutInfo,
		CoHDS target
	) {
		CoFace flSource = sourceEdge.getLeftFace();
		CoFace frSource = sourceEdge.getRightFace();
		if (flSource == null) {
			CoEdge idEdge = sourceCutInfo.edgeCutMap.get(sourceEdge);
			if (idEdge != null) {
				flSource = idEdge.getRightFace();
			}
		}
		if (frSource == null) {
			CoEdge idEdge = sourceCutInfo.edgeCutMap.get(sourceEdge);
			if (idEdge != null) {
				frSource = idEdge.getLeftFace();
			}
		}
		assert flSource != null || frSource != null;
		CoFace flTarget = null;
		CoFace frTarget = null;
		if (flSource != null) {
			flTarget = target.getFace(flSource.getIndex());
		}
		if (frSource != null) {
			frTarget = target.getFace(frSource.getIndex());
		}
		if (flTarget == null) {
			flTarget = frTarget;
			frTarget = null;
		}
		List<CoEdge> result = HalfEdgeUtils.findEdgesBetweenFaces(flTarget, frTarget);
		return result;
	}
	
	

	private static List<double[][][]> intersectTriangulation(CoHDS T, double[][] segment, int signature) {
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
						sp0 = getPointOnCorrespondingSegment(c0, ts0, ps0, signature);
					} else if (isOnSegment(c0, ts1)) {
						sp0 = getPointOnCorrespondingSegment(c0, ts1, ps1, signature);
					} else if (isOnSegment(c0, ts2)) {
						sp0 = getPointOnCorrespondingSegment(c0, ts2, ps2, signature);
					}
					if (isOnSegment(c1, ts0)) {
						sp1 = getPointOnCorrespondingSegment(c1, ts0, ps0, signature);
					} else if (isOnSegment(c1, ts1)) {
						sp1 = getPointOnCorrespondingSegment(c1, ts1, ps1, signature);
					} else if (isOnSegment(c1, ts2)) {
						sp1 = getPointOnCorrespondingSegment(c1, ts2, ps2, signature);
					}
					double[][][] curveSegment = {{c0, c1}, {sp0, sp1}};
					result.add(curveSegment);
				}
			}
		}
		return result;
	}
	
	static boolean isOnSegment(double[] p, double[][] s) {
		Pn.dehomogenize(p, p);
		Pn.dehomogenize(s, s);
		double[] ps0 = Rn.subtract(null, s[0], p);
		double[] ps1 = Rn.subtract(null, s[1], p);
		double d1 = Pn.norm(ps0, Pn.EUCLIDEAN);
		double d2 = Pn.norm(ps0, Pn.EUCLIDEAN);
		if (d1 < 0.0 || d2 < 0.0) return false;
		double[] s0s1 = Rn.subtract(null, s[0], s[1]);
		double[] cross = Rn.crossProduct(null, ps0, ps1);
		double dot = Rn.innerProduct(ps0, ps1);
		if (Rn.euclideanNorm(cross) > 1E-7) return false;
		if (dot > 0) return false;
		if (dot > Rn.euclideanNormSquared(s0s1)) return false;
	    return true;
	}
	
	
//	static boolean isOnSegment(double[] p, double[][] s) {
//		double l = Pn.distanceBetween(s[0], s[1], Pn.EUCLIDEAN);
//		double l1 = Pn.distanceBetween(s[0], p, Pn.EUCLIDEAN);
//		double l2 = Pn.distanceBetween(s[1], p, Pn.EUCLIDEAN);
//		return Math.abs(l1 + l2 - l) < 1E-7;
//	}
	
	
	static double[] getPointOnCorrespondingSegment(double[] p, double[][] source, double[][] target, int signature) {
		double l = Pn.distanceBetween(source[0], source[1], signature);
		double l1 = Pn.distanceBetween(source[0], p, signature) / l;
		double l2 = Pn.distanceBetween(source[1], p, signature) / l;
		if (Double.isNaN(l1)) {
			return target[0];
		}
		if (Double.isNaN(l2)) {
			return target[1];
		}
		double[] t0d = Pn.dehomogenize(null, target[0]);
		double[] t1d = Pn.dehomogenize(null, target[1]);
		return Rn.linearCombination(null, l1, t1d, l2, t0d);
	}
	
	
	static double getDistanceBetween(double[] p1, double[] p2, int signature) {
		double d = Pn.distanceBetween(p1, p2, signature);
		if (Double.isNaN(d)) {
			return 0.0;
		} else {
			return d;
		}
	}
	
}
