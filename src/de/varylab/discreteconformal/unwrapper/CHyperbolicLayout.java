package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;
import geom3d.Point;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.util.UniformizationUtility;

public class CHyperbolicLayout {

	
	public static Map<CoEdge, Double> getLengthMap(CoHDS hds, Vector u) {
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = getNewLength(e, u);
			lMap.put(e, l);
			lMap.put(e.getOppositeEdge(), l);
		}
		return lMap;
	}
	
	
	/** 
	 * Do flat layout for a HDS and a metric vector u
	 * @param hds mesh
	 * @param u new metric
	 * @param angleMapParam may be null
	 */
	public static Set<CoVertex> doLayout(CoHDS hds, Vector u) {
		Set<CoVertex> cutSet = new HashSet<CoVertex>();
		Map<CoEdge, Double> lMap = getLengthMap(hds, u);
		
		Set<CoVertex> visited = new HashSet<CoVertex>(hds.numVertices());
		Queue<CoVertex> Qv = new LinkedList<CoVertex>();
		Queue<CoEdge> Qe = new LinkedList<CoEdge>();
		// start
		CoEdge e0 = hds.getEdge(hds.numEdges() / 2);
		CoEdge e1 = e0.getOppositeEdge();
		CoVertex v1 = e0.getStartVertex();
		CoVertex v2 = e0.getTargetVertex();
		// queued data
		Qv.offer(v1);
		Qv.offer(v2);
		Qe.offer(e1);
		Qe.offer(e0);

		// vertices
		Double d = lMap.get(e0);
		
		Point p0 = new Point(0, 0, 1);
		v1.setTextureCoord(p0);
		Point p1 = normalize(new Point(sinh(d), 0, cosh(d)).asPoint());
		v2.setTextureCoord(p1);
		
		visited.add(v1);
		visited.add(v2);
		
		while (!Qv.isEmpty() && hds.numVertices() > visited.size()) {
			CoVertex v = Qv.poll();
			CoEdge inE = Qe.poll();
			CoEdge outE = inE.getOppositeEdge();
			
			CoEdge e = inE.getNextEdge();
			while (e != outE) {
				CoEdge next = e.getNextEdge();
				CoEdge prev = e.getPreviousEdge();
				CoVertex aVertex = next.getTargetVertex();
				CoVertex bVertex = prev.getTargetVertex();
				CoVertex cVertex = e.getTargetVertex();

				Double alpha = next.getAlpha();
				if (e.getLeftFace() == null) { // a boundary edge
					alpha = 2*PI - getAngleSum(v);
					e = e.getOppositeEdge().getNextEdge();
					continue;
				}
				if (!visited.contains(cVertex)) {
					d = lMap.get(e);
					double dCheck = lMap.get(next);
					Point A = aVertex.getTextureCoord();
					Point B = bVertex.getTextureCoord();
					
					Point C = layoutTriangle(A, B, alpha, d, dCheck);
					if (C != null) {
						cVertex.setTextureCoord(C);
						visited.add(cVertex);
						Qv.offer(cVertex);
						Qe.offer(e);
					} else {
						System.out.println("Vertex layout skipped for " + v + " checking triangle...");
						boolean check = UniformizationUtility.checkTriangle(e.getLeftFace(), lMap);
						System.out.println("triangle is valid: " + check);
						break;
					}
				} else {
					cutSet.add(cVertex);
				}
				e = e.getOppositeEdge().getNextEdge();
			}
		}
		
		System.out.println("Layout: Visited points: " + visited.size() + "/" + hds.numVertices());
		return cutSet;
	}
	
	
	
//	public static Set<CoEdge> makeCutGraph(Set<CoVertex> cutSet) {
//		Set<CoEdge> cutGraph = new HashSet<CoEdge>();
//		for (CoVertex v : cutSet) {
//			CoEdge e1 = null;
//			CoEdge e2 = null;
//			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
//				if (cutSet.contains(e.getStartVertex())) {
//					if (e1 == null) {
//						e1 = e;
//					} else {
//						e2 = e;
//						break;
//					}
//				}
//			}
//			if (e1 == null || e2 == null) {
//				throw new IllegalArgumentException("Invalid CutSet in makeCutGraph");
//			}
//			cutGraph.add(e1);
//			cutGraph.add(e1.getOppositeEdge());
//			cutGraph.add(e2);
//			cutGraph.add(e2.getOppositeEdge());
//		}
//		return cutGraph;
//	}
	
	
	
	private static Point layoutTriangle(Point A, Point B, double alpha, double d, double dP) {
		Point BHat = new Point(B.x(), B.y(), -B.z());
		Point AHat = new Point(A.x(), A.y(), -A.z());
		Point lAB = normalize(new Point(A).cross(B).asPoint());
		Point At = normalize(new Point(lAB).cross(BHat).asPoint());
		Point AtPerp = normalize(new Point(AHat).cross(BHat).asPoint());
		Point Ct = normalize(new Point(At).times(cos(alpha)).add(new Point(AtPerp).times(sin(alpha))).asPoint());
		Point C1 = normalize(new Point(B).times(Math.cosh(d)).add(new Point(Ct).times(Math.sinh(d))).asPoint());
		Point C2 = normalize(new Point(B).times(Math.cosh(d)).subtract(new Point(Ct).times(Math.sinh(d))).asPoint());
		double d1 = Double.MAX_VALUE;
		double d2 = Double.MAX_VALUE;
		try {
			d1 = Pn.distanceBetween(C1.get(), A.get(), Pn.HYPERBOLIC);
		} catch (IllegalArgumentException iae) {}
		try {
			d2 = Pn.distanceBetween(C2.get(), A.get(), Pn.HYPERBOLIC);
		} catch (IllegalArgumentException iae) {}
		double dif1 = Math.abs(d1 - dP);
		double dif2 = Math.abs(d2 - dP);
		Point C = dif1 < dif2 ? C1 : C2; 
		double dif = dif1 < dif2 ? dif1 : dif2;
		if (dif < 1E-5) {
			return C;
		} else {
			return null;
		}
	}
	
	
	
	
	private static Point normalize(Point p) {
		Pn.normalize(p.get(), p.get(), Pn.HYPERBOLIC);
		return p;
	}
	
	
	/**
	 * Calculate the angle sum at this vertex. Usually this will be 2PI, but at the boundary
	 * we sum only the inner angles
	 * @param v
	 * @return the angle sum
	 */
	public static Double getAngleSum(CoVertex v) {
		Double r = 0.0;
		List<CoEdge> star = incomingEdges(v);
		for (CoEdge e : star) {
			if (e.getLeftFace() != null) {
				r += e.getPreviousEdge().getAlpha();
			}
		}
		return r;
	}
	
	
	/**
	 * Calculate the edge length for the flat metric
	 * @param e
	 * @param u
	 * @return the new edge length
	 */
	public static Double getNewLength(CoEdge e, Vector u) {
		CoVertex v1 = e.getStartVertex();
		CoVertex v2 = e.getTargetVertex();
		Double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
		Double u2 = v2.getSolverIndex() >= 0 ? u.get(v2.getSolverIndex()) : 0.0;
		Double lambda = e.getLambda();
		Double lambdaNew = lambda + u1 + u2;
		return 2 * arsinh( exp(lambdaNew / 2) );
	}
	
	
	private static double arsinh(double x) {
		double r = x + sqrt(x*x + 1);
		return log(r);
	}
	
	
}
