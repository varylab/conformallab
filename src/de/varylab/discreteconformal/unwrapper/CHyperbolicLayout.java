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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.util.CuttingUtility;
import de.varylab.discreteconformal.heds.util.HomologyUtility;

public class CHyperbolicLayout {

	
	public static class HyperbolicLayoutContext {
		
		public Map<CoEdge, CoEdge>
			cutMap = new HashMap<CoEdge, CoEdge>();
		public List<Set<CoEdge>>
			paths = new ArrayList<Set<CoEdge>>();
		public Map<Set<CoEdge>, Set<CoEdge>>
			pathCutMap = new HashMap<Set<CoEdge>, Set<CoEdge>>();
		
	}
	
	
	
	/**
	 * Do flat layout for a HDS and a metric vector u
	 * @param hds mesh
	 * @param u new metric
	 * @param angleMapParam may be null
	 */
	public static HyperbolicLayoutContext doLayout(CoHDS hds, Vector u) {
		
		HyperbolicLayoutContext context = new HyperbolicLayoutContext();
		int X = hds.numVertices() - hds.numEdges() / 2 + hds.numFaces();
		int g = (2 - X) / 2;
		System.err.println("genus of the surface is " + g);
		if (g >= 2) {
			context.paths = HomologyUtility.getGeneratorPaths(hds.getVertex(0));
			Set<CoEdge> masterPath = new HashSet<CoEdge>();
			for (Set<CoEdge> path : context.paths) {
				masterPath.addAll(path);
			}
			for (CoEdge e : masterPath) {
				if (HalfEdgeUtils.isInteriorEdge(e)) {
					context.cutMap.put(e, e.getOppositeEdge());
					Map<CoVertex, CoVertex> vMap = CuttingUtility.cutAtEdge(e);
					for (CoVertex v : vMap.keySet()) {
						CoVertex newV = vMap.get(v);
						newV.setPosition(v.getPosition());
						newV.setSolverIndex(v.getSolverIndex());
					}
				}
			}
			X = hds.numVertices() - hds.numEdges() / 2 + hds.numFaces();
			g = (2 - X) / 2;
			System.err.println("genus of the surface after cutting is " + g);
			
			for (Set<CoEdge> path : context.paths) {
				Set<CoEdge> coPath = new HashSet<CoEdge>();
				context.pathCutMap.put(path, coPath);
				for (CoEdge e : path) {
					coPath.add(context.cutMap.get(e));
				}
			}

		}
		
		Set<CoVertex> visited = new HashSet<CoVertex>(hds.numVertices());
		Queue<CoVertex> Qv = new LinkedList<CoVertex>();
		Queue<CoEdge> Qe = new LinkedList<CoEdge>();
		// start
		CoEdge e0 = hds.getEdge(0);
		for (CoEdge e : hds.getEdges()) { // find an inner edge
			if (HalfEdgeUtils.isInteriorEdge(e)) {
				e0 = e;
			}
		}
		CoEdge e1 = e0.getOppositeEdge();
		CoVertex v1 = e0.getStartVertex();
		CoVertex v2 = e0.getTargetVertex();
		// queued data
		Qv.offer(v1);
		Qv.offer(v2);
		Qe.offer(e1);
		Qe.offer(e0);

		// vertices
		Double d = getNewLength(e0, u);
		
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

//				System.err.println("Visited: " + visited.contains(aVertex) + ", " + visited.contains(bVertex));
				
				if (!visited.contains(cVertex)) {
					d = getNewLength(e, u);
					
					Point A = aVertex.getTextureCoord();
					Point B = bVertex.getTextureCoord();
					
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
					double distACCalc = getNewLength(next, u);
					double dif1 = Math.abs(d1 - distACCalc);
					double dif2 = Math.abs(d2 - distACCalc);
					Point C = dif1 < dif2 ? C1 : C2; 
					double dif = dif1 < dif2 ? dif1 : dif2;
					if (dif < 1E-5) {
						cVertex.setTextureCoord(C);
						visited.add(cVertex);
						Qv.offer(cVertex);
						Qe.offer(e);	
//						System.err.println("Point is valid");
					} else {
//						System.err.println("Point is invalid");
					}
						
//					try {
//						double distAB = Pn.distanceBetween(A.get(), B.get(), Pn.HYPERBOLIC);
//						double distBC = Pn.distanceBetween(B.get(), C.get(), Pn.HYPERBOLIC);
//						double distAC = Pn.distanceBetween(C.get(), A.get(), Pn.HYPERBOLIC);
//						double distAB2 = getNewLength(e.getPreviousEdge(), u);
//						double distBC2 = getNewLength(e, u);
//						
//						System.err.println(e.getLeftFace() + " - (" + aVertex.getIndex() + "," + bVertex.getIndex() + "," + cVertex.getIndex() + ") ------------------------");
//						System.err.println("AB: (" + prev.getIndex() + "," + prev.getOppositeEdge().getIndex() + ")\t" + distAB + "\t" + distAB2);
//						System.err.println("BC: (" + e.getIndex() + "," + e.getOppositeEdge().getIndex() + ")\t" + distBC + "\t" + distBC2);
//						System.err.println("AC: (" + next.getIndex() + "," + next.getOppositeEdge().getIndex() + ")\t" + distAC + "\t" + distACCalc);
//						System.err.println("AC1: (" + next.getIndex() + "," + next.getOppositeEdge().getIndex() + ")\t" + d1 + "\t" + distACCalc);
//						System.err.println("AC2: (" + next.getIndex() + "," + next.getOppositeEdge().getIndex() + ")\t" + d2 + "\t" + distACCalc);
//					} catch (IllegalArgumentException iae) {
//						iae.printStackTrace();
//					}
				} 
				e = e.getOppositeEdge().getNextEdge();
			}
		}
		
		
//		List<CoEdge> eList = new LinkedList<CoEdge>(hds.getEdges());
//		for (CoEdge e : eList) {
//			if (e.isPositive()) {
//				continue;
//			}
//			Point s = e.getStartVertex().getTextureCoord();
//			Point t = e.getTargetVertex().getTextureCoord();
//			double d1 = Pn.distanceBetween(s.get(), t.get(), Pn.HYPERBOLIC);
//			double d2 = getNewLength(e, u);
//			if (Math.abs(d1 - d2) < 1E-3) {
//				continue;
//			}
//			
//			
//			if (e.getLeftFace() != null) {
//				hds.removeFace(e.getLeftFace());
//			}
//			if (e.getRightFace() != null) {
//				hds.removeFace(e.getRightFace());
//			}
//			hds.removeEdge(e.getOppositeEdge());
//			hds.removeEdge(e);
//		}
		
		
		System.err.println("Visited points: " + visited.size() + "/" + hds.numVertices());
		return context;
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
