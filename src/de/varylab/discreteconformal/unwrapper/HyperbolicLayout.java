package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.sin;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.adapter.LengthMapWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.NodeComparator;
import de.varylab.discreteconformal.util.PathUtility;
import de.varylab.discreteconformal.util.Search;

public class HyperbolicLayout {

	
	public static CoVertex guessRootVertex(
		CoHDS hds, 
		Map<CoEdge, Double> lMap,
		int mcSamples
	) {
		Random rnd = new Random();
		
		Set<CoVertex> boundary = new TreeSet<CoVertex>(new NodeComparator<CoVertex>());
		boundary.addAll(HalfEdgeUtils.boundaryVertices(hds));
		Set<CoVertex> mcBoundarySet = new HashSet<CoVertex>();
		Iterator<CoVertex> boundaryIterator = boundary.iterator();
		for (int i = 0; i < max(boundary.size() / 10, 10); i++) {
			mcBoundarySet.add(boundaryIterator.next());
		}
		
		LengthMapWeightAdapter wa = new LengthMapWeightAdapter(lMap);
		Map<CoVertex, Double> sMap = new HashMap<CoVertex, Double>();
		
		Set<CoVertex> mcSet = new HashSet<CoVertex>();
		for (int i = 0; i < mcSamples; i++) {
			int sampleIndex = rnd.nextInt(hds.numVertices());
			CoVertex sampleVertex = hds.getVertex(sampleIndex);
			if (!HalfEdgeUtils.isBoundaryVertex(sampleVertex)) {
				mcSet.add(sampleVertex);				
			}
		}
		
		for (CoVertex v : mcSet) {
			double mean = 0;
			Map<CoVertex, Double> distMap = new HashMap<CoVertex, Double>();
			Map<CoVertex, List<CoEdge>> pathMap = Search.getAllShortestPaths(v, boundary, wa, new HashSet<CoVertex>());
			for (CoVertex bv : mcBoundarySet) {
				List<CoEdge> path = pathMap.get(bv);
				double length = PathUtility.getTotalPathWeight(new HashSet<CoEdge>(path), wa);
				mean += length;
				distMap.put(bv, length);
			}
			mean /= boundary.size();
			double s = 0.0;
			for (CoVertex bv : distMap.keySet()) {
				double dist = distMap.get(bv);
				s += (mean - dist) * (mean - dist);
			}
			s /= boundary.size();
			sMap.put(v, s);
		}
		CoVertex root = mcSet.iterator().next();
		double minS = sMap.get(root); 
		for (CoVertex v : sMap.keySet()) {
			double s = sMap.get(v);
			if (s < minS) {
				minS = s;
				root = v;
			}
		}
		
		return root;
	}
	
	
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
	public static CoVertex doLayout(CoHDS hds, CoVertex root, Vector u) {
//		System.out.println("Layout --------------------");
//		for (CoVertex v : hds.getVertices()) {
//			System.out.println("sum " + v.getIndex() + ": " + getAngleSum(v));
//		}
		final Map<CoEdge, Double> lMap = getLengthMap(hds, u);
		
		final Set<CoVertex> visited = new HashSet<CoVertex>(hds.numVertices());
		final Queue<CoVertex> Qv = new LinkedList<CoVertex>();
		final Queue<CoEdge> Qe = new LinkedList<CoEdge>();
		// start
		final CoVertex v1 = root;//guessRootVertex(hds, lMap, 10);
		final CoEdge e1 = v1.getIncomingEdge();
		final CoEdge e0 = e1.getOppositeEdge();
		final CoVertex v2 = e0.getTargetVertex();
		// queued data
		Qv.offer(v1);
		Qv.offer(v2);
		Qe.offer(e1);
		Qe.offer(e0);

		// vertices
		Double d = lMap.get(e0);
		
		v1.T = new double[] {0, 0, 0, 1};
		v2.T = new double[] {sinh(d), 0, 0, cosh(d)};
		Pn.normalize(v2.T, v2.T, Pn.HYPERBOLIC);
		
		visited.add(v1);
		visited.add(v2);
		
		while (!Qv.isEmpty() && hds.numVertices() > visited.size()) {
			final CoVertex v = Qv.poll();
			final CoEdge inE = Qe.poll();
			final CoEdge outE = inE.getOppositeEdge();
			
			CoEdge e = inE.getNextEdge();
			while (e != outE) {
				final CoEdge next = e.getNextEdge();
				final CoEdge prev = e.getPreviousEdge();
				final CoVertex aVertex = prev.getStartVertex();
				final CoVertex bVertex = prev.getTargetVertex();
				final CoVertex cVertex = e.getTargetVertex();

				Double alpha = next.getAlpha();
				if (e.getLeftFace() == null) { // a boundary edge
					alpha = 2*PI - getAngleSum(v);
					e = e.getOppositeEdge().getNextEdge();
					continue;
				}
				if (!visited.contains(cVertex)) {
					d = lMap.get(e);
					double dCheck = lMap.get(next);
					double[] A = aVertex.T;
					double[] B = bVertex.T;
					
					double[] C = null;
					try {
						C = layoutTriangle(A, B, alpha, d, dCheck);
					} catch (Exception e2) {
						System.out.print(".");
						//e2.printStackTrace();
					}
					if (C != null) {
						cVertex.T = C;
						visited.add(cVertex);
						Qv.offer(cVertex);
						Qe.offer(e);
					}
				}
				e = e.getOppositeEdge().getNextEdge();
			}
		}
		
		return v1;
	}
	
	
	
	private static double[] layoutTriangle(double[] A, double[] B, double alpha, double d, double dP) throws Exception {
		// calculation is in RP2
		// project to RP2 
		double[] A3 = {A[0], A[1], A[3]};
		double[] B3 = {B[0], B[1], B[3]};
		// polarize
		double[] AHat = {A[0], A[1], -A[3]};
		double[] BHat = {B[0], B[1], -B[3]};
		double[] lAB = Rn.crossProduct(null, A3, B3);
		normalize(lAB);
		double[] At = Rn.crossProduct(null, lAB, BHat);
		normalize(At);
		double[] AtPerp = Rn.crossProduct(null, AHat, BHat);
		normalize(AtPerp);
		double[] Ct = Rn.linearCombination(null, cos(alpha), At, sin(alpha), AtPerp);
		normalize(Ct);
		double[] C1 = Rn.linearCombination(null, cosh(d), B3, sinh(d), Ct);
		normalize(C1);
		double[] C2 = Rn.linearCombination(null, cosh(d), B3, -sinh(d), Ct);
		normalize(C2);
		double d1 = Double.MAX_VALUE;
		double d2 = Double.MAX_VALUE;
		try {
			d1 = Pn.distanceBetween(C1, A3, Pn.HYPERBOLIC);
		} catch (IllegalArgumentException iae) {}
		try {
			d2 = Pn.distanceBetween(C2, A3, Pn.HYPERBOLIC);
		} catch (IllegalArgumentException iae) {}
		double dif1 = Math.abs(d1 - dP);
		double dif2 = Math.abs(d2 - dP);
		double[] C = dif1 < dif2 ? C1 : C2; 
		double dif = dif1 < dif2 ? dif1 : dif2;
		if (dif < 1E-5) {
			// lift to RP3
			return new double[] {C[0], C[1], 0, C[2]};
		} else {
			return null;
		}
	}
	
	
	
	
	private static double[] normalize(double[] p) {
		Pn.normalize(p, p, Pn.HYPERBOLIC);
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