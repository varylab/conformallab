package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryEdge;
import static java.lang.Math.cosh;
import static java.lang.Math.sinh;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.adapter.LengthMapWeightAdapter;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.util.NodeIndexComparator;
import de.varylab.discreteconformal.util.PathUtility;
import de.varylab.discreteconformal.util.Search;

public class HyperbolicLayout {

	private static Logger
		log = Logger.getLogger(HyperbolicLayout.class.getName());
	
	/**
	 * Implements a heuristic for finding a root vertex for the
	 * hyperbolic layout. Its numerically best to have equal distance to the boundary. 
	 * @param hds
	 * @param lMap
	 * @param mcSamples
	 * @return
	 */
	public static CoVertex guessRootVertex(
		CoHDS hds, 
		Map<CoEdge, Double> lMap,
		int mcSamples
	) {
		// deterministic random numbers
		Random rnd = new Random(0);
		rnd.setSeed(hds.numVertices());
		
		Set<CoVertex> boundary = new TreeSet<CoVertex>(new NodeIndexComparator<CoVertex>());
		boundary.addAll(HalfEdgeUtils.boundaryVertices(hds));
//		Set<CoVertex> mcBoundarySet = new HashSet<CoVertex>();
//		Iterator<CoVertex> boundaryIterator = boundary.iterator();
//		for (int i = 0; i < max(boundary.size() / 5, 5); i++) {
//			mcBoundarySet.add(boundaryIterator.next());
//		}
		
		LengthMapWeightAdapter wa = new LengthMapWeightAdapter(lMap);
		Map<CoVertex, Double> sMap = new HashMap<CoVertex, Double>();
		
		Set<CoVertex> mcSet = new HashSet<CoVertex>();// TreeSet<CoVertex>(new NodeIndexComparator<CoVertex>());
		for (int i = 0; i < Math.min(mcSamples, hds.numVertices()); i++) {
			int sampleIndex = rnd.nextInt(hds.numVertices());
			CoVertex sampleVertex = hds.getVertex(sampleIndex);
			if (!HalfEdgeUtils.isBoundaryVertex(sampleVertex)) {
				mcSet.add(sampleVertex);				
			}
		}
		
		if (mcSet.isEmpty()) {
			// no interior vertex
			return hds.getVertex(0);
		}
		
		for (CoVertex v : mcSet) {
			double mean = 0;
			Map<CoVertex, Double> distMap = new HashMap<CoVertex, Double>();
			Map<CoVertex, List<CoEdge>> pathMap = Search.getAllShortestPaths(v, boundary, wa, new HashSet<CoVertex>());
			for (CoVertex bv : boundary) {
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
		double minS = Double.MAX_VALUE; 
		for (CoVertex v : sMap.keySet()) {
			double s = sMap.get(v);
			if (s < minS) {
				minS = s;
				root = v;
			}
		}
		return root;
	}
	
	public static CoVertex doLayout(CoHDS hds, CoVertex root, ConformalFunctional<CoVertex, CoEdge, CoFace> fun, Vector u) {
		Map<CoEdge, Double> lMap = getLengthMap(hds, fun, u);
		return doLayout(hds, root, lMap);
	}
	
	public static CoVertex doLayout(CoHDS hds, CoVertex root, Map<CoEdge, Double> lMap) {
		final Set<CoFace> visitedFaces = new HashSet<CoFace>();
		final Set<CoVertex> visitedVertices = new HashSet<CoVertex>();
		final Queue<CoEdge> Qe = new LinkedList<CoEdge>();
		// start
		final CoVertex v1;
		if (root != null) {
			v1 = root;
		} else {
			v1 = guessRootVertex(hds, lMap, 200);
		}
		log.info("layout root is " + v1.getIndex());
		
		// get a deterministic first edge
		CoEdge e1 = v1.getIncomingEdge();
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v1)) {
			if (e.getIndex() < e1.getIndex()) {
				e1 = e;
			}
		}
		final CoEdge e0 = e1.getOppositeEdge();
		final CoVertex v2 = e0.getTargetVertex();
		// queued data
		Qe.offer(e1);
		Qe.offer(e0);
		visitedVertices.add(v1);
		visitedVertices.add(v2);

		// vertices
		double d = lMap.get(e0);
		v1.T = new double[] {0, 0, 0, 1};
		v2.T = new double[] {sinh(d), 0, 0, cosh(d)};
		
		while (!Qe.isEmpty()) {
			final CoEdge ab = Qe.poll();
			final CoEdge bc = ab.getNextEdge();
			final CoEdge ac = ab.getPreviousEdge();
			final CoFace f = ab.getLeftFace();
			if (visitedFaces.contains(f)) continue;
			final CoVertex a = ab.getStartVertex();
			final CoVertex b = ab.getTargetVertex();
			final CoVertex c = bc.getTargetVertex();
			final double alpha = ac.getAlpha();
			double dAB = lMap.get(ab);
			double dBC = lMap.get(bc);
			double dAC = lMap.get(ac);
			double[] C = layoutTriangle(a.T, b.T, alpha, dAB, dBC, dAC);
			if (C == null) {
				log.warning("layout at face " + f + " skipped");
				continue;
			}
			c.T = C;
			visitedFaces.add(f);
			visitedVertices.add(c);
			if (!(isBoundaryEdge(bc) || visitedFaces.contains(bc.getRightFace()))) {
				Qe.offer(bc.getOppositeEdge());
			}
			if (!(isBoundaryEdge(ac) || visitedFaces.contains(ac.getRightFace()))) {
				Qe.offer(ac.getOppositeEdge());
			}
		}
		
		if (visitedFaces.size() != hds.numFaces()) {
			log.warning("only " + visitedFaces.size() + " of " + hds.numFaces() + " faces have been visited");
		}
		if (visitedVertices.size() != hds.numVertices()) {
			log.warning("only " + visitedVertices.size() + " of " + hds.numVertices() + " vertices have been visited");
		}
		return v1;
	}
	
	
	private static double[]
		ZERO = {0,0,0,1};
	
	/*
	 * perform hyperbolic stretch rotation
	 */
	static double[] layoutTriangle(double[] A, double[] B, double alpha, double dAB, double dBC, double dAC) {
		MatrixBuilder mb = MatrixBuilder.hyperbolic();
		mb.translateFromTo(ZERO, B);
		mb.rotate(alpha, 0, 0, 1);
		mb.scale((sinh(dBC)/cosh(dBC)) / (sinh(dAB)/cosh(dAB)));
		mb.translate(B, ZERO);
		double[] C = A.clone();
		mb.getMatrix().transformVector(C);
		Pn.normalize(C, C, Pn.HYPERBOLIC);
		double dACcheck = Pn.distanceBetween(A, C, Pn.HYPERBOLIC);
		if (Math.abs(dACcheck - dAC) > 1E-5) {
			return null;
		} else {
			return C;
		}
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
	
	
	public static Map<CoEdge, Double> getLengthMap(CoHDS hds, ConformalFunctional<CoVertex, CoEdge, CoFace> fun, Vector uVec) {
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		MTJDomain u = new MTJDomain(uVec);
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = fun.getNewLength(e, u);
			lMap.put(e, l);
			lMap.put(e.getOppositeEdge(), l);
		}
		return lMap;
	}
	
	/**
	 * Calculate the edge length for the flat metric
	 * @param e
	 * @param u
	 * @return the new edge length
	 */
	public static Double getNewLength(CoEdge e, ConformalFunctional<CoVertex, CoEdge, CoFace> fun, Vector uVec) {
		MTJDomain u = new MTJDomain(uVec);
		return fun.getNewLength(e, u);
	}
	
}
