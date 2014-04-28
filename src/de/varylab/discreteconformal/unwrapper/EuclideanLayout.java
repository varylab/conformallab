package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.sin;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.logging.Logger;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;

public class EuclideanLayout {

	private static Logger
		log = Logger.getLogger(EuclideanLayout.class.getName());
	
	public static CoVertex doLayout(CoHDS hds, Map<CoEdge, Double> edgeLengths, Map<CoEdge, Double> triangleAngles) {
		Set<CoVertex> visited = new HashSet<CoVertex>(hds.numVertices());
		Queue<CoVertex> Qv = new LinkedList<CoVertex>();
		Queue<CoEdge> Qe = new LinkedList<CoEdge>();
		Queue<Double> Qa = new LinkedList<Double>();
		
		// start at boundary edge if there is one
		CoEdge e0 = hds.getEdge(0);
		Collection<CoEdge> boundary = HalfEdgeUtils.boundaryEdges(hds);
		if (!boundary.isEmpty()) {
			e0 = boundary.iterator().next();
			e0 = e0.getOppositeEdge();
		}
		CoEdge e1 = e0.getOppositeEdge();
		CoVertex v1 = e0.getStartVertex();
		CoVertex v2 = e0.getTargetVertex();
		// queued data
		Qv.offer(v1);
		Qv.offer(v2);
		Qe.offer(e1);
		Qe.offer(e0);
		Qa.offer(PI);
		Qa.offer(0.0);

		// vertices
		Double l = edgeLengths.get(e0);
		v1.T = new double[] {0,0,0,1};
		v2.T = new double[] {l,0,0,1};
		visited.add(v1);
		visited.add(v2);
		
		while (!Qv.isEmpty() && hds.numVertices() > visited.size()) {
			CoVertex v = Qv.poll();
			CoEdge inE = Qe.poll();
			Double a = Qa.poll();
			CoEdge outE = inE.getOppositeEdge();
			double[] tp = v.T;
			
			CoEdge e = inE.getNextEdge();
			Double globalAngle = a + PI;
			while (e != outE) {
				CoVertex nearVertex = e.getTargetVertex();
				
				CoEdge next = e.getNextEdge();
				Double alpha = triangleAngles.get(next);
				if (e.getLeftFace() == null) { // a boundary edge
					alpha = 2*PI - calculateAngleSum(v, triangleAngles);
				}
				
				globalAngle -= alpha;
				
				if (!visited.contains(nearVertex)) {
					visited.add(nearVertex);
					Qv.offer(nearVertex);
					Qe.offer(e);
					Qa.offer(globalAngle);

					l = edgeLengths.get(e);
					double[] dif = {cos(globalAngle), sin(globalAngle), 0.0, 1.0};
					Rn.times(dif, l, dif);
					double[] t = Rn.add(null, tp, dif);
					t[3] = 1.0;
					nearVertex.T = t;
				} 
				e = e.getOppositeEdge().getNextEdge();
			}
		}
		if (visited.size() != hds.numVertices()) {
			log.warning("only " + visited.size() + " of " + hds.numVertices() + " vertices habe been processed");
		}
		return v1;		
	}

	
	
	public static CoVertex doLayout(CoHDS hds, ConformalFunctional<CoVertex, CoEdge, CoFace> fun, Vector u) {
		Map<CoEdge, Double> lMap = getLengthMap(hds, fun, u);
		Map<CoEdge, Double> alphaMap = new HashMap<>();
		for (CoEdge e : hds.getEdges()) {
			alphaMap.put(e, e.getAlpha());
		}
		CoVertex root = doLayout(hds, lMap, alphaMap);
		// projective texture coordinates
		for (CoVertex v : hds.getVertices()) {
			double uv = v.getSolverIndex() < 0 ? 0.0 : u.get(v.getSolverIndex());
			double[] t = v.T;
//			double e = exp( -uv );
			double e = exp( -uv / 2 );
			Pn.dehomogenize(t, t);
			Rn.times(t, e, t);
		}
		return root;
	}
	
	
	
	/**
	 * Calculate the angle sum at this vertex. Usually this will be 2PI, but at the boundary
	 * we sum only the inner angles
	 * @param v
	 * @return the angle sum
	 */
	public static Double calculateAngleSum(CoVertex v) {
		Double r = 0.0;
		List<CoEdge> star = incomingEdges(v);
		for (CoEdge e : star) {
			if (e.getLeftFace() != null) {
				r += e.getPreviousEdge().getAlpha();
			}
		}
		return r;
	}
	
	public static Double calculateAngleSum(CoVertex v, Map<CoEdge, Double> angleMap) {
		Double r = 0.0;
		List<CoEdge> star = incomingEdges(v);
		for (CoEdge e : star) {
			if (e.getLeftFace() != null) {
				r += angleMap.get(e.getPreviousEdge());
			}
		}
		return r;
	}
	
	
	public static Map<CoEdge, Double> getLengthMap(CoHDS hds, ConformalFunctional<CoVertex, CoEdge, CoFace> fun, Vector uVec) {
		MTJDomain u = new MTJDomain(uVec);
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = fun.getNewLength(e, u);
			lMap.put(e, l);
			lMap.put(e.getOppositeEdge(), l);
		}
		return lMap;
	}
	
	public static double getNewLength(CoEdge e, ConformalFunctional<CoVertex, CoEdge, CoFace> fun, Vector uVec) {
		MTJDomain u = new MTJDomain(uVec);
		return fun.getNewLength(e, u);
	}

	
}
