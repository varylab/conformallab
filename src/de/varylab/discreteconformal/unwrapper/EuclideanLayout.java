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

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class EuclideanLayout {

	
	/**
	 * Do flat layout for a HDS and a metric vector u
	 * @param hds mesh
	 * @param u new metric
	 * @param angleMapParam may be null
	 */
	public static CoVertex doLayout(CoHDS hds, Vector u) {
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
		Double l = getNewLength(e0, u);
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
				Double alpha = next.getAlpha();
				if (e.getLeftFace() == null) { // a boundary edge
					alpha = 2*PI - calculateAngleSum(v);
				}
				
				globalAngle -= alpha;
				
				if (!visited.contains(nearVertex)) {
					visited.add(nearVertex);
					Qv.offer(nearVertex);
					Qe.offer(e);	
					Qa.offer(globalAngle);

					l = getNewLength(e, u);
					double[] dif = {cos(globalAngle), sin(globalAngle), 0.0, 1.0};
					Rn.times(dif, l, dif);
					double[] t = Rn.add(null, tp, dif);
					t[3] = 1.0;
					nearVertex.T = t;
				} 
				e = e.getOppositeEdge().getNextEdge();
			}
		}
		
		// projective texture coordinates
		for (CoVertex v : hds.getVertices()) {
			double uv = v.getSolverIndex() < 0 ? 0.0 : u.get(v.getSolverIndex());
			double[] t = v.T;
			double e = exp( -uv );
			Pn.dehomogenize(t, t);
			Rn.times(t, e, t);
		}
		
		assert (visited.size() == hds.numVertices());
		return v1;
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
		double l2 = e.getSolverIndex() >= 0 ? u.get(e.getSolverIndex()) : e.getLambda() + u1 + u2;
		return exp(l2 / 2);
	}
	
	
}
