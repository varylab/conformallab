package de.varylab.discreteconformal.heds.unwrap;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import geom3d.Point;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;

public class CDiskLayout {

	
	/**
	 * Do flat layout for a HDS and a metric vector u
	 * @param hds mesh
	 * @param u new metric
	 * @param angleMapParam may be null
	 */
	public static void doLayout(CHDS hds, Vector u) {
		Set<CVertex> visited = new HashSet<CVertex>(hds.numVertices());
		Queue<CVertex> Qv = new LinkedList<CVertex>();
		Queue<CEdge> Qe = new LinkedList<CEdge>();
		Queue<Double> Qa = new LinkedList<Double>();
		// start
		CEdge e0 = hds.getEdge(0);
		CEdge e1 = e0.getOppositeEdge();
		CVertex v1 = e0.getStartVertex();
		CVertex v2 = e0.getTargetVertex();
		// queued data
		Qv.offer(v1);
		Qv.offer(v2);
		Qe.offer(e1);
		Qe.offer(e0);
		Qa.offer(PI);
		Qa.offer(0.0);

		// vertices
		Double l = getNewLength(e0, u);
		v1.setTextureCoord(new Point(0, 0, 0));
		v2.setTextureCoord(new Point(l, 0, 0));
		visited.add(v1);
		visited.add(v2);
		
		while (!Qv.isEmpty() && hds.numVertices() > visited.size()) {
			CVertex v = Qv.poll();
			CEdge inE = Qe.poll();
			Double a = Qa.poll();
			CEdge outE = inE.getOppositeEdge();
			Point tp = v.getTextureCoord();
			
			CEdge e = inE.getNextEdge();
			Double globalAngle = a + PI;
			while (e != outE) {
				CVertex nearVertex = e.getTargetVertex();
				
				CEdge next = e.getNextEdge();
				Double alpha = next.getAlpha();
				if (e.getLeftFace() == null) { // a boundary edge
					alpha = 2*PI - getAngleSum(v);
				}
				
				globalAngle -= alpha;
				
				if (!visited.contains(nearVertex)) {
					visited.add(nearVertex);
					Qv.offer(nearVertex);
					Qe.offer(e);	
					Qa.offer(globalAngle);

					l = getNewLength(e, u);
					geom3d.Vector dif = new geom3d.Vector(cos(globalAngle), sin(globalAngle), 0.0).times(l);
					nearVertex.getTextureCoord().set(tp).add(dif);
				} 
				e = e.getOppositeEdge().getNextEdge();
			}
		}
		
		// projective texture coordinates
		for (CVertex v : hds.getVertices()) {
			double uv = v.getSolverIndex() < 0 ? 0.0 : u.get(v.getSolverIndex());
			Point t = v.getTextureCoord();
			double e = exp( -uv );
			t.set(e * t.x(), e * t.y(), e);
		}
		
		assert (visited.size() == hds.numVertices());
	}
	
	
	
	/**
	 * Calculate the angle sum at this vertex. Usually this will be 2PI, but at the boundary
	 * we sum only the inner angles
	 * @param v
	 * @return the angle sum
	 */
	public static Double getAngleSum(CVertex v) {
		Double r = 0.0;
		List<CEdge> star = incomingEdges(v);
		for (CEdge e : star) {
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
	public static Double getNewLength(CEdge e, Vector u) {
		CVertex v1 = e.getStartVertex();
		CVertex v2 = e.getTargetVertex();
		Double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
		Double u2 = v2.getSolverIndex() >= 0 ? u.get(v2.getSolverIndex()) : 0.0;
		Double lambda = e.getLambda();
		return exp(lambda + u1 + u2);
	}
	
	
}
