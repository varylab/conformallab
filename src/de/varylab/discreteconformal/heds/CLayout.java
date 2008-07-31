package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import geom3d.Point;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import no.uib.cipr.matrix.Vector;

public class CLayout {

	
	/**
	 * Do flat layout for a HDS and a scale vector u
	 * @param hds
	 * @param u
	 */
	public static void doLayout(CHDS hds, Vector u, Map<CEdge, Double>... aMap) {
		Map<CEdge, Double> alphaMap = aMap.length == 0 ? hds.calculateAlphas(u) : aMap[0];
		Set<CVertex> visited = new HashSet<CVertex>();
		Queue<CEdge> Qe = new LinkedList<CEdge>(); 
		Queue<CVertex> Qv = new LinkedList<CVertex>();
		Queue<Double> Qa = new LinkedList<Double>();
		
		// start
		CEdge e0 = hds.getEdge(0);
		CEdge e1 = e0.getOppositeEdge();
		CVertex v1 = e0.getStartVertex();
		CVertex v2 = e0.getTargetVertex();
		double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
		double u2 = v2.getSolverIndex() >= 0 ? u.get(v2.getSolverIndex()) : 0.0; 
		double lambda = e0.getLambda();
		double l = exp(lambda + u1 + u2 ); 
		// queued data
		Qv.offer(v1);
		Qe.offer(e1);
		Qv.offer(v2);
		Qe.offer(e0);
		Qa.offer(PI);
		Qa.offer(0.0);

		// vertices
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
				
				Double alpha = alphaMap.get(e.getNextEdge());
				if (alpha == null) {
					break;
//					alpha = 2*PI - getAngleSum(v, alphaMap);
				}
				
				globalAngle -= alpha;
				
				if (!visited.contains(nearVertex)) {
					visited.add(nearVertex);
					Qv.offer(nearVertex);
					Qe.offer(e);	
					Qa.offer(globalAngle);

					v1 = e.getStartVertex();
					v2 = e.getTargetVertex();
					u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
					u2 = v2.getSolverIndex() >= 0 ? u.get(v2.getSolverIndex()) : 0.0;
					lambda = e.getLambda();
					l = exp(lambda + u1 + u2); 
					geom3d.Vector dif = new geom3d.Vector(cos(globalAngle), sin(globalAngle), 0.0).times(l);
					nearVertex.getTextureCoord().set(tp).add(dif);
				} 
				e = e.getOppositeEdge().getNextEdge();
			}
		}
		assert (visited.size() == hds.numVertices());
	}
	
	
//	
//	/**
//	 * Calculate the angle sum at this vertex. Usually this will be 2PI, but at the boundary
//	 * we sum only the inner angles
//	 * @param v
//	 * @return
//	 */
//	private static Double getAngleSum(CVertex v, Map<CEdge, Double> aMap) {
//		Double r = 0.0;
//		List<CEdge> star = HalfEdgeUtils.incomingEdges(v);
//		for (CEdge e : star) {
//			Double a = aMap.get(e.getPreviousEdge());
//			if (a != null)
//				r += a;
//		}
//		return r;
//	}
//	
//	
//	
	
	
	
}
