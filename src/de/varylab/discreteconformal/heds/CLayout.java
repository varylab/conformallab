package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
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
		doLayout2(hds, u, aMap);
//		Map<CEdge, Double> alphaMap = aMap.length == 0 ? hds.calculateAlphas(u) : aMap[0];
//		Map<CEdge, Double> absoluteMap = new HashMap<CEdge, Double>();
//		Set<CVertex> readyVertices = new HashSet<CVertex>();
//		Set<CEdge> readyEdges = new HashSet<CEdge>();
//		Queue<CEdge> Q = new LinkedList<CEdge>(); 
//		
//		// start
//		CEdge e0 = hds.getEdge(0);
//		CEdge e1 = e0.getOppositeEdge();
//		CVertex v0 = e0.getStartVertex();
//		CVertex v1 = e0.getTargetVertex();
//		double u0 = v0.getSolverIndex() >= 0 ? u.get(v0.getSolverIndex()) : 0.0; 
//		double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
//		double l = exp(0.5 * (u0 + u1)) * e0.getLength(); 
//		//edges
//		Q.offer(e0);
//		Q.offer(e1);
//		absoluteMap.put(e0, 0.0);
//		absoluteMap.put(e1, PI);
//		readyEdges.add(e0);
//		readyEdges.add(e1);
//		// vertices
//		v0.setTextureCoord(new Point(0, 0, 0));
//		v1.setTextureCoord(new Point(l, 0, 0));
//		readyVertices.add(v0);
//		readyVertices.add(v1);
//		
//		// evolution
//		while (!Q.isEmpty()) {
//			CEdge e = Q.poll();
//			e0 = e.getNextEdge();
//			v0 = e0.getStartVertex();
//			v1 = e0.getTargetVertex();
//			if (e.getLeftFace() == null)
//				continue;
//			if (!readyEdges.contains(e0.getOppositeEdge()))
//				Q.offer(e0.getOppositeEdge());
//			if (readyVertices.contains(v1)) {
//				continue;
//			}
//			u0 = v0.getSolverIndex() >= 0 ? u.get(v0.getSolverIndex()) : 0.0; 
//			u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
//			l = exp(0.5 * (u0 + u1)) * e0.getLength(); 
//			double ae = absoluteMap.get(e);
//			double a = ae + PI - abs(alphaMap.get(e.getPreviousEdge()));
//			geom3d.Vector dif = new geom3d.Vector(cos(a), sin(a), 0).times(l);
//			v1.getTextureCoord().set(v0.getTextureCoord()).add(dif);
//			Q.offer(e0.getOppositeEdge());
//			absoluteMap.put(e0.getOppositeEdge(), a + PI);
//			readyEdges.add(e0);
//			readyEdges.add(e0.getOppositeEdge());
//			readyVertices.add(v1);
//			CEdge ne = e0.getNextEdge();
//			if (!readyEdges.contains(ne))
//				Q.offer(ne);
//		}
//		assert (readyVertices.size() == hds.numVertices());
//		
	}
	
	
	
	public static void doLayout2(CHDS hds, Vector u, Map<CEdge, Double>... aMap) {
		Map<CEdge, Double> alphaMap = aMap.length == 0 ? hds.calculateAlphas(u) : aMap[0];
		Set<CVertex> visited = new HashSet<CVertex>();
		Queue<CEdge> Qe = new LinkedList<CEdge>(); 
		Queue<CVertex> Qv = new LinkedList<CVertex>();
		Queue<Double> Qa = new LinkedList<Double>();
		
		// start
		CEdge e0 = hds.getEdge(0);
		CEdge e1 = e0.getOppositeEdge();
		CVertex v0 = e0.getStartVertex();
		CVertex v1 = e0.getTargetVertex();
		double u0 = v0.getSolverIndex() >= 0 ? u.get(v0.getSolverIndex()) : 0.0; 
		double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
		double l = exp(0.5 * (u0 + u1)) * e0.getLength(); 
		// queued data
		Qv.offer(v0);
		Qe.offer(e1);
		Qv.offer(v1);
		Qe.offer(e0);
		Qa.offer(PI);
		Qa.offer(0.0);

		// vertices
		v0.setTextureCoord(new Point(0, 0, 0));
		v1.setTextureCoord(new Point(l, 0, 0));
		visited.add(v0);
		visited.add(v1);
			
		while (!Qv.isEmpty()) {
			CVertex v = Qv.poll();
			CEdge inE = Qe.poll();
			Double a = Qa.poll();
			CEdge outE = inE.getOppositeEdge();
			Point tp = v.getTextureCoord();
			
			
			System.out.println("Reached Vertex " + v + ": " + tp);
			
			CEdge e = inE.getNextEdge();
			Double localAngle = a;
			while (e != outE) {
//				Double alpha = alphaMap.get(e.getNextEdge());
//				if (alpha == null) {
//					e = e.getOppositeEdge().getNextEdge();
//					continue;
//				}
//				localAngle += PI - alpha;
				CVertex nextVertex = e.getTargetVertex();
				if (visited.contains(nextVertex)) {
					e = e.getOppositeEdge().getNextEdge();
					continue;
				}
//				v0 = e.getTargetVertex();
//				v1 = e.getStartVertex();
//				u0 = v0.getSolverIndex() >= 0 ? u.get(v0.getSolverIndex()) : 0.0; 
//				u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
//				l = exp(0.5 * (u0 + u1)) * e.getLength(); 
//				geom3d.Vector dif = new geom3d.Vector(cos(localAngle), sin(localAngle), 0.0).times(l);
//				nextVertex.getTextureCoord().set(tp).add(dif);
//				
				
				visited.add(nextVertex);
				Qv.offer(nextVertex);
				Qe.offer(e.getOppositeEdge());	
				Qa.offer(localAngle);
				e = e.getOppositeEdge().getNextEdge();
			}
			
		}
		
	}
	
	
}
