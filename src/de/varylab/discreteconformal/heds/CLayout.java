package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import geom3d.Point;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
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
		Map<CEdge, Double> absoluteMap = new HashMap<CEdge, Double>();
		Set<CVertex> readyVertices = new HashSet<CVertex>();
		Set<CEdge> readyEdges = new HashSet<CEdge>();
		LinkedList<CEdge> Q = new LinkedList<CEdge>(); 
		
		// start
		CEdge e0 = hds.getEdge(0);
		CEdge e1 = e0.getOppositeEdge();
		CVertex v0 = e0.getStartVertex();
		CVertex v1 = e0.getTargetVertex();
		double u0 = v0.getSolverIndex() >= 0 ? u.get(v0.getSolverIndex()) : 0.0; 
		double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
		double l = exp(0.5 * (u0 + u1)) * e0.getLength(); 
		//edges
		Q.offer(e0);
		Q.offer(e1);
		absoluteMap.put(e0, 0.0);
		absoluteMap.put(e1, PI);
		readyEdges.add(e0);
		readyEdges.add(e1);
		// vertices
		v0.setTextureCoord(new Point(0, 0, 0));
		v1.setTextureCoord(new Point(l, 0, 0));
		readyVertices.add(v0);
		readyVertices.add(v1);
		
		// evolution
		while (!Q.isEmpty()) {
			CEdge e = Q.pop();
			e0 = e.getNextEdge();
			v0 = e0.getStartVertex();
			v1 = e0.getTargetVertex();
			if (e.getLeftFace() == null)
				continue;
			if (readyVertices.contains(v1)) {
				if (!readyEdges.contains(e0.getOppositeEdge()))
					Q.offer(e0.getOppositeEdge());
				continue;
			}
			u0 = v0.getSolverIndex() >= 0 ? u.get(v0.getSolverIndex()) : 0.0; 
			u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
			l = exp(0.5 * (u0 + u1)) * e0.getLength(); 
			double ae = absoluteMap.get(e);
			double a = ae + PI - abs(alphaMap.get(e.getPreviousEdge()));
			geom3d.Vector dif = new geom3d.Vector(cos(a), sin(a), 0).times(l);
			v1.getTextureCoord().set(v0.getTextureCoord()).add(dif);
			Q.offer(e0.getOppositeEdge());
			absoluteMap.put(e0.getOppositeEdge(), a + PI);
			readyEdges.add(e0);
			readyEdges.add(e0.getOppositeEdge());
			readyVertices.add(v1);
			CEdge ne = e0.getNextEdge();
			if (!readyEdges.contains(ne))
				Q.offer(ne);
		}
		assert (readyVertices.size() == hds.numVertices());
		
	}
	
	
}
