package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import geom3d.Point;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class CHyperbolicLayout {

	
	/**
	 * Do flat layout for a HDS and a metric vector u
	 * @param hds mesh
	 * @param u new metric
	 * @param angleMapParam may be null
	 */
	public static void doLayout(CoHDS hds, Vector u) {
		Set<CoVertex> visited = new HashSet<CoVertex>(hds.numVertices());
		Queue<CoVertex> Qv = new LinkedList<CoVertex>();
		Queue<CoEdge> Qe = new LinkedList<CoEdge>();
		Queue<Double> Qa = new LinkedList<Double>();
		// start
		CoEdge e0 = hds.getEdge(0);
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
		
		Point p0 = new Point(0, 0, 1);
		v1.setTextureCoord(p0);
		Matrix T0 = MatrixBuilder.hyperbolic().translate(l, 0, 0).getMatrix();
		Point p1 = new Point(p0.get());
		T0.transformVector(p1.get());
		v2.setTextureCoord(p1);
		
		visited.add(v1);
		visited.add(v2);
		
		while (!Qv.isEmpty() && hds.numVertices() > visited.size()) {
			CoVertex v = Qv.poll();
			CoEdge inE = Qe.poll();
			Double a = Qa.poll();
			CoEdge outE = inE.getOppositeEdge();
			Point tp = v.getTextureCoord();
			
			CoEdge e = inE.getNextEdge();
			Double globalAngle = a + PI;
			while (e != outE) {
				CoVertex nearVertex = e.getTargetVertex();
				
				CoEdge next = e.getNextEdge();
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
		for (CoVertex v : hds.getVertices()) {
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
