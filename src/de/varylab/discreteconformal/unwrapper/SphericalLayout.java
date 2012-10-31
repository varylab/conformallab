package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.tan;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import no.uib.cipr.matrix.Vector;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.geometry.ComplexProjective1;
import de.jtem.mfc.group.Moebius;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;

public class SphericalLayout {

	public static CoVertex doLayout(CoHDS hds, CoVertex root, ConformalFunctional<CoVertex, CoEdge, CoFace> fun, Vector u) {
		final Map<CoEdge, Double> lMap = getLengthMap(hds, fun, u);
		
		final Set<CoVertex> visited = new HashSet<CoVertex>(hds.numVertices());
		final Queue<CoVertex> Qv = new LinkedList<CoVertex>();
		final Queue<CoEdge> Qe = new LinkedList<CoEdge>();
		// start
		final CoVertex v1 = root;
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
		
		// complex homogeneous coordinates
		v1.T = new double[] {0, 0, 1, 0};
		v2.T = new double[] {tan(d/2), 0, 1, 0};
		
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
					double dd = lMap.get(prev);
					double[] A = aVertex.T;
					double[] B = bVertex.T;
					double[] C = layoutTriangle(A, B, alpha, d, dd, dCheck);
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
		
		if (visited.size() != hds.numVertices()) {
			throw new RuntimeException("Skipped vertices during layout.");
		}
		
		for (CoVertex v : hds.getVertices()){
			ComplexProjective1 tp = new ComplexProjective1(v.T[0], v.T[1], v.T[2], v.T[3]);
			Complex tc = new Complex();
			tp.projectTo(tc);
			double x = tc.re;
			double y = tc.im;
			double nx = 2 * x;
			double ny = x*x + y*y - 1;
			double nz = 2 * y;
			double nw = ny + 2;
			v.T = new double[] {nx / nw, ny / nw, nz / nw, 1.0};
		}
		
		return v1;
	}
	
	
	
	private static double[] layoutTriangle(double[] A, double[] B, double alpha, double d, double dd, double dP) {
		double logScale = Math.log(d / dd);
		Moebius M = new Moebius();
		ComplexProjective1 Ac = new ComplexProjective1(A[0], A[1], A[2], A[3]);
		ComplexProjective1 Bc = new ComplexProjective1(B[0], B[1], B[2], B[3]);
		M.assignSphericalLogScaleRotation(Bc, logScale, alpha);
		ComplexProjective1 Dc = new ComplexProjective1(); 
		M.applyTo(Ac, Dc);
		return new double[] {Dc.aRe, Dc.aIm, Dc.bRe, Dc.bIm};
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
