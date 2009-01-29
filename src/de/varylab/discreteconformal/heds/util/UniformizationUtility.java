package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.acos;
import static java.lang.Math.cosh;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import no.uib.cipr.matrix.Vector;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class UniformizationUtility {

	
	public static class UAdapter {
		
		private Vector
			u = null;
		
		public UAdapter(Vector u) {
			this.u = u;
		}
		
		public double getU(CoVertex v) {
			if (v.getSolverIndex() < 0 || v.getSolverIndex() > u.size() - 1) {
				return 0.0;
			} else {
				return u.get(v.getSolverIndex());
			}
		}
		
	}
	
	
	public static Map<CoEdge, Double> getLengthMap(		
		CoHDS hds, 
		UAdapter u
	) {
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = getNewLength(e, u);
			lMap.put(e, l);
			lMap.put(e.getOppositeEdge(), l);
		}
		return lMap;
	}
	

	public static void reduceToFundamentalPolygon(
		CoHDS hds, 
		CoVertex root,
		UAdapter u
	) {
		Map<CoEdge, Double> lMap = getLengthMap(hds, u);
		int tries = 0;
		while (hds.numVertices() > 50 && tries++ < 30) {
			int X = hds.numVertices() - hds.numEdges() / 2 + hds.numFaces();
			int g = (2 - X) / 2;
			System.out.println("genus of the surface is " + g);
			System.out.println("try " + tries + " - vertices: " + hds.numVertices());
			List<CoVertex> vSet = new ArrayList<CoVertex>(hds.getVertices());
			Collections.sort(vSet, new ValenceComparator());
			for (CoVertex v : vSet) {
				if (v == root) {
					continue;
				}
				List<CoEdge> star = incomingEdges(v);
				int degree = star.size();
				if (degree == 3
						&& star.get(0).getLeftFace() != star.get(1).getLeftFace()
						&& star.get(1).getLeftFace() != star.get(2).getLeftFace()
						&& star.get(2).getLeftFace() != star.get(1).getLeftFace()) {
					removeTrivalentVertex(v);
				} else if (degree > 3) {
					int newDegree = degree;
					for (int i = 0; i < degree - 3 && newDegree > 3; i++) {
						int index = i % star.size();
						flipLengthAndAlphas(star.get(index), lMap);
						star.remove(index);
						newDegree = incomingEdges(v).size();
					}
					if (newDegree == 3) {
						removeTrivalentVertex(v);
					} else {
						System.out.println("could not reduce vertex degree " + degree + " - " + newDegree);
					}
				}
			}
		}
	}
	
	
	private static class ValenceComparator implements Comparator<CoVertex> {
		
		@Override
		public int compare(CoVertex o1, CoVertex o2) {
			int v1 = HalfEdgeUtils.incomingEdges(o1).size();
			int v2 = HalfEdgeUtils.incomingEdges(o2).size();
			return v1 - v2;
		}
		
	}
	
	
	
	public static void flipLengthAndAlphas(
		CoEdge e, 
		Map<CoEdge, Double> lMap
	) {
		double a = lMap.get(e.getNextEdge());
		double aP = lMap.get(e.getOppositeEdge().getPreviousEdge());
		double b = lMap.get(e.getPreviousEdge());
		double bP = lMap.get(e.getOppositeEdge().getNextEdge());
		double alpha = e.getNextEdge().getAlpha();
		double alphaP = e.getOppositeEdge().getPreviousEdge().getAlpha();
		double beta = e.getPreviousEdge().getAlpha();
		double betaP = e.getOppositeEdge().getNextEdge().getAlpha();
		double gamma = e.getAlpha();
		double gammaP = e.getOppositeEdge().getAlpha();
		double betaSum = beta + betaP;
		double alphaSum = alpha + alphaP;
		double cP = arcosh(cosh(a)*cosh(aP) - sinh(a)*sinh(aP)*Math.cos(betaSum));
		double cPCheck = arcosh(cosh(b)*cosh(bP) - sinh(b)*sinh(bP)*Math.cos(alphaSum));
		System.out.println("l: " + cP + " - " + cPCheck);
		double l = (cP + cPCheck) / 2;
		lMap.put(e, l);
		lMap.put(e.getOppositeEdge(), l);
		
		flip(e);
		
		e.setAlpha(betaSum);
		e.getOppositeEdge().setAlpha(alphaSum);
		
		calculateAlpha(e.getNextEdge(), lMap);
		calculateAlpha(e.getPreviousEdge(), lMap);
		calculateAlpha(e.getOppositeEdge().getNextEdge(), lMap);
		calculateAlpha(e.getOppositeEdge().getPreviousEdge(), lMap);
		
		double gammaA = e.getNextEdge().getAlpha();
		double gammaAP = e.getPreviousEdge().getAlpha();
		double gammaB = e.getOppositeEdge().getPreviousEdge().getAlpha();
		double gammaBP = e.getOppositeEdge().getNextEdge().getAlpha();
		System.out.println("gamma: " + (gammaA + gammaB) + " - " + gamma);
		System.out.println("gammaP: " + (gammaAP + gammaBP) + " - " + gammaP);
	}
	
	
	private static void calculateAlpha(
		CoEdge e, 
		Map<CoEdge, Double> lMap
	) {
		double a = lMap.get(e);
		double b = lMap.get(e.getNextEdge());
		double c = lMap.get(e.getPreviousEdge());
		double alpha = acos((cosh(b)*cosh(c) - cosh(a)) / (sinh(b)*sinh(c)));
		e.setAlpha(alpha);
	}
	
	
	
	private static double getNewLength(
		CoEdge e, 
		UAdapter u
	) {
		CoVertex v1 = e.getStartVertex();
		CoVertex v2 = e.getTargetVertex();
		Double u1 = u.getU(v1); 
		Double u2 = u.getU(v2);
		Double la = e.getLambda();
		Double lambdaNew = la + u1 + u2;
		return 2 * arsinh( exp(lambdaNew / 2) );
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void removeTrivalentVertex(V v) {
		HalfEdgeDataStructure<V, E, F> hds = v.getHalfEdgeDataStructure();
		E in1 = v.getIncomingEdge();
		E out1 = in1.getOppositeEdge();
		E outer1 = in1.getPreviousEdge();
		E in2 = out1.getPreviousEdge();
		E out2 = in2.getOppositeEdge();
		E outer2 = in2.getPreviousEdge();
		E in3 = out2.getPreviousEdge();
		E out3 = in3.getOppositeEdge();
		E outer3 = in3.getPreviousEdge();
		
		F f1 = in1.getLeftFace();
		F f2 = in2.getLeftFace();
		F f3 = in3.getLeftFace();

		hds.removeVertex(v);
		hds.removeEdge(in1);
		hds.removeEdge(in2);
		hds.removeEdge(in3);
		hds.removeEdge(out1);
		hds.removeEdge(out2);
		hds.removeEdge(out3);
		hds.removeFace(f2);
		hds.removeFace(f3);
		
		outer1.setLeftFace(f1);
		outer2.setLeftFace(f1);
		outer3.setLeftFace(f1);
		
		outer1.linkNextEdge(outer2);
		outer2.linkNextEdge(outer3);
		outer3.linkNextEdge(outer1);
	}
	
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void flip(E e) {
		F leftFace = e.getLeftFace();
		F rightFace = e.getRightFace();
		if (leftFace == rightFace)
			return;
		E a1 = e.getOppositeEdge().getNextEdge();
		E a2 = a1.getNextEdge();
		E b1 = e.getNextEdge();
		E b2 = b1.getNextEdge();
		
		V v1 = e.getStartVertex();
		V v2 = e.getTargetVertex();
		V v3 = a1.getTargetVertex();
		V v4 = b1.getTargetVertex();

		//new connections
		e.linkNextEdge(a2);
		e.linkPreviousEdge(b1);
		e.getOppositeEdge().linkNextEdge(b2);
		e.getOppositeEdge().linkPreviousEdge(a1);
		e.setTargetVertex(v3);
		e.getOppositeEdge().setTargetVertex(v4);
		
		a2.linkNextEdge(b1);
		b2.linkNextEdge(a1);
		
		//set faces
		b2.setLeftFace(rightFace);
		a2.setLeftFace(leftFace);
		b2.setTargetVertex(v1);
		a2.setTargetVertex(v2);
		a1.setTargetVertex(v3);
		b1.setTargetVertex(v4);
	}
	
	
	
	private static double arsinh(double x) {
		double r = x + sqrt(x*x + 1);
		return log(r);
	}
	
	private static double arcosh(double x) {
		double r = x + sqrt(x*x - 1);
		return log(r);
	}
	
}
