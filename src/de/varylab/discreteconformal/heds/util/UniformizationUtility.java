package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cosh;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;
import static java.util.logging.Level.FINE;
import static java.util.logging.Level.FINEST;
import static java.util.logging.Level.OFF;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

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

	private static Logger 
		log = Logger.getLogger(UniformizationUtility.class.getName());
	
	static {
		log.setLevel(FINE);
	}
	
	
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
		while (hds.numVertices() > 1 && tries++ < 50) {
			int X = hds.numVertices() - hds.numEdges() / 2 + hds.numFaces();
			int g = (2 - X) / 2;
			
			log.log(FINE, "genus of the surface is " + g);
			log.log(FINE, "try " + tries + " - vertices: " + hds.numVertices());
			List<CoVertex> vSet = new ArrayList<CoVertex>(hds.getVertices());
			Collections.sort(vSet, new ValenceComparator());
			for (CoVertex v : vSet) {
				log.log(FINEST, "try to remove vertex " + v + " ---------------------");
				if (v == root) {
					log.log(FINEST, "is root, can't remove");
					continue; 
				}
				List<CoEdge> star = incomingEdges(v);
				int degree = star.size();
				log.log(FINEST, "vertex degree is " + degree);
				if (degree == 3) {
					assert (star.get(0).getLeftFace() != star.get(1).getLeftFace()
						&& star.get(1).getLeftFace() != star.get(2).getLeftFace()
						&& star.get(2).getLeftFace() != star.get(1).getLeftFace());
					removeTrivalentVertex(v);
				} else if (degree > 3) {
					int newDegree = degree;
					int flipTries = 0; 
					CoEdge flipEdge = v.getIncomingEdge();
					while (newDegree > 3 && flipTries++ < 30) {
						CoEdge nextFlipEdge = flipEdge.getNextEdge().getOppositeEdge();
						log.log(FINEST, "trying to flip edge " + flipEdge);
						if (!isFlippable(flipEdge)) {
							log.log(FINEST, "cannot flip: angle condition");
							flipEdge = nextFlipEdge;
							continue;
						}
						assert flipEdge.getLeftFace() != flipEdge.getRightFace();
						flipLengthAndAlphas(flipEdge, lMap);
						newDegree = incomingEdges(v).size();
						flipEdge = nextFlipEdge;
						log.log(FINEST, "new vertex degree is " + newDegree);
					}
					if (newDegree == 3) {
						star = incomingEdges(v);
						assert (star.get(0).getLeftFace() != star.get(1).getLeftFace()
							&& star.get(1).getLeftFace() != star.get(2).getLeftFace()
							&& star.get(2).getLeftFace() != star.get(1).getLeftFace());
						removeTrivalentVertex(v);
					} else {
						log.log(FINEST, "could not reduce vertex degree " + degree + " - " + newDegree);
					}
				}
			}
		}
		log.log(FINE, "tries: " + tries);
		for (CoVertex v : hds.getVertices()) {
			double alpha = 0.0;
			for (CoEdge e : incomingEdges(v)) {
				alpha += e.getPreviousEdge().getAlpha();
			}
			log.log(FINEST, "alpha at " + v + ": " + alpha);
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
	
	
	
	public static boolean isFlippable(CoEdge e) {
		double alpha = 0.0;
		alpha += e.getNextEdge().getAlpha();
		alpha += e.getOppositeEdge().getPreviousEdge().getAlpha();
		double beta = 0.0;
		beta += e.getPreviousEdge().getAlpha();
		beta += e.getOppositeEdge().getNextEdge().getAlpha();
		return alpha < PI && beta < PI;
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
//		log.log(FINEST, "l: " + cP + " - " + cPCheck);
		double l = (cP + cPCheck) / 2;
		lMap.put(e, l);
		lMap.put(e.getOppositeEdge(), l);
		
		log.log(FINEST, "flip edge " + e);
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
//		log.log(FINEST, "gamma: " + (gammaA + gammaB) + " - " + gamma);
//		log.log(FINEST, "gammaP: " + (gammaAP + gammaBP) + " - " + gammaP);
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
		log.log(FINEST, "remove vertex " + v);
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
