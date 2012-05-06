package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.varylab.discreteconformal.unwrapper.IsothermicUtility.calculateTriangleAngle;
import static java.lang.Math.sin;

import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddGrad;
import de.jtem.jtao.TaoAppAddObjective;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.unwrapper.IsothermicUtility;

public class SinConditionFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	HDS extends HalfEdgeDataStructure<V, E, F>
> extends TaoApplication implements TaoAppAddGrad, TaoAppAddObjective {

	private HDS
		hds = null;
	private Map<Integer, Integer>
		solverIndices = null;
	
	public SinConditionFunctional(HDS hds, Map<Integer, Integer> undirectedIndexMap) {
		this.hds = hds;
		this.solverIndices = undirectedIndexMap;
	}
	
	public void calculateAndSetInitionSolution(AdapterSet a) {
		Vec alpha = new Vec(hds.numEdges() / 2);
		for (E e : hds.getEdges()) {
			int index = solverIndices.get(e.getIndex());
			double[] N = a.getD(Normal.class, e);
			double[] Kmin = a.getD(CurvatureFieldMin.class, e);
			double[] E = a.getD(EdgeVector.class, e);
			double ae = IsothermicUtility.getSignedAngle(N, Kmin, E);
			alpha.setValue(index, ae, InsertMode.INSERT_VALUES);
		}
		alpha.assemble();
		setInitialSolutionVec(alpha);
	}
	
	@Override
	public void evaluateGradient(Vec x, Vec g) {
		defaultComputeGradient(x, g);
	}
	
	
//	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
//		for (V v : hds.getVertices()) {
//			if (HalfEdgeUtils.isBoundaryVertex(v)) {
//				continue;
//			}
//			// star alphas
//			for (E e : HalfEdgeUtils.incomingEdges(v)) {
//				double sl = getLeftSinProduct(v, x);
//				double sr = getRightSinProduct(v, x);
//				double eAlpha = 2 * (sl - sr);
//				double ol = 1.0;
//				double or = 1.0;
//				for (E e2 : HalfEdgeUtils.incomingEdges(v)) {
//					// betas left and right of e2
//					double betaL = getOppsiteBeta(e.getNextEdge(), x);
//					double betaR = getOppsiteBeta(e.getOppositeEdge().getPreviousEdge(), x);
//					if (e2 == e) {
//						ol *= cos(betaL);
//						or *= cos(betaR);
//					} else {
//						ol *= sin(betaL);
//						or *= sin(betaR);
//					}
//				}
//				eAlpha *= (ol - or);
//				int index = solverIndices.get(e.getIndex());
//				g.add(index, eAlpha);
//			}
//			// rim alphas
//			for (E e : HalfEdgeUtils.incomingEdges(v)) {
//				double sl = getLeftSinProduct(v, x);
//				double sr = getRightSinProduct(v, x);
//				double eAlpha = 2 * (sl - sr);
//				double ol = 1.0;
//				double or = 1.0;
//				for (E e2 : HalfEdgeUtils.incomingEdges(v)) {
//					// betas left and right of e2
//					double betaL = getOppsiteBeta(e.getNextEdge(), x);
//					double betaR = getOppsiteBeta(e.getOppositeEdge().getPreviousEdge(), x);
//					if (e2 == e) {
//						ol *= cos(betaL);
//						or *= cos(betaR);
//					} else {
//						ol *= sin(betaL);
//						or *= sin(betaR);
//					}
//				}
//				eAlpha *= (ol - or);
//				int index = solverIndices.get(e.getPreviousEdge().getIndex());
//				g.add(index, eAlpha);
//			}
//		}
//		return evaluateObjective(x);
//	}
	
	@Override
	public double evaluateObjective(Vec x) {
		double E = 0.0;
		for (V v : hds.getVertices()) {
			E += getVertexEnergy(v, x);
		}
		return E;
	}
	
	public double getVertexEnergy(V v, Vec aVec) {
		if (HalfEdgeUtils.isBoundaryVertex(v)) {
			return 0;
		}
		double sl = getLeftSinProduct(v, aVec);
		double sr = getRightSinProduct(v, aVec);
		return (sl - sr) * (sl - sr);
	}
	
	
	protected double getLeftSinProduct(V v, Vec aVec) {
		if (HalfEdgeUtils.isBoundaryVertex(v)) {
			return 0;
		}
		double sl = 1.0;
		for (E eIn : HalfEdgeUtils.incomingEdges(v)) {
			double alpha_ki = getAlpha(eIn, aVec);
			double alpha_ij = getAlpha(eIn.getNextEdge(), aVec);
			double alpha_jk = getAlpha(eIn.getPreviousEdge(), aVec);
			double betaLeft = calculateTriangleAngle(alpha_ij, alpha_jk, alpha_ki);
			sl *= sin(betaLeft);
		}
		return sl;
	}
	
	protected double getRightSinProduct(V v, Vec aVec) {
		if (HalfEdgeUtils.isBoundaryVertex(v)) {
			return 0;
		}
		double sr = 1.0;
		for (E eIn : HalfEdgeUtils.incomingEdges(v)) {
			double alpha_ki = getAlpha(eIn, aVec);
			double alpha_ij = getAlpha(eIn.getNextEdge(), aVec);
			double alpha_jk = getAlpha(eIn.getPreviousEdge(), aVec);
			double betaRight = calculateTriangleAngle(alpha_jk, alpha_ki, alpha_ij);
			sr *= sin(betaRight);
		}
		return sr;
	}
	
	
	protected double getOppsiteBeta(E e, Vec aVec) {
		double alpha_ki = getAlpha(e, aVec);
		double alpha_ij = getAlpha(e.getNextEdge(), aVec);
		double alpha_jk = getAlpha(e.getPreviousEdge(), aVec);
		return calculateTriangleAngle(alpha_ij, alpha_jk, alpha_ki);
	}
	
	
	protected double getAlpha(E e, Vec aVec) {
		int index = solverIndices.get(e.getIndex());
		return aVec.getValue(index);
	}
	
	
}
