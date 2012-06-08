package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility.calculateBeta;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.tan;

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
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;

public class SinConditionFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	HDS extends HalfEdgeDataStructure<V, E, F>
> extends TaoApplication implements TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {

	private HDS
		hds = null;
	private Map<Integer, Integer>
		solverIndices = null;
	private boolean recCall = false;
	
	public SinConditionFunctional(HDS hds, Map<Integer, Integer> undirectedIndexMap) {
		this.hds = hds;
		this.solverIndices = undirectedIndexMap;
	}
	
	public void calculateAndSetInitialSolution(AdapterSet a) {
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
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		// objective
		double E = 0.0;
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				continue;
			}
			double sl = getLeftLogSinSum(v, x);
			double sr = getRightLogSinSum(v, x);
			double e = sr - sl;
			E += e*e;
			if (g != null && !recCall) {
				for (E ein : incomingEdges(v)) {
					double bl = getOppsiteBeta(ein.getNextEdge(), x);
					double br = getOppsiteBeta(ein.getOppositeEdge().getPreviousEdge(), x);
					double br2 = getOppsiteBeta(ein, x);
					int iin = solverIndices.get(ein.getIndex());
					int iopp = solverIndices.get(ein.getPreviousEdge().getIndex());
					g.add(iin, 2*e*(-1/tan(br) - 1/tan(bl)));
					g.add(iopp, 2*e*(1/tan(br2) + 1/tan(bl)));
				}
			}
		}
		if (g != null && !recCall) {
			System.out.println("check --------------------------------");
			System.out.println("x: " + x);
			System.out.println("h: " + g);
			recCall = true;
			defaultComputeGradient(x, g);
			recCall = false;
			System.out.println("f: " + g);
		}
		return E;
	}
	
	
	@Override
	public PreconditionerType evaluateHessian(Vec x, Mat H, Mat Hpre) {
		return PreconditionerType.SAME_NONZERO_PATTERN;
	}

	
	public double getVertexProductEnergy(V v, Vec aVec) {
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
			double betaLeft = calculateBeta(alpha_ij, alpha_jk, alpha_ki);
			sl *= sin(betaLeft);
		}
		return sl;
	}
	
	protected double getLeftLogSinSum(V v, Vec aVec) {
		if (HalfEdgeUtils.isBoundaryVertex(v)) {
			return 0;
		}
		double sl = 0.0;
		for (E eIn : HalfEdgeUtils.incomingEdges(v)) {
			double alpha_ki = getAlpha(eIn, aVec);
			double alpha_ij = getAlpha(eIn.getNextEdge(), aVec);
			double alpha_jk = getAlpha(eIn.getPreviousEdge(), aVec);
			double betaLeft = calculateBeta(alpha_ij, alpha_jk, alpha_ki);
			sl += log(sin(betaLeft));
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
			double betaRight = calculateBeta(alpha_jk, alpha_ki, alpha_ij);
			sr *= sin(betaRight);
		}
		return sr;
	}
	
	protected double getRightLogSinSum(V v, Vec aVec) {
		if (HalfEdgeUtils.isBoundaryVertex(v)) {
			return 0;
		}
		double sr = 0.0;
		for (E eIn : HalfEdgeUtils.incomingEdges(v)) {
			double alpha_ki = getAlpha(eIn, aVec);
			double alpha_ij = getAlpha(eIn.getNextEdge(), aVec);
			double alpha_jk = getAlpha(eIn.getPreviousEdge(), aVec);
			double betaRight = calculateBeta(alpha_jk, alpha_ki, alpha_ij);
			sr += log(sin(betaRight));
		}
		return sr;
	}
	
	
	protected double getOppsiteBeta(E e, Vec aVec) {
		double alpha_ki = getAlpha(e, aVec);
		double alpha_ij = getAlpha(e.getNextEdge(), aVec);
		double alpha_jk = getAlpha(e.getPreviousEdge(), aVec);
		return calculateBeta(alpha_ij, alpha_jk, alpha_ki);
	}

	
	protected double getAlpha(E e, Vec aVec) {
		int index = solverIndices.get(e.getIndex());
		return aVec.getValue(index);
	}
	
	
}
