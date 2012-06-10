package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility.calculateBeta;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.tan;

import java.util.Map;

import sun.security.x509.AVA;

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
import de.jtem.jtao.TaoAppAddGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoAppAddObjective;
import de.jtem.jtao.TaoApplication;

public class SinConditionFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	HDS extends HalfEdgeDataStructure<V, E, F>
> extends TaoApplication implements TaoAppAddObjective, TaoAppAddGrad, TaoAppAddHess {

	private HDS
		hds = null;
	private Map<Integer, Integer>
		solverIndices = null;
	
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
	public double evaluateObjective(Vec x) {
		// objective
		double E = 0.0;
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				continue;
			}
			double sl = getLeftLogSinSum(v, x);
			double sr = getRightLogSinSum(v, x);
//			double sl = getLeftSinProduct(v, x);
//			double sr = getRightSinProduct(v, x);
			double e = sr - sl;
			E += e*e;
		}
		return E;
	}
	
	@Override
	public void evaluateGradient(Vec x, Vec g) {
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				continue;
			}
			double sl = getLeftLogSinSum(v, x);
			double sr = getRightLogSinSum(v, x);
			double e = sr - sl;
			if (g != null) {
				for (E ein : HalfEdgeUtils.incomingEdges(v)) {
					double bl = getOppsiteBeta(ein.getNextEdge(), x);
					double dblp = getOppBetaDerivativeWrtPrev(ein.getNextEdge(), x);
					double dbln = getOppBetaDerivativeWrtNext(ein.getNextEdge(), x);
					double br = getOppsiteBeta(ein.getOppositeEdge().getPreviousEdge(), x);
					double dbrn = getOppBetaDerivativeWrtNext(ein.getOppositeEdge().getPreviousEdge(), x);
					double bl2 = getOppsiteBeta(ein, x);
					double dbl2p = getOppBetaDerivativeWrtPrev(ein, x);
					int iin = solverIndices.get(ein.getIndex());
					int iopp = solverIndices.get(ein.getPreviousEdge().getIndex());
					g.add(iin, 2*e*(dbrn/tan(br) - dblp/tan(bl)));
					g.add(iopp, 2*e*(dbl2p/tan(bl2) - dbln/tan(bl)));
				}
			}
		}
	}
	
	@Override
	public PreconditionerType evaluateHessian(Vec x, Mat H, Mat Hpre) {
		return PreconditionerType.SAME_NONZERO_PATTERN;
	}

	
	protected double getLeftSinProduct(V v, Vec aVec) {
		if (HalfEdgeUtils.isBoundaryVertex(v)) {
			return 0;
		}
		double sl = 1.0;
		for (E eIn : HalfEdgeUtils.incomingEdges(v)) {
			double betaLeft = getOppsiteBeta(eIn, aVec);
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
			double betaLeft = getOppsiteBeta(eIn, aVec);
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
			double betaRight = getOppsiteBeta(eIn.getNextEdge(), aVec);
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
			double betaRight = getOppsiteBeta(eIn.getNextEdge(), aVec);
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
	
	protected double getOppBetaDerivativeWrtNext(E e, Vec aVec) {
		double alpha_ki = getAlpha(e, aVec);
		double alpha_ij = getAlpha(e.getNextEdge(), aVec);
		double alpha_jk = getAlpha(e.getPreviousEdge(), aVec);
		alpha_ki = IsothermicUtility.normalizeAngle(alpha_ki);
		alpha_ij = IsothermicUtility.normalizeAngle(alpha_ij);
		alpha_jk = IsothermicUtility.normalizeAngle(alpha_jk);
		double betaSign = Math.signum(alpha_jk - alpha_ij);
		if ((alpha_ki > alpha_jk && alpha_ki > alpha_ij) || (alpha_ki < alpha_jk && alpha_ki < alpha_ij)) {
			return -1*betaSign;
		} else {
			return 1*betaSign;
		}
	}
	
	protected double getOppBetaDerivativeWrtPrev(E e, Vec aVec) {
		double alpha_ki = getAlpha(e, aVec);
		double alpha_ij = getAlpha(e.getNextEdge(), aVec);
		double alpha_jk = getAlpha(e.getPreviousEdge(), aVec);
		alpha_ki = IsothermicUtility.normalizeAngle(alpha_ki);
		alpha_ij = IsothermicUtility.normalizeAngle(alpha_ij);
		alpha_jk = IsothermicUtility.normalizeAngle(alpha_jk);
		double betaSign = Math.signum(alpha_jk - alpha_ij);
		if ((alpha_ki > alpha_jk && alpha_ki > alpha_ij) || (alpha_ki < alpha_jk && alpha_ki < alpha_ij)) {
			return 1*betaSign;
		} else {
			return -1*betaSign;
		}
	}

	
	protected double getAlpha(E e, Vec aVec) {
		int index = solverIndices.get(e.getIndex());
		return aVec.getValue(index);
	}
	
	
}
