package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility.calculateBeta;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.tan;

import java.util.HashMap;
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
import de.jtem.jpetsc.MatStructure;
import de.jtem.jpetsc.SNES;
import de.jtem.jpetsc.Vec;
import de.jtem.jpetsc.SNES.FunctionEvaluator;
import de.jtem.jpetsc.SNES.JacobianEvaluator;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.Method;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;

public class SinConditionApplication <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	HDS extends HalfEdgeDataStructure<V, E, F>
> extends TaoApplication implements TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess, FunctionEvaluator, JacobianEvaluator {

	protected HDS
		hds = null;
	protected Map<E, Integer>
		solverEdgeIndices = new HashMap<E, Integer>();
	protected Map<V, Integer>
		solverVertexIndices = new HashMap<V, Integer>();
	protected Map<E, Double>
		initialAlphas = new HashMap<E, Double>();
	protected int dim = -1;
	
	public SinConditionApplication(HDS hds) {
		this.hds = hds;
	}
	
	public void initialize(AdapterSet a, boolean excludeBoundary) {
		if (!a.isAvailable(CurvatureFieldMin.class, hds.getEdge(0).getClass(), double[].class)) {
			throw new IllegalArgumentException("No curvature directions found!");
		}
		for (E e : hds.getEdges()) {
			double[] N = a.getD(Normal.class, e);
			double[] Kmin = a.getD(CurvatureFieldMin.class, e);
			double[] E = a.getD(EdgeVector.class, e);
			double ae = IsothermicUtility.getSignedAngle(N, Kmin, E);
			initialAlphas.put(e, ae);
			initialAlphas.put(e.getOppositeEdge(), ae);
		}
		initialize(initialAlphas, excludeBoundary);
	}
	
	public void initialize(Map<E, Double> initAlphas, boolean excludeBoundary) {
		this.solverEdgeIndices = IsothermicUtility.createSolverEdgeIndexMap(hds, excludeBoundary);
		this.solverVertexIndices = IsothermicUtility.createSolverVertexIndexMap(hds);
		this.initialAlphas = initAlphas;
		int maxIndex = -1;
		for (E e : hds.getEdges()) {
			int index = solverEdgeIndices.get(e);
			maxIndex = maxIndex < index ? index : maxIndex;
		}
		this.dim = maxIndex + 1;
		Vec alpha = new Vec(dim);
		for (E e : hds.getEdges()) {
			int index = solverEdgeIndices.get(e);
			if (index >= 0) {
				double ae = initialAlphas.get(e);
				alpha.setValue(index, ae, InsertMode.INSERT_VALUES);
			}
		}	
		alpha.assemble();
		setInitialSolutionVec(alpha);
//		Mat H = createHessianTemplate();
//		setHessianMat(H, H);
		System.out.println("SinFunctional number of variables: " + alpha.getSize());
	}
	
	
	public void solveCG(int maxIterations, double tol) {
		Tao tao = new Tao(Method.CG);
		tao.setFromOptions();
		tao.setApplication(this);
		tao.setMaximumIterates(1000);
		tao.setTolerances(tol, tol, tol, tol);
		tao.setGradientTolerances(tol, tol, tol);
		System.out.println("energy before optimization: " + evaluateObjective(getSolutionVec()));
		tao.solve();
		System.out.println(tao.getSolutionStatus());
		System.out.println("energy after optimization: " + evaluateObjective(getSolutionVec()));
	}
	
	public void solveSNES(int maxIterations, double tol) {
		SNES snes = SNES.create();
		snes.setFromOptions();
		Vec f = new Vec(dim);
		Mat J = new Mat(dim, dim);
		for (int i = 0; i < dim; i++) {
			J.setValue(i, i, 0.0, InsertMode.INSERT_VALUES);
		}
		J.assemble();
		snes.setFunction(this, f);
		snes.setJacobian(this, J, J);
		snes.setTolerances(tol, tol, tol, maxIterations, 10000000);
		snes.getKSP().setInitialGuessNonzero(false);
		System.out.println("residual: " + snes.getFunctionNorm());
		System.out.println("guess\n" + getSolutionVec());
		snes.solve(null, getSolutionVec());
		System.out.println("residual: " + snes.getFunctionNorm());
		System.out.println("solution\n" + getSolutionVec());
		System.out.println(snes.getConvergedReason());
	}
	
	
	@Override
	public MatStructure evaluateJacobian(Vec x, Mat J, Mat Jpre) {
		J.zeroEntries();
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			int vIndex = solverVertexIndices.get(v);
			for (E ein : HalfEdgeUtils.incomingEdges(v)) {
				double bl = getOppsiteBeta(ein.getNextEdge(), x);
				double dblp = getOppBetaDerivativeWrtPrev(ein.getNextEdge(), x);
				double dbln = getOppBetaDerivativeWrtNext(ein.getNextEdge(), x);
				double br = getOppsiteBeta(ein.getOppositeEdge().getPreviousEdge(), x);
				double dbrn = getOppBetaDerivativeWrtNext(ein.getOppositeEdge().getPreviousEdge(), x);
				double bl2 = getOppsiteBeta(ein, x);
				double dbl2p = getOppBetaDerivativeWrtPrev(ein, x);
				int iin = solverEdgeIndices.get(ein);
				int iopp = solverEdgeIndices.get(ein.getPreviousEdge());
				J.setValue(vIndex, iin, dbrn/tan(br) - dblp/tan(bl), InsertMode.ADD_VALUES);
				J.setValue(vIndex, iopp, dbl2p/tan(bl2) - dbln/tan(bl), InsertMode.ADD_VALUES);
			}
		}
		J.assemble();
		return MatStructure.DIFFERENT_NONZERO_PATTERN;
	}

	@Override
	public void evaluateFunction(Vec x, Vec f) {
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			int index = solverVertexIndices.get(v);
			double sl = getLeftLogSinSum(v, x);
			double sr = getRightLogSinSum(v, x);
			f.setValue(index, sl - sr, INSERT_VALUES);
		}
	}
	
	
	public DBFSolution<V, E, F, HDS> getDBFSolution() {
		DBFSolution<V, E, F, HDS> s = new DBFSolution<V, E, F, HDS>();
		s.hds = hds;
		s.initialAlphaMap = initialAlphas;
		Vec vec = getSolutionVec();
		s.solutionAlphaMap = new HashMap<E, Double>();
		for (E e : hds.getEdges()) {
			int index = solverEdgeIndices.get(e);
			if (index >= 0) {
				s.solutionAlphaMap.put(e, vec.getValue(index));
			} else {
				Double initAlpha = initialAlphas.get(e);
				s.solutionAlphaMap.put(e, initAlpha);
			}
		}
		return s;
	}
	
	public double evaluateObjective(Vec x) {
		return evaluateObjectiveAndGradient(x, null);
	}
	
	
	@Override
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		double E = 0.0;
		if (g != null) g.zeroEntries();
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			// energy
			double sl = getLeftLogSinSum(v, x);
			double sr = getRightLogSinSum(v, x);
			double e = sr - sl;
			E += e*e;
			// gradient
			if (g == null) continue;
			for (E ein : HalfEdgeUtils.incomingEdges(v)) {
				double bl = getOppsiteBeta(ein.getNextEdge(), x);
				double dblp = getOppBetaDerivativeWrtPrev(ein.getNextEdge(), x);
				double dbln = getOppBetaDerivativeWrtNext(ein.getNextEdge(), x);
				double br = getOppsiteBeta(ein.getOppositeEdge().getPreviousEdge(), x);
				double dbrn = getOppBetaDerivativeWrtNext(ein.getOppositeEdge().getPreviousEdge(), x);
				double bl2 = getOppsiteBeta(ein, x);
				double dbl2p = getOppBetaDerivativeWrtPrev(ein, x);
				int iin = solverEdgeIndices.get(ein);
				int iopp = solverEdgeIndices.get(ein.getPreviousEdge());
				g.add(iin, -2*e*(dbrn/tan(br) - dblp/tan(bl)));
				g.add(iopp, -2*e*(dbl2p/tan(bl2) - dbln/tan(bl)));
			}
		}
		return E;
	}
	
	
	@Override
	public PreconditionerType evaluateHessian(Vec x, Mat H, Mat Hpre) {
//		H.zeroEntries();
//		for (V v : hds.getVertices()) {
//			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
//			for (E ein : HalfEdgeUtils.incomingEdges(v)) {
//				int iin = solverIndices.get(ein);
//				int iopp = solverIndices.get(ein.getPreviousEdge());
//				
//			}
//		}
		return PreconditionerType.SAME_NONZERO_PATTERN;
	}
	
	
	public Mat createHessianTemplate() {
		int dim = hds.numEdges() / 2;
		int[] nz = new int[dim]; 
		for (E e : hds.getPositiveEdges()) {
			int i = solverEdgeIndices.get(e);
			if (HalfEdgeUtils.isBoundaryEdge(e)) {
				nz[i] = 3;
			} else {
				nz[i] = 5;
			}
		}
		Mat H = Mat.createSeqAIJ(dim, dim, -1, nz);
		H.zeroEntries();
		H.assemble();
		return H;
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
		int index = solverEdgeIndices.get(e);
		if (index >= 0) {
			return aVec.getValue(index);
		} else {
			return initialAlphas.get(e);
		}
	}
	
	public int getDimension() {
		return dim;
	}
	
}
