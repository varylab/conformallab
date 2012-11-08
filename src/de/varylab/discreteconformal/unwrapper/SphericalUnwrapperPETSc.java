package de.varylab.discreteconformal.unwrapper;

import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CSphericalApplication;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;

public class SphericalUnwrapperPETSc implements Unwrapper {

	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	private CoVertex
		layoutRoot = null;
	
	static {
		Tao.Initialize();		
	}
	
	@Override
	public void unwrap(CoHDS hds, int g, AdapterSet a) throws Exception {
		double maxLength = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = a.get(Length.class, e, Double.class);
			maxLength = maxLength < l ? l : maxLength;
		}
		double scale = 1.0/2.0;
		
		CSphericalApplication opt = new CSphericalApplication(hds);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, a, zeroU, scale);

		double[] u = calculateConformalFactors(hds, a, opt); 
		Vector uVec = new DenseVector(u);
		
		layoutRoot = hds.getVertex(0);
		SphericalLayout.doLayout(hds, layoutRoot, opt.getFunctional(), uVec);
	}
	
	double[] calculateConformalFactors(CoHDS surface, AdapterSet aSet, CSphericalApplication app) {

		int n = app.getDomainDimension(); 
		Vec u = new Vec(n);
		// set variable lambda start values
		boolean hasCircularEdges = false;
		for (CoEdge e : surface.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				u.setValue(e.getSolverIndex(), e.getLambda(), InsertMode.INSERT_VALUES);
				hasCircularEdges = true;
			}
		}
		app.setInitialSolutionVec(u);
		if (!hasCircularEdges) {
			Mat H = app.getHessianTemplate();
			app.setHessianMat(H, H);
		}
		
		Vec G = new Vec(n);
		app.computeGradient(u, G);
		System.out.println(G);
		
		Tao optimizer = new Tao(hasCircularEdges ? Tao.Method.LMVM : Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(gradTolerance, gradTolerance, gradTolerance); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(maxIterations);
		System.out.println("Using grad tolerance " + gradTolerance);
		optimizer.solve();
		if (optimizer.getSolutionStatus().reason != ConvergenceFlags.CONVERGED_ATOL) {
			throw new RuntimeException("Optinizer did not converge: \n" + optimizer.getSolutionStatus());
		}
		System.out.println(optimizer.getSolutionStatus());
		double[] uVec = u.getArray();
		u.restoreArray();
		return uVec;
	}
	
	

	@Override
	public void setGradientTolerance(double tol) {
		gradTolerance = tol;
	}
	@Override
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}

	@Override
	public void setCutRoot(CoVertex root) {
	}

	@Override
	public CuttingInfo<CoVertex, CoEdge, CoFace> getCutInfo() {
		return null;
	}

	@Override
	public Map<CoEdge, Double> getlengthMap() {
		return null;
	}

	@Override
	public CoVertex getLayoutRoot() {
		return layoutRoot;
	}
	
	
}
