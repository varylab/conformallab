package de.varylab.discreteconformal.unwrapper;

import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CSphericalOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;

public class SphericalUnwrapper implements Unwrapper {

	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	private CoVertex
		layoutRoot = null;
	
	@Override
	public void unwrap(CoHDS hds, int g, AdapterSet a) throws Exception {
		double maxLength = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = a.get(Length.class, e, Double.class);
			maxLength = maxLength < l ? l : maxLength;
		}
		double scale = 1.0/2.0;
		
		CSphericalOptimizable opt = new CSphericalOptimizable(hds);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, a, zeroU, scale);

		// optimization
		int n = opt.getDomainDimension();
		DenseVector u = new DenseVector(n);
		Matrix H = opt.getHessianTemplate();
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
//		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.BiCGstab); 
		optimizer.setError(gradTolerance);
		optimizer.setMaxIterations(maxIterations);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}

		layoutRoot = hds.getVertex(0);
		SphericalLayout.doLayout(hds, layoutRoot, opt.getFunctional(), u);
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
