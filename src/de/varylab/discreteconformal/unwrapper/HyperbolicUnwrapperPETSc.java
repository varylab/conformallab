package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.SparseUtility.getPETScNonZeros;
import static de.varylab.jpetsc.PETSc.PETSC_DEFAULT;
import static de.varylab.jtao.ConvergenceFlags.CONVERGED_ATOL;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicApplication;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;

public class HyperbolicUnwrapperPETSc implements Unwrapper{

	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	
	
	@Override
	public Vector unwrap(CoHDS surface, AdapterSet aSet) throws Exception {
		UnwrapUtility.prepareInvariantDataHyperbolic(surface, aSet);
		
		Tao.Initialize();
		CHyperbolicApplication app = new CHyperbolicApplication(surface);
		int n = app.getDomainDimension(); 
		
		Vec u = new Vec(n);
		Mat H = Mat.createSeqAIJ(n, n, PETSC_DEFAULT, getPETScNonZeros(surface));
		H.assemble();
		
		app.setInitialSolutionVec(u);
		app.setHessianMat(H, H);	
		
		Tao optimizer = new Tao(Tao.Method.NLS);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(gradTolerance, gradTolerance, gradTolerance); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(maxIterations);
		System.out.println("Using grad tolerance " + gradTolerance);
		optimizer.solve();
		if (optimizer.getSolutionStatus().reason != CONVERGED_ATOL) {
			throw new RuntimeException("Optinizer did not converge: \n" + optimizer.getSolutionStatus());
		}
		System.out.println(optimizer.getSolutionStatus());
		DenseVector result = new DenseVector(u.getArray());
		u.restoreArray();
		return result;
	}

	
	@Override
	public void setGradientTolerance(double tol) {
		gradTolerance = tol;
	}
	
	@Override
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}
	
}
