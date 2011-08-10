package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.SparseUtility.getPETScNonZeros;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicApplication;
import de.varylab.discreteconformal.util.UnwrapUtility;

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
			Mat H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(surface));
			H.assemble();
			app.setHessianMat(H, H);
		}
		
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
