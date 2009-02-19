package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.heds.util.SparseUtility.makeNonZeros;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicOptimizable;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class CHyperbolicUnwrapper implements CUnwrapper{

	
	public void unwrap(CoHDS hds) throws UnwrapException {
		Vector u = getConformalFactors(hds);
		CHyperbolicLayout.doLayout(hds, u);
	}
	
	
	public Vector getConformalFactors(CoHDS hds) throws UnwrapException {
		hds.prepareInvariantDataHyperbolic();
		return getConformalFactorsUnprepared(hds);
	}
	
	
	public Vector getConformalFactorsUnprepared(CoHDS hds) throws UnwrapException {
		CHyperbolicOptimizable opt = new CHyperbolicOptimizable(hds);
		int n = opt.getDomainDimension();
		
		// optimization
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n,makeNonZeros(hds));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.BiCGstab); 
		optimizer.setError(1E-7);
		optimizer.setMaxIterations(100);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		return u;
	}

}
