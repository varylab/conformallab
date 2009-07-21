package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.SparseUtility.makeNonZeros;
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

public class CHyperbolicUnwrapper implements Unwrapper{

	
	public Vector unwrap(CoHDS surface) throws Exception {
		surface.prepareInvariantDataHyperbolic();
		CHyperbolicOptimizable opt = new CHyperbolicOptimizable(surface);
		int n = opt.getDomainDimension();
		
		// optimization
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n,makeNonZeros(surface));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.BiCGstab); 
		optimizer.setError(1E-8);
		optimizer.setMaxIterations(150);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
//			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		return u;
	}
	

}
