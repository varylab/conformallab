package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.SparseUtility.makeNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.ConesUtility;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class EuclideanUnwrapper implements Unwrapper{

	
	public Vector unwrap(CoHDS surface, int numCones, boolean quantizeCones) throws Exception {
		surface.prepareInvariantDataEuclidean();
		
		// cones
		Collection<CoVertex> cones = null;
		if (numCones > 0) {
			cones = ConesUtility.setUpMesh(surface, numCones);
		}
		
		CEuclideanOptimizable opt = new CEuclideanOptimizable(surface);
		int n = opt.getDomainDimension();
		
		// optimization
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n,makeNonZeros(surface));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(1E-8);
		optimizer.setMaxIterations(150);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		
		if (quantizeCones && numCones > 0) {
			cones = ConesUtility.quantizeCones(surface, cones);
			n = opt.getDomainDimension();
			u = new DenseVector(n);
			H = new CompRowMatrix(n,n,makeNonZeros(surface));
			optimizer.setHessianTemplate(H);
			try {
				optimizer.minimize(u, opt);
			} catch (NotConvergentException e) {
				throw new UnwrapException("Cone quantization did not succeed: " + e.getMessage());
			}
		}
		// layout
		if (numCones > 0) {
			ConesUtility.cutMesh(surface, cones, u);
		}
		return u;
	}

	
}
