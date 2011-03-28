package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.SparseUtility.makeNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.util.UnwrapUtility.QuantizationMode;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class EuclideanUnwrapper implements Unwrapper{

	private QuantizationMode
		conesMode = QuantizationMode.AllAngles,
		boundaryQuantMode = QuantizationMode.AllAngles;
	private BoundaryMode
		boundaryMode = BoundaryMode.Isometric;
	private int 
		maxIterations = 150,
		numCones = 0;
	private double
		gradTolerance = 1E-8;
	
	
	@Override
	public Vector unwrap(CoHDS surface, AdapterSet aSet) throws Exception {
		UnwrapUtility.prepareInvariantDataEuclidean(surface, boundaryMode, boundaryQuantMode, aSet);
		// cones
		Collection<CoVertex> cones = ConesUtility.setUpCones(surface, numCones); 
		// optimization
		CEuclideanOptimizable opt = new CEuclideanOptimizable(surface);
		int n = opt.getDomainDimension();
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n,makeNonZeros(surface));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(gradTolerance);
		optimizer.setMaxIterations(maxIterations);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		if (!cones.isEmpty()) {
			if (conesMode != QuantizationMode.AllAngles) {
				cones = ConesUtility.quantizeCones(surface, cones, conesMode);
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
		}
		return u;
	}
	
	public void setNumCones(int numCones) {
		this.numCones = numCones;
	}

	public void setConeMode(QuantizationMode quantizationMode) {
		this.conesMode = quantizationMode;
	}

	@Override
	public void setGradientTolerance(double tol) {
		gradTolerance = tol;
	}
	
	@Override
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}

	public void setBoundaryQuantMode(QuantizationMode boundaryQuantMode) {
		this.boundaryQuantMode = boundaryQuantMode;
	}
	
	public void setBoundaryMode(BoundaryMode boundaryMode) {
		this.boundaryMode = boundaryMode;
	}
	
}
