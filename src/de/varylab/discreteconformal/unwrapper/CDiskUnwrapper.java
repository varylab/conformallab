package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.heds.util.SparseUtility.makeNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

import org.eclipse.core.runtime.IProgressMonitor;

import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.util.ConesUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class CDiskUnwrapper implements CUnwrapper{

	
	private int
		numCones = 0;
	private boolean
		quantizeCones = true;
	
	public CDiskUnwrapper(int numCones, boolean quantizeCones) {
		this.numCones = numCones;
		this.quantizeCones = quantizeCones;
	}
	
	
	public void unwrap(CoHDS hds, IProgressMonitor mon) throws UnwrapException {
		mon.beginTask("Unwrapping", 2 + (quantizeCones ? 2 : 0));
		hds.prepareInvariantDataEuclidean();
		
		// cones
		Collection<CoVertex> cones = null;
		if (numCones > 0) {
			mon.subTask("Processing " + numCones + " cones");
			cones = ConesUtility.setUpMesh(hds, numCones);
			mon.worked(1);
		}
		
		CEuclideanOptimizable opt = new CEuclideanOptimizable(hds);
		int n = opt.getDomainDimension();
		
		// optimization
		mon.subTask("Minimizing");
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n,makeNonZeros(hds));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(1E-5);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			mon.setCanceled(true);
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		mon.worked(1);
		
		if (quantizeCones && numCones > 0) {
			mon.subTask("Quantizing Cone Singularities");
			cones = ConesUtility.quantizeCones(hds, cones);
			n = opt.getDomainDimension();
			u = new DenseVector(n);
			H = new CompRowMatrix(n,n,makeNonZeros(hds));
			optimizer.setHessianTemplate(H);
			try {
				optimizer.minimize(u, opt);
			} catch (NotConvergentException e) {
				mon.setCanceled(true);
				throw new UnwrapException("Cone quantization did not succeed: " + e.getMessage());
			}
			mon.worked(1);
		}
		
		// layout
		mon.subTask("Layout");
		if (numCones > 0) {
			ConesUtility.cutMesh(hds, cones, u);
		}
		CEuclideanLayout.doLayout(hds, u);
		mon.worked(1);
		mon.done();
	}

	
	
	
	public int getNumCones() {
		return numCones;
	}



	public void setNumCones(int numCones) {
		this.numCones = numCones;
	}



	public boolean isQuantizeCones() {
		return quantizeCones;
	}



	public void setQuantizeCones(boolean quantizeCones) {
		this.quantizeCones = quantizeCones;
	}
	
}
