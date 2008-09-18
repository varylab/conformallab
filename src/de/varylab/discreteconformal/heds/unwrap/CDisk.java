package de.varylab.discreteconformal.heds.unwrap;

import static de.varylab.discreteconformal.heds.util.SparseUtility.makeNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

import org.eclipse.core.runtime.IProgressMonitor;

import de.varylab.discreteconformal.heds.CCones;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.math.optimization.NotConvergentException;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer.Solver;
import de.varylab.discreteconformal.math.optimization.stepcontrol.ArmijoStepController;

public class CDisk implements CUnwrapper{

	
	private int
		numCones = 0;
	private boolean
		quantizeCones = true;
	
	public CDisk(int numCones, boolean quantizeCones) {
		this.numCones = numCones;
		this.quantizeCones = quantizeCones;
	}
	
	
	public void unwrap(CHDS hds, IProgressMonitor mon) throws UnwrapException {
		mon.beginTask("Unwrapping", 2 + (quantizeCones ? 2 : 0));
		hds.prepareInvariantData();
		
		// cones
		Collection<CVertex> cones = null;
		if (numCones > 0) {
			mon.subTask("Processing " + numCones + " cones");
			cones = CCones.setUpMesh(hds, numCones);
			mon.worked(1);
		}
		int n = hds.getDomainDimension();
		
		// optimization
		mon.subTask("Minimizing");
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n,makeNonZeros(hds));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(1E-5);
		try {
			optimizer.minimize(u, hds);
		} catch (NotConvergentException e) {
			mon.setCanceled(true);
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		mon.worked(1);
		
		if (quantizeCones && numCones > 0) {
			mon.subTask("Quantizing Cone Singularities");
			cones = CCones.quantizeCones(hds, cones, u, hds.calculateAlphas(u));
			n = hds.getDomainDimension();
			u = new DenseVector(n);
			H = new CompRowMatrix(n,n,makeNonZeros(hds));
			optimizer.setHessianTemplate(H);
			try {
				optimizer.minimize(u, hds);
			} catch (NotConvergentException e) {
				mon.setCanceled(true);
				throw new UnwrapException("Cone quantization did not succeed: " + e.getMessage());
			}
			mon.worked(1);
		}
		
		// layout
		mon.subTask("Layout");
		if (numCones > 0) {
			CCones.cutMesh(hds, cones, u);
		}
		CDiskLayout.doLayout(hds, u, hds.calculateAlphas(u));
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