package de.varylab.discreteconformal.heds.unwrap;

import static de.varylab.discreteconformal.heds.util.SparseUtility.getPETScNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;

import org.eclipse.core.runtime.IProgressMonitor;

import de.varylab.discreteconformal.heds.CConesUtility;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.math.CEuclideanApplication;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.PETSc;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;
import de.varylab.jtao.Tao.GetSolutionStatusResult;

public class CDiskPETSc implements CUnwrapper{

	private int
		numCones = 0;
	private boolean
		quantizeCones = true;
	
	public CDiskPETSc(int numCones, boolean quantizeCones) {
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
			cones = CConesUtility.setUpMesh(hds, numCones);
			mon.worked(1);
		}

		// optimization
		mon.subTask("Minimizing");
		Vec u;
		Mat H;
		Tao optimizer;
		Tao.Initialize();
		CEuclideanApplication app = new CEuclideanApplication(hds);
		int n = app.getDomainDimension();
		
		u = new Vec(n);
		H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(hds));
		H.assemble();
		
		app.setInitialSolutionVec(u);
		app.setHessianMat(H, H);	
		
		optimizer = new Tao(Tao.Method.NLS);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(1E-5, 0, 0);
		optimizer.solve();
		
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		System.out.println("Minimization: " + status);
		if (status.reason.cvalue() < 0) {
			mon.setCanceled(true);
			throw new UnwrapException("Optimization did not succeed: " + status);
		}
		mon.worked(1);

		if (quantizeCones && numCones > 0) {
			mon.subTask("Quantizing Cone Singularities");
			// calculating cones
			cones = CConesUtility.quantizeCones(hds, cones);
			
			// optimizing conformal structure
			CEuclideanApplication app2 = new CEuclideanApplication(hds);
			n = app2.getDomainDimension();

			u = new Vec(n);
			H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(hds));
			H.assemble();
			
			app2.setInitialSolutionVec(u);
			app2.setHessianMat(H, H);
			
			optimizer.setApplication(app2);
			optimizer.setGradientTolerances(1E-5, 0, 0);
			optimizer.solve();
			
			status = optimizer.getSolutionStatus();
			System.out.println("Cone Quantization: " + status);
			if (status.reason.cvalue() < 0) {
				mon.setCanceled(true);
				throw new UnwrapException("Cone quantization did not succeed: " + status);
			}
			mon.worked(1);
		}
		
		// layout
		mon.subTask("Layout");
		double [] uValues = u.getArray();
		if (numCones > 0) {
			CConesUtility.cutMesh(hds, cones, new DenseVector(uValues));
		}
		CDiskLayout.doLayout(hds,  new DenseVector(uValues));
		u.restoreArray();
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
