package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.heds.util.SparseUtility.getPETScNonZeros;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicApplication;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.PETSc;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;
import de.varylab.jtao.Tao.GetSolutionStatusResult;

public class CHyperbolicUnwrapperPETSc implements CUnwrapper{

	
	public void unwrap(CoHDS hds) throws UnwrapException {
		Vector u = getConformalFactors(hds);
		CHyperbolicLayout.doLayout(hds, u);
	}

	public Vector getConformalFactors(CoHDS hds) {
		hds.prepareInvariantDataHyperbolic();
		
		Tao.Initialize();
		CHyperbolicApplication app = new CHyperbolicApplication(hds);
		int n = app.getDomainDimension();
		
		Vec u = new Vec(n);
		Mat H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(hds));
		H.assemble();
		
		app.setInitialSolutionVec(u);
		app.setHessianMat(H, H);	
		
		Tao optimizer = new Tao(Tao.Method.LMVM);
		optimizer.setApplication(app);
		optimizer.setTolerances(1E-10, 0, 0, 0); 
		 
		optimizer.solve();
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		System.out.println("Minimization: " + status);
		return new DenseVector(u.getArray());
	}
	
	
}
