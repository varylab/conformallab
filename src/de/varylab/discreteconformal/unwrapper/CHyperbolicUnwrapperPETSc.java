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
		
		Tao optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(1E-5, 0, 0);
		
		optimizer.solve();
		return new DenseVector(u.getArray());
	}
	
	
}
