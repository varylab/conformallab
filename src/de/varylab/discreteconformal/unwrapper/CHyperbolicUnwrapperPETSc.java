package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.SparseUtility.getPETScNonZeros;
import static de.varylab.jpetsc.PETSc.PETSC_DEFAULT;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicApplication;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;

public class CHyperbolicUnwrapperPETSc implements Unwrapper{

	
	public Vector unwrap(CoHDS surface, int numCones, boolean quantizeCones) throws Exception {
		surface.prepareInvariantDataHyperbolic();
		Tao.Initialize();
		CHyperbolicApplication app = new CHyperbolicApplication(surface);
		int n = app.getDomainDimension();
		
		Vec u = new Vec(n);
		Mat H = Mat.createSeqAIJ(n, n, PETSC_DEFAULT, getPETScNonZeros(surface));
		H.assemble();
		
		app.setInitialSolutionVec(u);
		app.setHessianMat(H, H);	
		
		Tao optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(1E-8, 0, 0); 
		optimizer.setMaximumIterates(150);
		optimizer.solve();
		return new DenseVector(u.getArray());
	}

}
