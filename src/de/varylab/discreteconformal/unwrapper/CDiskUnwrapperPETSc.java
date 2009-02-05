package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.heds.util.SparseUtility.getPETScNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.util.ConesUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.PETSc;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;
import de.varylab.jtao.Tao.GetSolutionStatusResult;

public class CDiskUnwrapperPETSc implements CUnwrapper{

	private int
		numCones = 0;
	private boolean
		quantizeCones = true;
	
	public CDiskUnwrapperPETSc(int numCones, boolean quantizeCones) {
		this.numCones = numCones;
		this.quantizeCones = quantizeCones;
	}
	
	
	public void unwrap(CoHDS hds) throws UnwrapException {
		// cones
		Collection<CoVertex> cones = null;
		if (numCones > 0) {
			cones = ConesUtility.setUpMesh(hds, numCones);
		}

		// optimization
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
		
		optimizer = new Tao(Tao.Method.LMVM);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(1E-5, 0, 0);
		optimizer.solve();
		
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		System.out.println("Minimization: " + status);
		if (status.reason.cvalue() < 0) {
			throw new UnwrapException("Optimization did not succeed: " + status);
		}

		if (quantizeCones && numCones > 0) {
			// calculating cones
			cones = ConesUtility.quantizeCones(hds, cones);
			
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
				throw new UnwrapException("Cone quantization did not succeed: " + status);
			}
		}
		
		// layout
		double [] uValues = u.getArray();
		if (numCones > 0) {
			ConesUtility.cutMesh(hds, cones, new DenseVector(uValues));
		}
		CEuclideanLayout.doLayout(hds,  new DenseVector(uValues));
		u.restoreArray();
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
