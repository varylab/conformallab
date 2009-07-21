package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.SparseUtility.getPETScNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.tasks.Unwrap.QuantizationMode;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.discreteconformal.util.ConesUtility;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.PETSc;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;
import de.varylab.jtao.Tao.GetSolutionStatusResult;

public class EuclideanUnwrapperPETSc implements Unwrapper {

	private QuantizationMode
		quantizationMode = QuantizationMode.Quads;
	private int 
		numCones = 0;
	private boolean
		quantizeCones = false;
	
	
	public Vector unwrap(CoHDS surface) throws Exception {
		surface.prepareInvariantDataEuclidean();
		// cones
		Collection<CoVertex> cones = null;
		if (numCones > 0) {
			cones = ConesUtility.setUpMesh(surface, numCones);
		}

		// optimization
		Vec u;
		Mat H;
		Tao optimizer;
		Tao.Initialize();
		CEuclideanApplication app = new CEuclideanApplication(surface);
		int n = app.getDomainDimension();
		
		u = new Vec(n);
		H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(surface));
		H.assemble();
		
		app.setInitialSolutionVec(u);
		app.setHessianMat(H, H);	
		
		optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(1E-8, 0, 0);
		optimizer.setMaximumIterates(150);
		optimizer.solve();
		
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		System.out.println("Minimization: " + status);
		if (status.reason.cvalue() < 0) {
			throw new UnwrapException("Optimization did not succeed: " + status);
		}

		if (quantizeCones && numCones > 0) {
			// calculating cones
			cones = ConesUtility.quantizeCones(surface, cones, quantizationMode);
			
			// optimizing conformal structure
			CEuclideanApplication app2 = new CEuclideanApplication(surface);
			n = app2.getDomainDimension();

			u = new Vec(n);
			H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(surface));
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
			ConesUtility.cutMesh(surface, cones, new DenseVector(uValues));
		}
		DenseVector result = new DenseVector(uValues);
		u.restoreArray();
		return result; 
	}
	
	public void setNumCones(int numCones) {
		this.numCones = numCones;
	}
	
	public void setQuantizeCones(boolean quantizeCones) {
		this.quantizeCones = quantizeCones;
	}
	
	public void setQuantizationMode(QuantizationMode quantizationMode) {
		this.quantizationMode = quantizationMode;
	}

	
}