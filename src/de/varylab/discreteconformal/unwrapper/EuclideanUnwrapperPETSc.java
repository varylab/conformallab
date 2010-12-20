package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.SparseUtility.getPETScNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.util.UnwrapUtility.QuantizationMode;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.PETSc;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;
import de.varylab.jtao.Tao.GetSolutionStatusResult;

public class EuclideanUnwrapperPETSc implements Unwrapper {

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

	public static double
		lastGNorm = 0;
	
	@Override
	public Vector unwrap(CoHDS surface, AdapterSet aSet) throws Exception {
		UnwrapUtility.prepareInvariantDataEuclidean(surface, boundaryMode, boundaryQuantMode, aSet);
		// cones
		Collection<CoVertex> cones = ConesUtility.setUpCones(surface, numCones); 
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
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setGradientTolerances(gradTolerance, gradTolerance, gradTolerance);
		optimizer.setMaximumIterates(maxIterations);
		optimizer.solve();
		
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		lastGNorm = status.gnorm;
		if (status.reason.cvalue() < 0) {
			throw new UnwrapException("Optimization did not succeed: " + status);
		}

		if (!cones.isEmpty()) {
			if (conesMode != QuantizationMode.AllAngles) {
				// calculating cones
				cones = ConesUtility.quantizeCones(surface, cones, conesMode);
				
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
			double [] uValues = u.getArray();
			ConesUtility.cutMesh(surface, cones, new DenseVector(uValues));
		}
		
		// layout
		double [] uValues = u.getArray();
		DenseVector result = new DenseVector(uValues);
		u.restoreArray();
		return result; 
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
