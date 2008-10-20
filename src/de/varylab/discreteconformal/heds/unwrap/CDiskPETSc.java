package de.varylab.discreteconformal.heds.unwrap;

import static de.varylab.discreteconformal.heds.util.SparseUtility.makeNonZeros;

import java.util.Collection;

import no.uib.cipr.matrix.DenseVector;

import org.eclipse.core.runtime.IProgressMonitor;

import de.varylab.discreteconformal.heds.CCones;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.heds.PETScCHDSEvaluator;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.PETSc;
import de.varylab.jpetsc.PETScException;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;

public class CDiskPETSc implements CUnwrapper{

	
	private int
		numCones = 0;
	private boolean
		quantizeCones = true;
	
	public CDiskPETSc(int numCones, boolean quantizeCones) {
		this.numCones = numCones;
		this.quantizeCones = quantizeCones;
	}
	
	
	public int[] getPETScNonZeros(CHDS hds){
		int n = hds.getDomainDimension();
		int [] nnz = new int[n];
		{
			int [][] sparseStucture = makeNonZeros(hds);
			for(int i = 0; i < n; i++){
				nnz[i] = sparseStucture[i].length;
			}
		}
		return nnz;
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
		Vec u;
		Mat H;
		Tao optimizer;
		try {
			Tao.Initialize();
		} catch (PETScException e1) {
			e1.printStackTrace();
			throw(new UnwrapException("Tao could not be loaded:"+e1.getMessage()));
		}
		PETScCHDSEvaluator app = new PETScCHDSEvaluator(hds);
		
		try{
			u = new Vec(n);
			H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(hds));
			H.assemblyBegin(Mat.AssemblyType.FINAL_ASSEMBLY);
			H.assemblyEnd(Mat.AssemblyType.FINAL_ASSEMBLY);


			optimizer = new Tao("tao_nls");
			app.setInitialSolutionVec(u);
			app.setHessianMat(H, H);			
			optimizer.setApplication(app);
	//		TODO: optimizer.setStepController(new ArmijoStepController());
	//		optimizer.setTolerances(1E-5, 1E-5, 0, 0);
			optimizer.setGradientTolerances(1E-5, 0, 0);
			optimizer.solve();
//			GetSolutionStatusResult res = optimizer.getSolutionStatus();
//			System.out.println(res);
		}
		catch(PETScException e){
			mon.setCanceled(true);
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		mon.worked(1);

		double [] uValues; 
		
		if (quantizeCones && numCones > 0) {
			mon.subTask("Quantizing Cone Singularities");
			try {
				app = new PETScCHDSEvaluator(hds);
				uValues = u.getArray();
				cones = CCones.quantizeCones(hds, cones, new DenseVector(uValues), app.calculateAlphas(u));
				u.restoreArray();
				n = hds.getDomainDimension();
				u = new Vec(n);
				u.assemblyBegin();
				u.assemblyEnd();
				app.setInitialSolutionVec(u);
				H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(hds));
				H.assemblyBegin(Mat.AssemblyType.FINAL_ASSEMBLY);
				H.assemblyEnd(Mat.AssemblyType.FINAL_ASSEMBLY);
				app.setHessianMat(H, H);
				optimizer.setApplication(app);
				optimizer.solve();
			} catch (PETScException e) {
				mon.setCanceled(true);
				throw new UnwrapException("Cone quantization did not succeed: " + e.getMessage());
			}
			mon.worked(1);
		}
		
		// layout
		mon.subTask("Layout");
		try {
			uValues = u.getArray();
			if (numCones > 0) {
				CCones.cutMesh(hds, cones, new DenseVector(uValues));
			}
			CDiskLayout.doLayout(hds,  new DenseVector(uValues), app.calculateAlphas(u));
			u.restoreArray();
		} catch (PETScException e) {
			e.printStackTrace();
		}
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
