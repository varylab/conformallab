package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.heds.util.SparseUtility.getPETScNonZeros;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import no.uib.cipr.matrix.DenseVector;

import org.eclipse.core.runtime.IProgressMonitor;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.util.CuttingUtility;
import de.varylab.discreteconformal.heds.util.HomologyUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.PETSc;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;

public class CHyperbolicUnwrapperPETSc implements CUnwrapper{

	
	public void unwrap(CoHDS hds, IProgressMonitor mon) throws UnwrapException {
		mon.beginTask("Unwrapping", 3);
		mon.subTask("Cut to disk");
		int X = hds.numVertices() - hds.numEdges() / 2 + hds.numFaces();
		if (X >= 3) {
			List<Set<CoEdge>> paths = HomologyUtility.getGeneratorPaths(hds.getVertex(0));
			Set<CoEdge> masterPath = new HashSet<CoEdge>();
			for (Set<CoEdge> path : paths) {
				masterPath.addAll(path);
			}
			for (CoEdge e : masterPath) {
				CuttingUtility.cutAtEdge(e);
			}
		}
		
		mon.subTask("Minimizing");
		hds.prepareInvariantDataHyperbolic();
		
		Tao.Initialize();
		CEuclideanApplication app = new CEuclideanApplication(hds);
		int n = app.getDomainDimension();
		
		Vec u = new Vec(n);
		Mat H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(hds));
		H.assemble();
		
		app.setInitialSolutionVec(u);
		app.setHessianMat(H, H);	
		
		Tao optimizer = new Tao(Tao.Method.CG);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(1E-5, 0, 0);
		optimizer.solve();
		mon.worked(1);
		
		double [] uValues = u.getArray();
		// layout
		mon.subTask("Layout");
		CHyperbolicLayout.doLayout(hds, new DenseVector(uValues)); 
		mon.worked(1);
		mon.done();
	}

}
