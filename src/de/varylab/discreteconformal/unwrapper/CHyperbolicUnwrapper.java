package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.heds.util.SparseUtility.makeNonZeros;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

import org.eclipse.core.runtime.IProgressMonitor;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.util.CuttingUtility;
import de.varylab.discreteconformal.heds.util.HomologyUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicOptimizable;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class CHyperbolicUnwrapper implements CUnwrapper{

	
	public void unwrap(CoHDS hds, IProgressMonitor mon) throws UnwrapException {
		mon.beginTask("Unwrapping", 3);
		//TODO maybe we dont need to cut
//		int X = hds.numVertices() - hds.numEdges() / 2 + hds.numFaces();
//		int g = (2 - X) / 2;
//		System.err.println("genus of the surface is " + g);
//		if (g >= 2) {
//			System.err.println("Cut to disk...");
//			mon.subTask("Cut to disk");
//			List<Set<CoEdge>> paths = HomologyUtility.getGeneratorPaths(hds.getVertex(0));
//			System.err.println("Found " + paths.size() + " paths");
//			Set<CoEdge> masterPath = new HashSet<CoEdge>();
//			for (Set<CoEdge> path : paths) {
//				masterPath.addAll(path);
//			}
//			Set<CoEdge> oppEdges = new HashSet<CoEdge>();
//			for (CoEdge e : masterPath) {
//				if (masterPath.contains(e.getOppositeEdge())) {
//					oppEdges.add(e.getOppositeEdge());
//				}
//			}
//			masterPath.removeAll(oppEdges);
//			for (CoEdge e : masterPath) {
//				CuttingUtility.cutAtEdge(e);
//			}
//			X = hds.numVertices() - hds.numEdges() / 2 + hds.numFaces();
//			g = (2 - X) / 2;
//			System.err.println("genus of the surface after cutting is " + g);
//		}
		
		mon.subTask("Minimizing");
		hds.prepareInvariantDataHyperbolic();
		
		CHyperbolicOptimizable opt = new CHyperbolicOptimizable(hds);
		int n = opt.getDomainDimension();
		
		// optimization
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n,makeNonZeros(hds));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(1E-5);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			mon.setCanceled(true);
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		mon.worked(1);
		
		System.err.println("U: " + u);
		
		// layout
		mon.subTask("Layout");
		CHyperbolicLayout.doLayout(hds, u); 
		mon.worked(1);
		mon.done();
	}

}
