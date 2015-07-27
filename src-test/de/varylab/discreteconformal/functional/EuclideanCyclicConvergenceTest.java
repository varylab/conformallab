package de.varylab.discreteconformal.functional;

import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

import org.junit.Assert;
import org.junit.Test;

import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.functional.hds.MyConformalAdapters.CPhi;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomEdgeInfo;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.TestUtility;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class EuclideanCyclicConvergenceTest  {

	private CTheta
		theta = new CTheta();
	private CPhi
		phi = new CPhi();
	private CVariable
		variable = new CVariable();
	private CLambda
		lambda = new CLambda();
	private CAlpha
		alpha = new CAlpha();
	private CInitialEnergy
		energy = new CInitialEnergy();
	public EuclideanCyclicFunctional<CoVertex, CoEdge, CoFace>
		functional = new EuclideanCyclicFunctional<CoVertex, CoEdge, CoFace>(variable, theta, phi, lambda, alpha, energy);
	
	@Test
	public void testEuclideanConvergence() {
		CoHDS hds = TestUtility.readOBJ(EuclideanCyclicConvergenceTest.class, "cathead.obj"); 
		
//		 one edge is circular
		CoEdge circularEdge = null;
		for (CoEdge e : hds.getPositiveEdges()) {
			CoVertex s = e.getStartVertex();
			CoVertex t = e.getTargetVertex();
			if (isBoundaryVertex(s) || isBoundaryVertex(t)) {
				continue;
			}
			circularEdge = e;
			CustomEdgeInfo info = new CustomEdgeInfo();
			info.circularHoleEdge = true;
			info.phi = Math.PI - 0.1; // with modified angle sum phi
			e.info = info;
			e.getOppositeEdge().info = info;
			break;
		}
		
		// optimization
		CEuclideanOptimizable opt = new CEuclideanOptimizable(hds);
		AdapterSet a = new ConformalAdapterSet();
		int n = UnwrapUtility.prepareInvariantDataEuclidean(opt.getFunctional(), hds, a);
		DenseVector u = new DenseVector(n);
		// set variable lambda start values
		for (CoEdge e : hds.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				u.set(e.getSolverIndex(), e.getLambda());
			}
		}
		
		Matrix H = new CompRowMatrix(n, n, functional.getNonZeroPattern(hds));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.BiCGstab);
		optimizer.setError(1E-11);
		optimizer.setMaxIterations(10);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			Assert.fail(e.getLocalizedMessage());
		}
		Assert.assertEquals(Math.PI - 0.1, circularEdge.getAlpha() + circularEdge.getOppositeEdge().getAlpha(), 1E-12);
	}
	
	
}
