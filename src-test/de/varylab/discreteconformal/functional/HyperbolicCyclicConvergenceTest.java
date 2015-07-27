package de.varylab.discreteconformal.functional;

import static java.lang.Math.PI;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

import org.junit.Assert;
import org.junit.Test;

import de.jtem.halfedge.util.HalfEdgeUtils;
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
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicOptimizable;
import de.varylab.discreteconformal.util.TestUtility;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class HyperbolicCyclicConvergenceTest  {

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
		CoHDS hds = TestUtility.readOBJ(HyperbolicCyclicConvergenceTest.class, "hyperbolic_convergence_model.obj"); 
		
//		one triangle of edges is circular
		for (CoFace f : hds.getFaces()) {
			if (!HalfEdgeUtils.isInteriorFace(f)) continue;
			CoEdge e1 = f.getBoundaryEdge();
			CoEdge e2 = e1.getNextEdge();
			CoEdge e3 = e2.getNextEdge();
			CustomEdgeInfo info = new CustomEdgeInfo();
			info.circularHoleEdge = true;
			e1.info = info;
			e2.info = info;
			e3.info = info;
			e1.getOppositeEdge().info = info;
			e2.getOppositeEdge().info = info;
			e3.getOppositeEdge().info = info;
			break;
		}
		
		// optimization
		CHyperbolicOptimizable opt = new CHyperbolicOptimizable(hds);
		AdapterSet a = new ConformalAdapterSet();
		ZeroU zeroU = new ZeroU();
		int n = UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, a, zeroU);
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
		// check flatness
		for (CoVertex v : hds.getVertices()) {
			double alpha = 0.0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				alpha += e.getPreviousEdge().getAlpha();
			}
			Assert.assertEquals(2*PI, alpha, 1E-8);
		}
	}
	
	
}
