package de.varylab.discreteconformal.functional;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static java.lang.Math.PI;

import java.util.Random;
import java.util.logging.Logger;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperIdealApplication;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class HyperIdealConvergenceTest {

	private Logger
		log = Logger.getLogger(HyperIdealConvergenceTest.class.getName());
	private double 
		tolerance = 1E-7;
	
	@BeforeClass
	public static void initPetsc() {
		NativePathUtility.set("native");
		Tao.Initialize();
		LoggingUtility.initLogging();
		PETSc.optionsSetValue("-tao_lmm_vectors", "5");
		PETSc.optionsSetValue("-tao_lmm_scale_type", "broyden");
		PETSc.optionsSetValue("-tao_lmm_broyden_phi", "0.125");
		PETSc.optionsSetValue("-tao_lmm_rescale_type", "scalar");
		PETSc.optionsSetValue("-tao_lmm_rescale_history", "5");
		PETSc.optionsSetValue("-tao_lmm_rescale_alpha", "5.0");
		PETSc.optionsSetValue("-tao_lmm_rescale_beta", "0.5");
		PETSc.optionsSetValue("-tao_lmm_limit_type", "relative");
		PETSc.optionsSetValue("-tao_lmm_limit_mu", "1.0");
		PETSc.optionsSetValue("-tao_lmm_limit_nu", "1.0");
	}
	
	private void checkSolution(CoHDS hds) throws Exception {
		for (CoVertex v :hds.getVertices()) {
			double sum = 0.0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				sum += e.getPreviousEdge().getBeta();
			}
			Assert.assertEquals(2*PI, sum, tolerance);
		}
	}
	
	@Test
	public void testHyperIdealConvergence() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiled();
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		int n = app.getDomainDimension();
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vec u = new Vec(n);
		for (int i = 0; i < n; i++) {
			u.setValue(i, 1E-12, INSERT_VALUES);//1.1 + 0.1*Math.abs(rnd.nextDouble()), INSERT_VALUES);
		}
		app.setInitialSolutionVec(u);
		Vec lowerBounds = new Vec(n);
		Vec upperBounds = new Vec(n);
		lowerBounds.set(1E-12);
		upperBounds.set(Double.MAX_VALUE);
		app.setVariableBounds(lowerBounds, upperBounds);
		log.info("start   : " + u.toString());
		
		Tao optimizer = new Tao(Tao.Method.BLMVM);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(tolerance, 0, 0); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(50);
		optimizer.solve();
		Assert.assertEquals(ConvergenceFlags.CONVERGED_ATOL, optimizer.getSolutionStatus().reason);
		UnwrapUtility.logSolutionStatus(optimizer, log);
		double[] uVec = u.getArray();
		double[] expectedSolution = {
			1.1462158341786262, 1.1462158341786262, 1.1462158341786262, 1.1462158341786262, 
			1.7627471737467797, 1.7627471737467797, 1.7627471737467866, 1.7627471737467866, 
			1.7627471737467797, 1.7627471737467797, 1.7627471737467866, 1.7627471737467866, 
			1.7627471737467797, 1.7627471737467866, 1.7627471737467797, 1.7627471737467866, 
			2.633915794495759, 2.633915794495759, 2.633915794495759, 2.633915794495759, 
			2.633915794495759, 2.633915794495759
		};
		Assert.assertArrayEquals(expectedSolution, uVec, 1E-6);
		u.restoreArray();
		log.info("solution: " + u.toString());
		app.evaluateObjectiveAndGradient(u, null);
		checkSolution(hds);
	}
	
	
	@Test
	public void testHyperIdealConvergenceWithBranchPoints() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiledWithBranchPoints();
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		app.setFromOptions();
		int n = app.getDomainDimension();
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vec u = new Vec(n);
		for (int i = 0; i < n; i++) {
			u.setValue(i, 0.1 + 0.01*Math.abs(rnd.nextDouble()), INSERT_VALUES);
		}
		app.setInitialSolutionVec(u);
		Vec lowerBounds = new Vec(n);
		Vec upperBounds = new Vec(n);
		lowerBounds.set(-Double.MAX_VALUE);
		for (int i = 0; i < 16; i++) {
			lowerBounds.setValue(i, 1E-12, INSERT_VALUES);
		}
		upperBounds.set(Double.MAX_VALUE);
		app.setVariableBounds(lowerBounds, upperBounds);
		log.info("start   : " + u.toString());
		
		Tao optimizer = new Tao(Tao.Method.BLMVM);
		optimizer.setFromOptions();
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(tolerance, 0, 0); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(50);
		optimizer.solve();
		log.info(optimizer.getSolutionStatus().toString());
		Assert.assertEquals(ConvergenceFlags.CONVERGED_ATOL, optimizer.getSolutionStatus().reason);
		UnwrapUtility.logSolutionStatus(optimizer, log);
		double[] uVec = u.getArray();
		double[] expectedSolution = {
			1.3287091283689778, 1.328709144060893, 1.1828942665778752, 1.182894274711472, 
			2.1314267110516862, 2.1314266886793214, 2.145393699070382, 2.2924316728851686, 
			2.2362567280260284, 2.2362567650314595, 2.2924316732803294, 2.1453936986802735, 
			2.1314267097029105, 2.2672397657673327, 2.1314266893826184, 1.0416102421811257, 
			-0.4502780224399055, -0.45027802846560555, -0.0411278384808546, -0.0411278214456359, 
			-0.08996021141665035, -0.08996020110009696, -0.030869897541421235, -0.03086991142974196, 
			0.010918417011581107, 0.01091843405536372, -0.04817186915701951, -0.04817188236012954, 
			0.010918434042524154, 0.010918417025070244, -0.04817188210816582, -0.04817186915473929, 
			-0.08996020096657777, -0.0899602107928726, -0.03086991097650844, -0.030869902607411796, 
			-0.450278027759061, -0.450278021732768, -0.04112782131503768, -0.04112783768419421
		};
		Assert.assertArrayEquals(expectedSolution, uVec, 1E-6);
		u.restoreArray();
		app.evaluateObjectiveAndGradient(u, null);
		checkSolution(hds);
		
		for (CoEdge e : HalfEdgeUtils.incomingEdges(hds.getVertex(2))) {
			log.info("Beta: " + e.getPreviousEdge().getBeta());
		}
	}
	
}
