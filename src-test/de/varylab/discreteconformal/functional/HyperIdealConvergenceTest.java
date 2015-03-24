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
		PETSc.optionsSetValue("-tao_lmm_vectors", "20");
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
		for (CoVertex v : hds.getVertices()) {
			double sum = 0.0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				sum += e.getPreviousEdge().getBeta();
				log.info("beta@" + e + ": " + (e.getPreviousEdge().getBeta() / PI));
			}
			Assert.assertEquals(2*PI, sum, tolerance);
		}
		for (CoEdge e : hds.getPositiveEdges()) {
			double sum = e.getAlpha() + e.getOppositeEdge().getAlpha();
			log.info("alphaSum@" + e + ": " + (sum / PI));
			Assert.assertEquals(e.getTheta(), sum, tolerance);
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
			u.setValue(i, 1.0 + 0.01*Math.abs(rnd.nextDouble()), INSERT_VALUES);
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
			1.3169578891731424, 1.3169578927823975, 1.3169579011731716, 1.3169579020905602, 
			2.2924316583333115, 2.2924316392049833, 2.292431678989581, 2.292431643089526, 
			2.2924316941659173, 2.292431674031334, 2.2924316940930196, 2.2924316825837465, 
			2.292431682108128, 2.292431635825035, 2.292431663416895, 2.2924316882259688, 
			1.1959012101177978E-8, -2.0563476436963373E-8, -4.456976728393993E-8, -1.8429237951612702E-8, 
			-6.6205808102341355E-9, 1.5717234379207197E-8, -1.7293648067715866E-8, -4.483256556006885E-8, 
			-9.573468060377547E-10, -3.1141227023802106E-8, 1.0908622136947316E-8, 4.776674403327454E-8, 
			7.312334494445497E-9, 3.7181879656685486E-8, 2.14898248315642E-8, 2.8691886984531554E-9, 
			6.47831816268519E-9, -3.331867297470249E-9, 1.2954234834539E-8, 3.051358129457958E-8, 
			4.344545817244695E-9, 1.2223017352538969E-8, -1.6138776412832104E-8, -2.7556701966816714E-8
		};
		Assert.assertArrayEquals(expectedSolution, uVec, 1E-6);
		u.restoreArray();
		app.evaluateObjectiveAndGradient(u, null);
		checkSolution(hds);
	}
	
}
