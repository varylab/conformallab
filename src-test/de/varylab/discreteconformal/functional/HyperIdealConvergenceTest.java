package de.varylab.discreteconformal.functional;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;

import java.util.Random;
import java.util.logging.Logger;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoHDS;
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
	}
	
	@Test
	public void testHyperIdealConvergence() {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiled();
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		int n = app.getDomainDimension();
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vec u = new Vec(n);
		for (int i = 0; i < n; i++) {
			u.setValue(i, 1.0 + 0.1*Math.abs(rnd.nextDouble()), INSERT_VALUES);
		}
		app.setInitialSolutionVec(u);
		Vec lowerBounds = new Vec(n);
		Vec upperBounds = new Vec(n);
		lowerBounds.set(-Double.MAX_VALUE);
		lowerBounds.setValue(0, 1E-12, INSERT_VALUES);
		lowerBounds.setValue(1, 1E-12, INSERT_VALUES);
		lowerBounds.setValue(2, 1E-12, INSERT_VALUES);
		lowerBounds.setValue(3, 1E-12, INSERT_VALUES);
		upperBounds.set(Double.MAX_VALUE);
		app.setVariableBounds(lowerBounds, upperBounds);
		log.info("start   : " + u.toString());
		
		Tao optimizer = new Tao(Tao.Method.BLMVM);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(tolerance, 0, 0); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(300);
		optimizer.solve();
		Assert.assertEquals(ConvergenceFlags.CONVERGED_ATOL, optimizer.getSolutionStatus().reason);
		UnwrapUtility.logSolutionStatus(optimizer, log);
		double[] uVec = u.getArray();
		double[] expectedSolution = {
			1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 
			1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 
			1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 
			1.7627471360523435, 1.7627471360523428, 1.7627471360523435, 1.7627471360523428, 
			2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531, 
			2.633915759978531, 2.633915759978531
		};
		Assert.assertArrayEquals(expectedSolution, uVec, 1E-6);
		u.restoreArray();
		log.info("solution: " + u.toString());
	}
	
	
	@Test
	public void testHyperIdealConvergenceWithBranchPoints() {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiledWithBranchPoints();
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		int n = app.getDomainDimension();
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vec u = new Vec(n);
		for (int i = 0; i < n; i++) {
			u.setValue(i, 10.0 + 0.1*Math.abs(rnd.nextDouble()), INSERT_VALUES);
		}
		app.setInitialSolutionVec(u);
		Vec lowerBounds = new Vec(n);
		Vec upperBounds = new Vec(n);
		lowerBounds.set(-Double.MAX_VALUE);
		lowerBounds.setValue(0, 1E-10, INSERT_VALUES);
		lowerBounds.setValue(1, 1E-10, INSERT_VALUES);
		lowerBounds.setValue(2, 1E-10, INSERT_VALUES);
		lowerBounds.setValue(3, 1E-10, INSERT_VALUES);
		upperBounds.set(Double.MAX_VALUE);
		app.setVariableBounds(lowerBounds, upperBounds);
		log.info("start   : " + u.toString());
		
		Tao optimizer = new Tao(Tao.Method.BLMVM);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(tolerance, 0, 0); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(300);
		optimizer.solve();
		log.info(optimizer.getSolutionStatus().toString());
		Assert.assertEquals(ConvergenceFlags.CONVERGED_ATOL, optimizer.getSolutionStatus().reason);
		UnwrapUtility.logSolutionStatus(optimizer, log);
		double[] uVec = u.getArray();
		double[] expectedSolution = new double[28];
		Assert.assertArrayEquals(expectedSolution, uVec, 1E-6);
		u.restoreArray();
		log.info("solution: " + u.toString());
	}
	
}
