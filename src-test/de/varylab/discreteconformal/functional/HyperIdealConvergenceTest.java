package de.varylab.discreteconformal.functional;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static de.varylab.discreteconformal.functional.HyperIdealFunctionalTest.createLawsonsSurface;

import java.util.Random;
import java.util.logging.Logger;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.jtem.jtao.TaoVec;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperIdealApplication;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class HyperIdealConvergenceTest {

	private Logger
		log = Logger.getLogger(HyperIdealConvergenceTest.class.getName());
	private double[]
		expectedSolution = {1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 1.7627471360523435, 1.7627471360523428, 1.7627471360523435, 1.7627471360523428, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531};
	
	@BeforeClass
	public static void initPetsc() {
		NativePathUtility.set("native");
		Tao.Initialize();
	}
	
	@Test
	public void testHyperIdealConvergence() {
		LoggingUtility.initLogging();
		CoHDS hds = createLawsonsSurface();
		
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		int n = app.getDomainDimension();
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vec u = new Vec(n);
		for (int i = 0; i < n; i++) {
//			u.setValue(i, 0.5 + abs(rnd.nextDouble()), INSERT_VALUES);
			u.setValue(i, 1.0 + 0.1*Math.abs(rnd.nextDouble()), INSERT_VALUES);
		}
		app.setInitialSolutionVec(u);
		TaoVec lowerBounds = new TaoVec(n);
		TaoVec upperBounds = new TaoVec(n);
		lowerBounds.setToConstant(0.0);
		upperBounds.setToConstant(Double.MAX_VALUE);
		log.info("start   : " + u.toString());
		
		Tao optimizer = new Tao(Tao.Method.LMVM);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(1E-9, 1E-9, 1E-9); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(300);
		//optimizer.setVariableBounds(lowerBounds, upperBounds);
		optimizer.solve();
		Assert.assertEquals(ConvergenceFlags.CONVERGED_ATOL, optimizer.getSolutionStatus().reason);
		UnwrapUtility.logSolutionStatus(optimizer, log);
		double[] uVec = u.getArray();
		Assert.assertArrayEquals(expectedSolution, uVec, 1E-6);
		u.restoreArray();
		log.info("solution: " + u.toString());
	}
	
	
}
