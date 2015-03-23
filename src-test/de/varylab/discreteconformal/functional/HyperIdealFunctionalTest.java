package de.varylab.discreteconformal.functional;

import static java.lang.Math.abs;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CBeta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

public class HyperIdealFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	public static final Double
		eps = 1E-5,
		error = 1E-4;
	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CAlpha
		alpha = new CAlpha();
	private CBeta
		beta = new CBeta();
	private HyperIdealFunctional<CoVertex, CoEdge, CoFace>
		functional = new HyperIdealFunctional<>(variable, theta, alpha, beta);
	
	@BeforeClass
	public static void beforeClass() {
		LoggingUtility.initLogging();
	}
	
	@Override
	public void init() {
		setEps(eps);
		setError(error);
		setFunctional(functional);
		
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiled();
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, 0.5 + abs(rnd.nextDouble()));
		}
		MyDomainValue u = new MyDomainValue(x);
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
	}
	
	@Test@Ignore
	public void testHessian() throws Exception {
	}

	@Test
	public void testGradientWithHyperIdealAndIdealPoints() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiledWithBranchPoints();
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, 0.5 + abs(rnd.nextDouble()));
		}
		MyDomainValue u = new MyDomainValue(x);
		setHDS(hds);
		setXGradient(u);
		super.testGradient();
	}
	
	
	@Test
	public void testGradientInTheExtendedDomain() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiledWithBranchPoints();
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, 1.2 + abs(rnd.nextDouble()));
		}
		MyDomainValue u = new MyDomainValue(x);
		setHDS(hds);
		setXGradient(u);
		super.testGradient();
	}
	
}
