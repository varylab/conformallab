package de.varylab.discreteconformal.functional;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.math.P3;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.GetSolutionStatusResult;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.util.SparseUtility;
import de.varylab.mtjoptimization.Optimizable;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;

public class MobiusCenteringFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	private Random
		rnd = new Random();
	
	private AdapterSet 
		aSet = new ConformalAdapterSet();
	private MobiusCenteringFunctional<CoVertex, CoEdge, CoFace, Position4d>
		fun = new MobiusCenteringFunctional<CoVertex, CoEdge, CoFace, Position4d>(Position4d.class, aSet);
	private CoHDS
		hds = new CoHDS();
	
	@BeforeClass
	public static void staticInit() throws Exception {
		NativePathUtility.set("native");
		Tao.Initialize();
		LoggingUtility.initLogging();
	}
	
	@Before
	public void create() {
		rnd.setSeed(0);
//		for (int i = 0; i < 100; i++) {
//			CoVertex v = hds.addNewVertex();
//			v.P[0] = rnd.nextGaussian();
//			v.P[1] = rnd.nextGaussian();
//			v.P[2] = rnd.nextGaussian();
//			v.P[3] = 1.0;
//			Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
//		}
		hds.addNewVertex().P = new double[]{0.6666666666667581, 1.0, -0.2357022603955113, 1.0};
		hds.addNewVertex().P = new double[]{0.7179487179486832, 0.3846153846153386, 0.7252377242939296, 1.0};
		hds.addNewVertex().P = new double[]{-0.6666666666665687, 1.0, 0.2357022603955759, 1.0};
		hds.addNewVertex().P = new double[]{-0.10256410256407844, 0.38461538461540973, 1.0153328140114315, 1.0};
		hds.addNewVertex().P = new double[]{0.8000000000000451, -0.6000000000000179, -1.1313708498982586, 1.0};
		hds.addNewVertex().P = new double[]{0.7999999999999837, -0.5999999999999704, 0.5656854249492858, 1.0};
		hds.addNewVertex().P = new double[]{-1.3333333333333077, -0.5999999999995608, -0.37712361663256727, 1.0};
		hds.addNewVertex().P = new double[]{-0.2666666666667714, -0.6000000000000041, 0.94280904158208, 1.0};
	}
	
	
	@Override
	public void init() {
		MyDomainValue x = new MyDomainValue(new DenseVector(3));
		x.set(0, rnd.nextDouble() - 0.5);
		x.set(1, rnd.nextDouble() - 0.5);
		x.set(2, rnd.nextDouble() - 0.5);
		
		setHDS(hds);
		setFunctional(fun);
		setXGradient(x);
		setXHessian(x);
	}
	
	@Override
	public void testGradient() throws Exception {
		super.testGradient();
	}
	
	@Override
	public void testHessian() throws Exception {
		super.testHessian();
	}
	
	@Test
	public void testConvergence() throws Exception {
		NewtonOptimizer min = new NewtonOptimizer();
		Vector x = new DenseVector(new double[] {0,0,0});
		
		double[] cm = getCenterOfMass(hds);
		double normCm = Pn.norm(cm, Pn.EUCLIDEAN);
		Assert.assertTrue("center of mass is not the origin before normalization", normCm > 0.05);
		
		Optimizable opt = fun.getOptimizatble(hds); 
		DenseVector g = new DenseVector(3);
		opt.evaluate(x, g);
		System.out.println("x " + x);
		System.out.println(g.norm(Norm.Two));
		
		min.setMaxIterations(100);
		min.setError(1E-13);
		min.minimize(x, opt);
		
		opt.evaluate(x, g);
		System.out.println("x " + x);
		System.out.println(g.norm(Norm.Two));

		// transform via hyperbolic motion
		double[] xp = {x.get(0), x.get(1), x.get(2), 1.0};
		Pn.dehomogenize(xp, xp);
		double[] TInv = P3.makeTranslationMatrix(null, xp, Pn.HYPERBOLIC);
		double[] T = Rn.inverse(null, TInv);
		for (CoVertex v : hds.getVertices()) {
			Rn.matrixTimesVector(v.P, T, v.P);
		}
		
		cm = getCenterOfMass(hds);
		normCm = Pn.norm(cm, Pn.EUCLIDEAN);
		Assert.assertEquals("center of mass is in the origin", 0.0, normCm, 1E-12);
	}
	
	
	@Test
	public void testConvergencePETSc() throws Exception {
		double[] cm = getCenterOfMass(hds);
		double normCm = Pn.norm(cm, Pn.EUCLIDEAN);
		Assert.assertTrue("center of mass is not the origin before normalization", normCm > 0.05);
		
		TaoApplication app = fun.getTaoApplication(hds);
		Vec x = new Vec(fun.getDimension(hds));
		x.set(0.0);
		app.setInitialSolutionVec(x);
		Mat H = SparseUtility.getHessianTemplate(fun, hds);
		app.setHessianMat(H, H);
		
		Tao optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setGradientTolerances(1E-13, 0, 0);
		optimizer.setMaximumIterates(100);
		optimizer.solve();
		
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		System.out.println(status);
		Assert.assertTrue(status.reason.ordinal() <= 4);
		
		// transform via hyperbolic motion
		double[] xp = x.getArray(); x.restoreArray();
		xp = Pn.homogenize(null, xp);
		double[] TInv = P3.makeTranslationMatrix(null, xp, Pn.HYPERBOLIC);
		double[] T = Rn.inverse(null, TInv);
		for (CoVertex v : hds.getVertices()) {		
			Rn.matrixTimesVector(v.P, T, v.P);
		}
		
		cm = getCenterOfMass(hds);
		normCm = Pn.norm(cm, Pn.EUCLIDEAN);
		Assert.assertEquals("center of mass is in the origin", 0.0, normCm, 1E-15);
	}
	
	
	private double[] getCenterOfMass(CoHDS hds) {
		double[] cm = new double[4];
		double[] tmp = new double[4];
		for (CoVertex v : hds.getVertices()) {
			Rn.add(cm, Pn.dehomogenize(tmp, v.P), cm);
		}
		return Pn.dehomogenize(cm, cm);
	}
	
}
