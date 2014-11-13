package de.varylab.discreteconformal.functional;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import de.jreality.math.P3;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
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
	
	@Before
	public void create() {
		rnd.setSeed(0);
		for (int i = 0; i < 100; i++) {
			CoVertex v = hds.addNewVertex();
			v.P[0] = rnd.nextGaussian();
			v.P[1] = rnd.nextGaussian();
			v.P[2] = rnd.nextGaussian();
			v.P[3] = 1.0;
			Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
		}
	}
	
	
	@Override
	public void init() {
		MyDomainValue x = new MyDomainValue(new DenseVector(4));
		x.set(0, rnd.nextDouble() - 0.5);
		x.set(1, rnd.nextDouble() - 0.5);
		x.set(2, rnd.nextDouble() - 0.5);
		x.set(3, 1.0);
		
		setHDS(hds);
		setFunctional(fun);
		setXGradient(x);
		setXHessian(x);
	}
	
	
	@Test
	public void testConvergence() throws Exception {
		NewtonOptimizer min = new NewtonOptimizer();
		Vector x = new DenseVector(new double[] {0,0,0,1});
		
		double[] cm = getCenterOfMass(hds);
		double normCm = Pn.norm(cm, Pn.EUCLIDEAN);
		Assert.assertTrue("center of mass is not the origin before normalization", normCm > 0.05);
		
		Optimizable opt = fun.getOptimizatble(hds); 
		DenseVector g = new DenseVector(4);
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
		double[] xp = {x.get(0), x.get(1), x.get(2), x.get(3)};
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
