package de.varylab.discreteconformal.functional;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

import org.junit.Before;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedgetools.adapter.AdapterSet;
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
	private MobiusCenteringFunctional 
		fun = new MobiusCenteringFunctional(aSet);
	private CoHDS
		hds = new CoHDS();
	
	@Before
	public void create() {
		for (int i = 0; i < 100; i++) {
			CoVertex v = hds.addNewVertex();
			v.P[0] = rnd.nextGaussian();
			v.P[1] = rnd.nextGaussian();
			v.P[2] = rnd.nextGaussian() * 2; // unbalance
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
		
		Optimizable opt = fun.getOptimizatble(hds); 
		DenseVector g = new DenseVector(4);
		opt.evaluate(x, g);
		System.out.println("x " + x);
		System.out.println(g.norm(Norm.Two));
		
		min.setMaxIterations(100);
		min.setError(1E-13);
		min.minimize(x, opt);
		
		double xp[] = {x.get(0), x.get(1), x.get(2), x.get(3)};
		double scale = 1 / Math.sqrt(-Pn.innerProduct(xp, xp, Pn.HYPERBOLIC));
		Rn.times(xp, scale, xp);
		System.out.println("sqrt(-<x,x>): " + Pn.norm(xp, Pn.HYPERBOLIC));
		x = new DenseVector(xp);
		
		opt.evaluate(x, g);
		System.out.println("x " + x);
		System.out.println(g.norm(Norm.Two));
	}
	
}
