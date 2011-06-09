package de.varylab.discreteconformal.functional;

import java.util.Random;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

import org.junit.Before;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.ConformalEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.unwrapper.numerics.MTJGradient;
import de.varylab.discreteconformal.unwrapper.numerics.MTJHessian;
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
	public void ceate() {
		for (int i = 0; i < 100; i++) {
			CoVertex v = hds.addNewVertex();
			v.P[0] = rnd.nextGaussian();
			v.P[1] = rnd.nextGaussian();
			v.P[2] = rnd.nextGaussian() * 5; // unbalance
			v.P[3] = 1.0;
			Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
			System.out.println(Pn.norm(v.P, Pn.EUCLIDEAN));
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
		Optimizable opt = new Optimizable() {
			
			@Override
			public Matrix getHessianTemplate() {
				return new DenseMatrix(4, 4);
			}
			
			@Override
			public Integer getDomainDimension() {
				return 4;
			}
			
			@Override
			public Double evaluate(Vector x, Vector gradient, Matrix hessian) {
				MTJDomain u = new MTJDomain(x);
				MTJGradient G = new MTJGradient(gradient);
				MTJHessian H = new MTJHessian(hessian);
				ConformalEnergy E = new ConformalEnergy();
				fun.evaluate(hds, u, E, G, H);
				return E.get();
			}

			@Override
			public Double evaluate(Vector x, Vector gradient) {
				MTJDomain u = new MTJDomain(x);
				MTJGradient G = new MTJGradient(gradient);
				ConformalEnergy E = new ConformalEnergy();
				fun.evaluate(hds, u, E, G, null);
				return E.get();
			}

			@Override
			public Double evaluate(Vector x, Matrix hessian) {
				MTJDomain u = new MTJDomain(x);
				MTJHessian H = new MTJHessian(hessian);
				ConformalEnergy E = new ConformalEnergy();
				fun.evaluate(hds, u, E, null, H);
				return E.get();
			}

			@Override
			public Double evaluate(Vector x) {
				MTJDomain u = new MTJDomain(x);
				ConformalEnergy E = new ConformalEnergy();
				fun.evaluate(hds, u, E, null, null);
				return E.get();
			}
		}; 
		min.minimize(x, opt);
		System.out.println("x " + x);
	}
	
	
	
}
