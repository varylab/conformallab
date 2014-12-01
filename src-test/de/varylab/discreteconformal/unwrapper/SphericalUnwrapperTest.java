package de.varylab.discreteconformal.unwrapper;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import de.jreality.math.Rn;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.jpetsc.NormType;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CSphericalApplication;
import de.varylab.discreteconformal.unwrapper.numerics.CSphericalOptimizable;
import de.varylab.discreteconformal.util.TestUtility;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;

public class SphericalUnwrapperTest {

	private static Random
		rnd = new Random();
	
	static {
		NativePathUtility.set("native");
		Tao.Initialize();
	}
	
	@Test
	public void testPETScAtCriticalPoint() throws Exception {
		rnd.setSeed(2);
		
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		hds.addNewVertices(12);
		for (CoVertex v : hds.getVertices()) {
			double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			Rn.setEuclideanNorm(pos, 1.0, pos);
			aSet.set(Position.class, v, pos);
		}
		ConvexHull.convexHull(hds, aSet);
		
		CSphericalApplication opt = new CSphericalApplication(hds);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, aSet, zeroU, 0.5);
		
		SphericalUnwrapperPETSc uw = new SphericalUnwrapperPETSc();
		double[] u = uw.calculateConformalFactors(hds, aSet, opt);
		
		for (double ui : u) {
			Assert.assertEquals(0, ui, 1E-8);
		}
	}
	
	
	@Test
	public void testPETScAtGenericPoint() throws Exception {
		rnd.setSeed(2);
		
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		hds.addNewVertices(12);
		for (CoVertex v : hds.getVertices()) {
			double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			Rn.setEuclideanNorm(pos, 0.5 + rnd.nextDouble(), pos);
			aSet.set(Position.class, v, pos);
		}
		ConvexHull.convexHull(hds, aSet);
		
		CSphericalApplication opt = new CSphericalApplication(hds);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, aSet, zeroU, 0.5);
		
		SphericalUnwrapperPETSc uw = new SphericalUnwrapperPETSc();
		uw.setMaxIterations(200);
		uw.setGradientTolerance(1E-13);
		double[] u = uw.calculateConformalFactors(hds, aSet, opt);
		
		Vec uVec = new Vec(u.length);
		for (int i = 0; i < u.length; i++) {
			uVec.setValue(i, u[i], INSERT_VALUES);
		}
		Vec G = new Vec(u.length);
		opt.evaluateObjectiveAndGradient(uVec, G);
		double gNorm = G.norm(NormType.NORM_FROBENIUS);
		System.out.println("|G|: " + gNorm);
		Assert.assertEquals(0, gNorm, 1E-6);
	}
	
	
	@Test@Ignore("I still cannot get this to run with Petsc :-(")
	public void testPETScAtCathead() throws Exception {
		AdapterSet aSet = new ConformalAdapterSet();
		CoHDS hds = TestUtility.readOBJ(SphericalUnwrapperTest.class, "cathead_sphere.obj");
		
		CSphericalApplication opt = new CSphericalApplication(hds);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, aSet, zeroU, 0.5);
		
		SphericalUnwrapperPETSc uw = new SphericalUnwrapperPETSc();
		uw.setMaxIterations(1000);
		uw.setGradientTolerance(1E-12);
		double[] u = uw.calculateConformalFactors(hds, aSet, opt);
		
		Vec uVec = new Vec(u.length);
		for (int i = 0; i < u.length; i++) {
			uVec.setValue(i, u[i], INSERT_VALUES);
		}
		Vec G = new Vec(u.length);
		opt.evaluateObjectiveAndGradient(uVec, G);
		
		double gNorm = G.norm(NormType.NORM_FROBENIUS);
		System.out.println("|G|: " + gNorm);
		Assert.assertEquals(0, gNorm, 1E-6);
	}
	
	
	@Test
	public void testMTJAtGenericPoint() throws Exception {
		rnd.setSeed(2);
		
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		hds.addNewVertices(12);
		for (CoVertex v : hds.getVertices()) {
			double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			Rn.setEuclideanNorm(pos, 0.5 + rnd.nextDouble(), pos);
			aSet.set(Position.class, v, pos);
		}
		ConvexHull.convexHull(hds, aSet);
		
		CSphericalOptimizable opt = new CSphericalOptimizable(hds);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, aSet, zeroU, 0.5);
		
		SphericalUnwrapper uw = new SphericalUnwrapper();
		uw.setMaxIterations(200);
		uw.setGradientTolerance(1E-10);
		Vector u = uw.calculateConformalFactors(opt);
		
		Vector G = new DenseVector(u.size());
		opt.evaluate(u, G);
		
		double gNorm = G.norm(Norm.Two);
		System.out.println("|G|: " + gNorm);
		Assert.assertEquals(0, gNorm, 1E-6);
	}
	
	
	@Test
	public void testMTJAtCathead() throws Exception {
		AdapterSet aSet = new ConformalAdapterSet();
		CoHDS hds = TestUtility.readOBJ(SphericalUnwrapperTest.class, "cathead_sphere.obj");
		
		CSphericalOptimizable opt = new CSphericalOptimizable(hds);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, aSet, zeroU, 0.5);
		
		SphericalUnwrapper uw = new SphericalUnwrapper();
		uw.setMaxIterations(200);
		uw.setGradientTolerance(1E-7);
		Vector u = uw.calculateConformalFactors(opt);
		
		Vector G = new DenseVector(u.size());
		opt.evaluate(u, G);
		
		double gNorm = G.norm(Norm.Two);
		System.out.println("|G|: " + gNorm);
		Assert.assertEquals(0, gNorm, 1E-6);
	}

	
}
