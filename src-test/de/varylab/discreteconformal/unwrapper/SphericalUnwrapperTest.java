package de.varylab.discreteconformal.unwrapper;

import java.util.Random;

import junit.framework.Assert;

import org.junit.Test;

import de.jreality.math.Rn;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CSphericalApplication;
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
	public void testUnwrap() throws Exception {
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
	
}
