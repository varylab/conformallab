package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;
import static java.lang.Math.asin;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;

import java.io.FileOutputStream;

import junit.framework.Assert;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import org.junit.Test;

import de.jreality.geometry.Primitives;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.writer.WriterOBJ;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.LengthTex;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition3d;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.jreality.ConverterHeds2JR;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.geometry.ComplexProjective1;
import de.jtem.mfc.group.Moebius;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.functional.SphericalFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

public class SphericalLayoutTest {

	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CLambda
		lambda = new CLambda();
	private CAlpha
		alpha = new CAlpha();
	private CInitialEnergy
		energy = new CInitialEnergy();
	private SphericalFunctional<CoVertex, CoEdge, CoFace>
		functional = new SphericalFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	
	@Test
	public void testSphericalLayout() {
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		FunctionalTest.createOctahedron(hds, aSet);
		
		for (CoEdge e : hds.getEdges()) {
			e.setAlpha(PI / 2);
			e.setLambda(2 * log(sqrt(2)/2));
		}
		
		Vector u = new DenseVector(6);
		SphericalLayout.doLayout(hds, hds.getVertex(0), functional, u);
		
		for (CoEdge e : hds.getPositiveEdges()) {
			double[] t1 = aSet.getD(TexturePosition3d.class, e.getStartVertex());
			double[] t2 = aSet.getD(TexturePosition3d.class, e.getTargetVertex());
			double l = Rn.euclideanDistance(t1, t2);
			Assert.assertEquals(sqrt(2), l, 1E-8);
		}
		
	}
	
	@Test
	public void testSphericalLayout2() throws Exception {
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		ConverterJR2Heds c2 = new ConverterJR2Heds();
		c2.ifs2heds(Primitives.icosahedron(), hds, aSet);
		for (CoVertex v : hds.getVertices()) {
			Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
		}
		
		CoEdge e0 = hds.getEdge(1);
		TopologyAlgorithms.flipEdge(e0);

		double lEuclidean = 1.0 / sin(2*PI / 5);
		double aSpherical = 6.0 * PI / 15.0;
		double lSpherical = 2 * asin(lEuclidean / 2);
		double lSphericalFlip = 2 * asin(sin(lSpherical) * sin(aSpherical));
		double lEuclideanFlip = 2 * sin(lSphericalFlip / 2);
		
		for (CoEdge e : hds.getEdges()) {
			if (e == e0 || e.getOppositeEdge() == e0) {
				// data for the flipped edge
				double l = aSet.get(Length.class, e, Double.class);
				Assert.assertEquals(l, lEuclideanFlip, 1E-6);
				e.setLambda(2 * log(lEuclideanFlip / 2));
				e.setAlpha(2 * aSpherical);
				continue;
			}
			double l = aSet.get(Length.class, e, Double.class);
			Assert.assertEquals(lEuclidean, l, 1E-6);
			e.setLambda(2 * log(lEuclidean / 2));
			if (e == e0.getNextEdge() || e == e0.getPreviousEdge() || 
				e == e0.getOppositeEdge().getNextEdge() || e == e0.getOppositeEdge().getPreviousEdge()) {
				e.setAlpha(aSpherical / 2);
			} else {
				e.setAlpha(aSpherical);
			}
		}
		
		Vector u = new DenseVector(6);
		SphericalLayout.doLayout(hds, hds.getVertex(0), functional, u);
		
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = aSet.get(LengthTex.class, e, Double.class);
			if (e == e0 || e.getOppositeEdge() == e0) {
				Assert.assertEquals(lEuclideanFlip, l, 1E-8);
			} else {
				Assert.assertEquals(lEuclidean, l, 1E-8);
			}
		}
		
		ConverterHeds2JR c = new ConverterHeds2JR();
		IndexedFaceSet ifs = c.heds2ifs(hds, aSet);
		WriterOBJ.write(ifs, new FileOutputStream("test.obj"));
	}
	
	
	@Test
	public void testLayoutTriangle1() throws Exception {
		double[] A = {0, 0, 1, 0};
		double[] B = {1, 0, 1, 0};
		
		double[] C = SphericalLayout.layoutTriangle(A, B, PI/2, PI/2, PI/2);
		ComplexProjective1 Cc = new ComplexProjective1(C[0], C[1], C[2], C[3]);
		Complex z = new Complex();
		Cc.projectTo(z);
		Assert.assertEquals(0, z.re, 1E-8);
		Assert.assertEquals(-1, z.im, 1E-8);
	}
	
	@Test
	public void testLayoutTriangle2() throws Exception {
		double[] A = {0, 0, 1, 0};
		double[] B = {1, 0, 1, 0};
		
		double[] C = SphericalLayout.layoutTriangle(A, B, PI/2, PI, PI/2);
		ComplexProjective1 Cc = new ComplexProjective1(C[0], C[1], C[2], C[3]);
		System.out.println(Cc);
		Complex z = new Complex();
		Cc.projectTo(z);
		System.out.println(z);
		Assert.assertEquals(1.0, z.abs(), 1E-8);
		Assert.assertEquals(-1, z.re, 1E-8);
		Assert.assertEquals(0, z.im, 1E-8);
	}
	
	@Test
	public void testAssignSphericalLogScaleRotation() throws Exception {
		ComplexProjective1 center = new ComplexProjective1(0, 0, 1, 0);
		Moebius M = new Moebius();
		M.assignSphericalLogScaleRotation(center, log(tan(PI/4)), 0);
		ComplexProjective1 R = new ComplexProjective1(1, 0, 1, 0);
		M.applyTo(R);
		Complex z = new Complex();
		R.projectTo(z);
		System.out.println(z);
		Assert.assertEquals(tan(PI/4), z.re, 1E-8);
	}
	
}

