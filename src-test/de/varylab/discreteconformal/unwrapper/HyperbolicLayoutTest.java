package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.cosh;
import static java.lang.Math.sinh;

import java.io.IOException;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.plugin.JRViewer;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.data.Attribute;
import de.jreality.scene.data.StringArray;
import de.jreality.util.Input;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterHeds2JR;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;

public class HyperbolicLayoutTest {

	private static CoHDS 	
		hds = null;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		System.out.println("CLayoutTest.setUpBeforeClass()");
		System.out.println("CHDSTest.setUpBeforeClass()");
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			Input in = new Input("Obj File", HyperbolicLayoutTest.class.getResourceAsStream("tetraflat.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds converter = new ConverterJR2Heds();
			hds = new CoHDS();
			AdapterSet a = new AdapterSet(new CoPositionAdapter());
			converter.ifs2heds(ifs, hds, a, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Test
	public void testLayoutTriangle() throws Exception {
		double[] ZERO = {0,0,0,1};
		double[] A = {-sinh(1.0)/cosh(1.0), 0.0, 0.0, 1.0};
		double[] B = { 0.0, 0.0, 0.0, 1.0};
		double alpha = Math.PI / 4;

		// move points
		double[] O = {0.2, 0.2, 0.0, 1};
		MatrixBuilder Om = MatrixBuilder.hyperbolic();
		Om.rotate(0.2, 0, 0, 1);
		Om.translateFromTo(ZERO, O);
		
		Om.getMatrix().transformVector(A);
		Om.getMatrix().transformVector(B);
		
		// check distance
		double dAB = Pn.distanceBetween(A, B, Pn.HYPERBOLIC);
		Assert.assertEquals(1.0, dAB, 1E-8);
		
		// stretch rotation
		MatrixBuilder mb = MatrixBuilder.hyperbolic();
		mb.translateFromTo(ZERO, B);
		mb.rotate(alpha, 0, 0, 1);
		mb.scale((sinh(1.1)/cosh(1.1)) / (sinh(1.0)/cosh(1.0)));
		mb.translate(B, ZERO);
		
		// check B identity
		double[] Bcheck = B.clone();
		mb.getMatrix().transformVector(Bcheck);
		Assert.assertArrayEquals(B, Bcheck, 1E-8);
		
		// create C
		double[] C = A.clone();
		mb.getMatrix().transformVector(C);
		double dBC = Pn.distanceBetween(B, C, Pn.HYPERBOLIC);
		double dAC = Pn.distanceBetween(A, C, Pn.HYPERBOLIC);
		
		Assert.assertEquals(1.1, dBC, 1E-8);
		
		double[] CP = HyperbolicLayout.layoutTriangle(A, B, alpha, dBC, dAC);
		Pn.dehomogenize(C, C);
		Pn.dehomogenize(CP, CP);
		Assert.assertArrayEquals(C, CP, 1e-10);
	}
	
	@Test@Ignore
	public void testDoLayout() {
		throw new RuntimeException("Not implemented correctly yet!");
//		System.out.println("CLayoutTest.testDoLayout()");
//		AdapterSet a = AdapterSet.createGenericAdapters();
//		a.add(new CoPositionAdapter());
//		a.add(new CoPositionAdapter());
//		a.add(new CoTexturePositionAdapter());
//		int n = UnwrapUtility.prepareInvariantDataHyperbolic(hds, a);
//		Vector u = new SparseVector(n);
//		HyperbolicLayout.doLayout(hds, hds.getVertex(0), u);
//		
//		for (CoEdge e : hds.getPositiveEdges()) {
//			double[] st = e.getStartVertex().T;
//			double[] tt = e.getTargetVertex().T;
//			double l1 = HyperbolicLayout.getNewLength(e, u);
//			double l2 = Pn.distanceBetween(st, tt, Pn.HYPERBOLIC);
//			Assert.assertEquals(l1, l2, 1E-3);
//		}
	}

	public static void main(String[] args) throws Exception{
		setUpBeforeClass();
		ConverterHeds2JR converter = new ConverterHeds2JR();
		AdapterSet a = new AdapterSet(new CoPositionAdapter());
		IndexedFaceSet ifs = converter.heds2ifs(hds, a);
		String[] vertexLabels = new String[hds.numVertices()];
		for (int  i = 0; i < hds.numVertices(); i++) {
			vertexLabels[i] = "" + i;
		}
		ifs.setVertexAttributes(Attribute.LABELS, new StringArray(vertexLabels));
		SceneGraphComponent c = new SceneGraphComponent();
		c.setGeometry(ifs);
		JRViewer.display(c);
	}
	
	
}
