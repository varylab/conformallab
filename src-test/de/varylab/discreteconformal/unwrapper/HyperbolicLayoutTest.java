package de.varylab.discreteconformal.unwrapper;

import java.io.IOException;

import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import cern.colt.Arrays;
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
		double[] A = {0,0,0,1};
		double[] B = {0.1,0,0,1};
		double[] C = HyperbolicLayout.layoutTriangle(A, B, Math.PI/5, 1, 1);
		System.out.println(Arrays.toString(C));
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
