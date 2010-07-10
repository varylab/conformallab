package de.varylab.discreteconformal.heds;

import java.io.IOException;

import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.SparseVector;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

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
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.unwrapper.CHyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.UnwrapUtility;

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
			AdapterSet a = new AdapterSet(new PositionAdapter());
			converter.ifs2heds(ifs, hds, a, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
		UnwrapUtility.prepareInvariantDataHyperbolic(hds);
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		System.out.println("CLayoutTest.tearDownAfterClass()");
	}

	@Before
	public void setUp() throws Exception {
		System.out.println("CLayoutTest.setUp()");
	}

	@After
	public void tearDown() throws Exception {
		System.out.println("CLayoutTest.tearDown()");
	}

	@Test
	public void testDoLayout() {
		System.out.println("CLayoutTest.testDoLayout()");
		int n = UnwrapUtility.prepareInvariantDataHyperbolic(hds);
		Vector u = new SparseVector(n);
		CHyperbolicLayout.doLayout(hds, hds.getVertex(0), u);
		
//		for (CoEdge e : hds.getPositiveEdges()) {
			//TODO figure out how a reasonable test looks like
//			Point s = e.getStartVertex().getTextureCoord();
//			Point t = e.getTargetVertex().getTextureCoord();
			
//			double l1 = Pn.distanceBetween(null, null, n)
//			double l2 = s.getTextureCoord().distanceTo(t.getTextureCoord());
//			Assert.assertEquals(l1, l2, 1E-3);
//		}
	}

	public static void main(String[] args) throws Exception{
		setUpBeforeClass();
		ConverterHeds2JR converter = new ConverterHeds2JR();
		AdapterSet a = new AdapterSet(new PositionAdapter());
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
