package de.varylab.discreteconformal.heds;

import static de.jreality.scene.data.Attribute.LABELS;

import java.io.File;
import java.io.IOException;

import junit.framework.Assert;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.SparseVector;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.data.StringArray;
import de.jreality.ui.viewerapp.ViewerApp;
import de.jtem.halfedge.jReality.converter.ConverterHeds2JR;
import de.jtem.halfedge.jReality.converter.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.unwrap.CDiskLayout;

public class CLayoutTest {

	private static CHDS 	
		hds = null;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		System.out.println("CLayoutTest.setUpBeforeClass()");
		System.out.println("CHDSTest.setUpBeforeClass()");
		File file = new File("data/planar04.obj");
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			c =reader.read(file);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds<CVertex, CEdge, CFace> converter = new ConverterJR2Heds<CVertex, CEdge, CFace>(CVertex.class, CEdge.class, CFace.class);
			hds = new CHDS();
			converter.ifs2heds(ifs, hds, new PositionAdapter());
		} catch (IOException e) {
			e.printStackTrace();
		}
		hds.prepareInvariantData();
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
		int n = hds.getDomainDimension();
		Vector u = new SparseVector(n);
		CDiskLayout.doLayout(hds, u, null);
		
		for (CEdge e : hds.getPositiveEdges()) {
			CVertex s = e.getStartVertex();
			CVertex t = e.getTargetVertex();
			double l1 = s.getPosition().distanceTo(t.getPosition());
			double l2 = s.getTextureCoord().distanceTo(t.getTextureCoord());
			Assert.assertEquals(l1, l2, 1E-8);
		}
	}

	public static void main(String[] args) throws Exception{
		setUpBeforeClass();
		ConverterHeds2JR<CVertex, CEdge, CFace> converter = new ConverterHeds2JR<CVertex, CEdge, CFace>();
		IndexedFaceSet ifs = converter.heds2ifs(hds, new PositionAdapter());
		String[] vertexLabels = new String[hds.numVertices()];
		for (int  i = 0; i < hds.numVertices(); i++) {
			vertexLabels[i] = "" + i;
		}
		ifs.setVertexAttributes(LABELS, new StringArray(vertexLabels));
		SceneGraphComponent c = new SceneGraphComponent();
		c.setGeometry(ifs);
		ViewerApp.display(c);
	}
	
	
}
