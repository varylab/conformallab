package de.varylab.discreteconformal.heds;

import static de.jreality.scene.data.Attribute.LABELS;

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
import de.jreality.util.Input;
import de.jtem.halfedge.jreality.ConverterHeds2JR;
import de.jtem.halfedge.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.unwrapper.CEuclideanLayout;

public class EuclideanLayoutTest {

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
			Input in = new Input("Obj File", EuclideanLayoutTest.class.getResourceAsStream("planar01.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds<CoVertex, CoEdge, CoFace> converter = new ConverterJR2Heds<CoVertex, CoEdge, CoFace>(CoVertex.class, CoEdge.class, CoFace.class);
			hds = new CoHDS();
			converter.ifs2heds(ifs, hds, new PositionAdapter());
		} catch (IOException e) {
			e.printStackTrace();
		}
		hds.prepareInvariantDataEuclidean();
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
		int n = hds.prepareInvariantDataEuclidean();
		Vector u = new SparseVector(n);
		CEuclideanLayout.doLayout(hds, u);
		
		for (CoEdge e : hds.getPositiveEdges()) {
			CoVertex s = e.getStartVertex();
			CoVertex t = e.getTargetVertex();
			double l1 = s.getPosition().distanceTo(t.getPosition());
			double l2 = s.getTextureCoord().distanceTo(t.getTextureCoord());
			Assert.assertEquals(l1, l2, 1E-6);
		}
	}

	public static void main(String[] args) throws Exception{
		setUpBeforeClass();
		ConverterHeds2JR<CoVertex, CoEdge, CoFace> converter = new ConverterHeds2JR<CoVertex, CoEdge, CoFace>();
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