package de.varylab.discreteconformal.heds;

import java.io.IOException;

import junit.framework.Assert;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.SparseVector;

import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class EuclideanLayoutTest {

	private static CoHDS 	
		hds = null;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		ReaderOBJ reader = new ReaderOBJ(); 
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			Input in = new Input("Obj File", EuclideanLayoutTest.class.getResourceAsStream("planar01.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds converter = new ConverterJR2Heds();
			hds = new CoHDS();
			AdapterSet a = new AdapterSet(new CoPositionAdapter());
			converter.ifs2heds(ifs, hds, a, null);;
		} catch (IOException e) {
			e.printStackTrace();
		}
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		UnwrapUtility.prepareInvariantDataEuclidean(hds, a);
	}

	@Test
	public void testDoLayout() {
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		int n = UnwrapUtility.prepareInvariantDataEuclidean(hds, a);
		Vector u = new SparseVector(n);
		EuclideanLayout.doLayout(hds, u);
		
		for (CoEdge e : hds.getPositiveEdges()) {
			CoVertex s = e.getStartVertex();
			CoVertex t = e.getTargetVertex();
			double l1 = Pn.distanceBetween(s.P, t.P, Pn.EUCLIDEAN);
			double l2 = Pn.distanceBetween(s.T, t.T, Pn.EUCLIDEAN);
			Assert.assertEquals(l1, l2, 1E-6);
		}
	}
	
}
