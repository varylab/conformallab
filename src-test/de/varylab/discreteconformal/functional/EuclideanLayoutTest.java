package de.varylab.discreteconformal.functional;

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
import de.jreality.util.SceneGraphUtility;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class EuclideanLayoutTest {

	private static CoHDS 	
		hds = null;
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
	public EuclideanFunctional<CoVertex, CoEdge, CoFace>
		fun = new EuclideanFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		ReaderOBJ reader = new ReaderOBJ(); 
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		try {
			Input in = new Input("Obj File", EuclideanLayoutTest.class.getResourceAsStream("tetraflat.obj"));
			c = reader.read(in);
			ifs = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c);
			ConverterJR2Heds converter = new ConverterJR2Heds();
			hds = new CoHDS();
			converter.ifs2heds(ifs, hds, a, null);;
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Test
	public void testDoLayout() {
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		hds.normalizeCoordinates();
		int n = UnwrapUtility.prepareInvariantDataEuclidean(fun, hds, a);
		Vector u = new SparseVector(n);
		EuclideanLayout.doLayout(hds, fun, u);
		for (CoEdge e : hds.getEdges()) {
			CoVertex s = e.getStartVertex();
			CoVertex t = e.getTargetVertex();
			double l1 = Pn.distanceBetween(s.P, t.P, Pn.EUCLIDEAN);
			double l2 = Pn.distanceBetween(s.T, t.T, Pn.EUCLIDEAN);
			Assert.assertEquals(l1, l2, 1E-11);
		}
	}
	
}
