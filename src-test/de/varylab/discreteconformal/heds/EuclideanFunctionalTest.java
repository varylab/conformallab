package de.varylab.discreteconformal.heds;

import java.io.IOException;

import junit.framework.Assert;

import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;
import de.varylab.jtao.Tao.Method;


public class EuclideanFunctionalTest {

	private static CoHDS 	
		hds = null;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			Input in = new Input("Obj File", EuclideanLayoutTest.class.getResourceAsStream("tetrahedron.obj"));
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
	public void testGradient() throws Exception {
		UnwrapUtility.prepareInvariantDataEuclidean(hds, new AdapterSet());
		
		Tao.Initialize();
		Tao optimizer = new Tao(Method.CG);
		
		CEuclideanApplication app = new CEuclideanApplication(hds);
		optimizer.setApplication(app);
		
		Vec testPoint = new Vec(app.getDomainDimension());
		double grad = optimizer.testGradient(testPoint, true);
		Assert.assertEquals("Gradient", 0.0, grad, 1E-5);

	}
	
	@Test
	public void testHessian() throws Exception {
		UnwrapUtility.prepareInvariantDataEuclidean(hds, new AdapterSet());
		
		Tao.Initialize();
		Tao optimizer = new Tao(Method.CG);
		
		CEuclideanApplication app = new CEuclideanApplication(hds);
		optimizer.setApplication(app);
		
		Vec testPoint = new Vec(app.getDomainDimension());
		double hess = optimizer.testHessian(testPoint, true);
		Assert.assertEquals("Hessian", 0.0, hess, 1E-5);
	}
	
}
