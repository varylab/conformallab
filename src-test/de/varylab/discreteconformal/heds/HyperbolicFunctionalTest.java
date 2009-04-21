package de.varylab.discreteconformal.heds;

import java.io.IOException;

import junit.framework.Assert;

import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedge.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.plugin.adapter.PositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicApplication;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;
import de.varylab.jtao.Tao.Method;


public class HyperbolicFunctionalTest {

	private static CoHDS 	
		hds = null;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			Input in = new Input("Obj File", EuclideanLayoutTest.class.getResourceAsStream("torusCoarse.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds<CoVertex, CoEdge, CoFace> converter = new ConverterJR2Heds<CoVertex, CoEdge, CoFace>(CoVertex.class, CoEdge.class, CoFace.class);
			hds = new CoHDS();
			converter.ifs2heds(ifs, hds, new PositionAdapter());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	@Test
	public void testGradient() throws Exception {
		hds.prepareInvariantDataHyperbolic();
		
		Tao.Initialize();
		Tao optimizer = new Tao(Method.CG);
		
		CHyperbolicApplication app = new CHyperbolicApplication(hds);
		optimizer.setApplication(app);
		
		Vec testPoint = new Vec(app.getDomainDimension());
		double grad = optimizer.testGradient(testPoint, true);
		Assert.assertEquals("Gradient", 0.0, grad, 1E-5);

	}
	
	@Test
	public void testHessian() throws Exception {
		hds.prepareInvariantDataHyperbolic();
		
		Tao.Initialize();
		Tao optimizer = new Tao(Method.CG);
		
		CHyperbolicApplication app = new CHyperbolicApplication(hds);
		optimizer.setApplication(app);
		
		Vec testPoint = new Vec(app.getDomainDimension());
		double hess = optimizer.testHessian(testPoint, true);
		Assert.assertEquals("Hessian", 0.0, hess, 1E-5);

	}
	
}
