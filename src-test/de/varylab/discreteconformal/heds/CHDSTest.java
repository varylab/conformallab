/**
 * 
 */
package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;

import java.io.IOException;

import junit.framework.Assert;
import no.uib.cipr.matrix.DenseVector;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedge.jreality.ConverterJR2Heds;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.math.CEuclideanOptimizable;
import de.varylab.discreteconformal.math.optimization.NotConvergentException;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer.Solver;
import de.varylab.discreteconformal.math.optimization.stepcontrol.ArmijoStepController;

/**
 * @author sechel
 *
 */
public class CHDSTest {

	private static CHDS 	
		hds = null;
	
	
	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass 
	public static void setUpBeforeClass() throws Exception {
		System.out.println("CHDSTest.setUpBeforeClass()");
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		try {
			Input in = new Input("Obj File", CLayoutTest.class.getResourceAsStream("cathead.obj"));
			c =reader.read(in);
			IndexedFaceSet ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds<CVertex, CEdge, CFace> converter = new ConverterJR2Heds<CVertex, CEdge, CFace>(CVertex.class, CEdge.class, CFace.class);
			hds = new CHDS();
			converter.ifs2heds(ifs, hds, new PositionAdapter());
		} catch (IOException e) {
			e.printStackTrace();
		}
		hds.prepareInvariantData();
	}

	/**
	 * @throws java.lang.Exception
	 */
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		System.out.println("CHDSTest.tearDownAfterClass()");
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		System.out.println("CHDSTest.setUp()");
	}

	/**
	 * @throws java.lang.Exception
	 */
	@After
	public void tearDown() throws Exception {
		System.out.println("CHDSTest.tearDown()");
	}

	/**
	 * Test method for {@link de.varylab.discreteconformal.heds.CHDS#conformalEnergy(no.uib.cipr.matrix.Vector, double[], no.uib.cipr.matrix.Vector, no.uib.cipr.matrix.Matrix)}.
	 */
	@Test
	public void testConformalEnergy() throws Exception {
		System.out.println("CHDSTest.testConformalEnergy()");
		CEuclideanOptimizable opt = new CEuclideanOptimizable(hds);
		int n = opt.getDomainDimension();
		DenseVector u = new DenseVector(n);
		NewtonOptimizer optimizer = new NewtonOptimizer();
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(1E-10);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			e.printStackTrace();
		}

		
		for (CVertex v : hds.getVertices()) {
			if (v.getSolverIndex() < 0) {
				continue;
			}
			double aSum = 0.0;
			for (CEdge e : HalfEdgeUtils.incomingEdges(v)) {
				aSum += e.getPreviousEdge().getAlpha();
			}
			Assert.assertEquals(2 * PI, aSum, 1E-8);
		}
		
	}

	
}
