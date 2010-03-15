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
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

/**
 * @author sechel
 *
 */
public class CHDSTest {

	private static CoHDS 	
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
			Input in = new Input("Obj File", EuclideanLayoutTest.class.getResourceAsStream("cathead.obj"));
			c =reader.read(in);
			IndexedFaceSet ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds<CoVertex, CoEdge, CoFace> converter = new ConverterJR2Heds<CoVertex, CoEdge, CoFace>(CoVertex.class, CoEdge.class, CoFace.class);
			hds = new CoHDS();
			converter.ifs2heds(ifs, hds, new PositionAdapter());
		} catch (IOException e) {
			e.printStackTrace();
		}
		hds.prepareInvariantDataEuclidean();
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
	 * Test method for {@link de.varylab.discreteconformal.heds.CoHDS#conformalEnergy(no.uib.cipr.matrix.Vector, double[], no.uib.cipr.matrix.Vector, no.uib.cipr.matrix.Matrix)}.
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

		
		for (CoVertex v : hds.getVertices()) {
			if (v.getSolverIndex() < 0) {
				continue;
			}
			double aSum = 0.0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				aSum += e.getPreviousEdge().getAlpha();
			}
			Assert.assertEquals(2 * PI, aSum, 1E-8);
		}
		
	}

	
}
