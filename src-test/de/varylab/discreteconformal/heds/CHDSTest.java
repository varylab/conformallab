/**
 * 
 */
package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import geom3d.Point;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseVector;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jtem.halfedge.jReality.converter.ConverterJR2Heds;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.math.optimization.NotConvergentException;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer.Solver;
import de.varylab.discreteconformal.math.optimization.stepcontrol.ArmijoStepController;

/**
 * @author sechel
 *
 */
public class CHDSTest extends TestCase{

	private static CHDS 	
		hds = null;
	
	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		File file = new File("data/test01.obj");
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c =reader.read(file);
		IndexedFaceSet ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
		ConverterJR2Heds<CVertex, CEdge, CFace> converter = new ConverterJR2Heds<CVertex, CEdge, CFace>(CVertex.class, CEdge.class, CFace.class);
		hds = new CHDS();
		converter.ifs2heds(ifs, hds, new PositionAdapter());
	}

	/**
	 * @throws java.lang.Exception
	 */
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {

	}

	/**
	 * @throws java.lang.Exception
	 */
	@After
	public void tearDown() throws Exception {
	}

	/**
	 * Test method for {@link de.varylab.discreteconformal.heds.CHDS#conformalEnergy(no.uib.cipr.matrix.Vector, double[], no.uib.cipr.matrix.Vector, no.uib.cipr.matrix.Matrix)}.
	 */
	@Test
	public void testConformalEnergy() throws Exception {
		hds.prepareInvariantData();
		int n = hds.getDomainDimension();
		DenseVector u = new DenseVector(n);
		NewtonOptimizer optimizer = new NewtonOptimizer();
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(1E-10);
		try {
			optimizer.minimize(u, hds);
		} catch (NotConvergentException e) {
			e.printStackTrace();
		}

		
		Map<CEdge, Double> aMpa = hds.calculateAlphas(u);
		for (CVertex v : hds.getVertices()) {
			if (!hds.isVariable(v)) 
				continue;
			double aSum = 0.0;
			for (CEdge e : HalfEdgeUtils.incomingEdges(v))
				aSum += aMpa.get(e.getPreviousEdge());
			assertEquals(2 * PI, aSum, 1E-8);
		}
		
	}

	
	public void testTriangleEnergyAndAlphas() throws Exception {
		CHDS hds = new CHDS();
		CVertex v0 = hds.addNewVertex();
		CVertex v1 = hds.addNewVertex();
		CVertex v2 = hds.addNewVertex();
		HalfEdgeUtils.constructFaceByVertices(hds, v0, v1, v2);
		v0.setPosition(new Point(0,0,0));
		v1.setPosition(new Point(1,0,0));
		v2.setPosition(new Point(1,1,0));
		hds.prepareInvariantData();
		Map<CEdge, Double> aMap = new HashMap<CEdge, Double>();
		hds.triangleEnergyAndAlphas(null, hds.getFace(0), aMap);
		for (CEdge e : aMap.keySet()) {
			System.err.println("Vertex " + e.getNextEdge().getTargetVertex() + ": " + aMap.get(e));
		}
	}
	
	
}
