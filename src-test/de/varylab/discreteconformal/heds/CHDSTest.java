/**
 * 
 */
package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;

import java.io.IOException;

import junit.framework.Assert;
import no.uib.cipr.matrix.DenseVector;

import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.unwrapper.UnwrapUtility;
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
			ConverterJR2Heds converter = new ConverterJR2Heds();
			hds = new CoHDS();
			AdapterSet a = new AdapterSet(new CoPositionAdapter());
			converter.ifs2heds(ifs, hds, a, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
		UnwrapUtility.prepareInvariantDataEuclidean(hds, new AdapterSet());
	}


	/**
	 * Test method for {@link de.varylab.discreteconformal.heds.CoHDS#conformalEnergy(no.uib.cipr.matrix.Vector, double[], no.uib.cipr.matrix.Vector, no.uib.cipr.matrix.Matrix)}.
	 */
	@Test
	public void testConformalEnergy() throws Exception {
		CEuclideanOptimizable opt = new CEuclideanOptimizable(hds);
		int n = opt.getDomainDimension();
		DenseVector u = new DenseVector(n);
		NewtonOptimizer optimizer = new NewtonOptimizer();
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CGS);
		optimizer.setError(1E-13);
		optimizer.setMaxIterations(20);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			Assert.fail(e.getMessage());
		}
		for (CoVertex v : hds.getVertices()) {
			if (v.getSolverIndex() < 0) {
				continue;
			}
			double aSum = 0.0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				aSum += e.getPreviousEdge().getAlpha();
			}
			Assert.assertEquals(2 * PI, aSum, 1E-13);
		}
	}

	
}
