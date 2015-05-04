package de.varylab.discreteconformal.plugin;

import java.util.logging.Logger;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.functional.HyperIdealGenerator;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperIdealApplication;

public class HyperIdealPluginTest {

	private static Logger
		log = Logger.getLogger(HyperIdealPluginTest.class.getName());
	
	@BeforeClass
	public static void initPetsc() {
		NativePathUtility.set("native");
		Tao.Initialize();
		LoggingUtility.initLogging();
	}
	
	@Test
	public void testSetSolverIndicesForCurve() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonHyperelliptic();
		Selection s = new Selection(hds.getVertex(0), hds.getVertex(1), hds.getVertex(2), hds.getVertex(3), hds.getVertex(6), hds.getVertex(7));
		HyperIdealPlugin.setSolverIndicesForCurve(hds, s);
		for (CoVertex v : hds.getVertices()) {
			int i = v.getIndex();
			if (i == 0 || i == 1 || i == 2 || i == 3 || i == 6 || i == 7) {
				Assert.assertTrue(v.getSolverIndex() >= 0);
			} else {
				Assert.assertTrue(v.getSolverIndex() < 0);
			}
		}
		for (CoEdge e : hds.getPositiveEdges()) {
			Assert.assertTrue(e.getSolverIndex() >= 0);
			Assert.assertTrue(e.getOppositeEdge().getSolverIndex() >= 0);
		}
	}
	
	
	@Test
	public void testCreateTaoApplicationForCurve() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonHyperelliptic();
		Selection s = new Selection(hds.getVertex(0), hds.getVertex(1), hds.getVertex(2), hds.getVertex(3), hds.getVertex(6), hds.getVertex(7));
		HyperIdealPlugin.setSolverIndicesForCurve(hds, s);
		CHyperIdealApplication app =  HyperIdealPlugin.createTaoApplication(hds);
		Assert.assertNotNull(app.getSolutionVec());
		Assert.assertEquals(6 + hds.numEdges()/2, app.getDomainDimension());
	}
	
	@Test
	public void testOptimizeHyperIdealApplication() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonHyperelliptic();
		Selection s = new Selection(hds.getVertex(0), hds.getVertex(1), hds.getVertex(2), hds.getVertex(3), hds.getVertex(6), hds.getVertex(7));
		HyperIdealPlugin.setSolverIndicesForCurve(hds, s);
		CHyperIdealApplication app =  HyperIdealPlugin.createTaoApplication(hds);
		Tao optimizer = HyperIdealPlugin.optimizeHyperIdealApplication(app, 1E-6);
		log.info(optimizer.getSolutionStatus().toString());		
		Assert.assertEquals(ConvergenceFlags.CONVERGED_ATOL, optimizer.getSolutionStatus().reason);
	}
	
}
