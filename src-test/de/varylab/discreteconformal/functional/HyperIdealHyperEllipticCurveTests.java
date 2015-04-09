package de.varylab.discreteconformal.functional;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static java.lang.Math.PI;

import java.io.InputStream;
import java.util.Random;
import java.util.logging.Logger;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.HalfedgeEmbedding;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperIdealApplication;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class HyperIdealHyperEllipticCurveTests {

	private static Logger
		log = Logger.getLogger(HyperIdealHyperellipticUtility.class.getName());
	
	@BeforeClass
	public static void initPetsc() {
		NativePathUtility.set("native");
		Tao.Initialize();
		LoggingUtility.initLogging();
		PETSc.optionsSetValue("-tao_lmm_vectors", "20");
		PETSc.optionsSetValue("-tao_lmm_scale_type", "broyden");
		PETSc.optionsSetValue("-tao_lmm_broyden_phi", "0.125");
		PETSc.optionsSetValue("-tao_lmm_rescale_type", "scalar");
		PETSc.optionsSetValue("-tao_lmm_rescale_history", "5");
		PETSc.optionsSetValue("-tao_lmm_rescale_alpha", "5.0");
		PETSc.optionsSetValue("-tao_lmm_rescale_beta", "0.5");
		PETSc.optionsSetValue("-tao_lmm_limit_type", "relative");
		PETSc.optionsSetValue("-tao_lmm_limit_mu", "1.0");
		PETSc.optionsSetValue("-tao_lmm_limit_nu", "1.0");
	}
	
	@Test
	public void testCreateDataOnHyperEllipticCurveLawson() throws Exception {
		InputStream in = getClass().getResourceAsStream("lawson_curve_source.xml");
		HalfedgeEmbedding he = DataIO.readConformalData(HalfedgeEmbedding.class, in);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
		CoHDS hds = new CoHDS();
		AdapterSet a = new ConformalAdapterSet();
		DataUtility.toHalfedge(he, a, Position.class, hds, cutInfo);
		HyperIdealHyperellipticUtility.calculateCircleIntersections(hds);
		
		// vertex data
		int index = 0;
		for (CoVertex v : hds.getVertices()) {
			v.setTheta(2*PI);
			double thetaSum = 0.0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				thetaSum += e.getTheta();
			}
			switch (v.getIndex()) {
			case 0: case 1: case 2: case 3: case 6: case 7:
				Assert.assertEquals(4*PI, thetaSum, 1e-8);
				v.setSolverIndex(index++);
				break;
			default:
				Assert.assertEquals(2*PI, thetaSum, 1e-8);
				v.setSolverIndex(-1);	
			}
		}
		
		// optimize
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		app.setFromOptions();
		int n = app.getDomainDimension();
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vec u = new Vec(n);
		for (int i = 0; i < n; i++) {
			u.setValue(i, 0.1 + 0.01*Math.abs(rnd.nextDouble()), INSERT_VALUES);
		}
		app.setInitialSolutionVec(u);
		Vec lowerBounds = new Vec(n);
		Vec upperBounds = new Vec(n);
		lowerBounds.set(-Double.MAX_VALUE);
		for (int i = 0; i < 6; i++) {
			lowerBounds.setValue(i, 1E-12, INSERT_VALUES);
		}
		upperBounds.set(Double.MAX_VALUE);
		app.setVariableBounds(lowerBounds, upperBounds);
		log.info("start   : " + u.toString());
		
		Tao optimizer = new Tao(Tao.Method.BLMVM);
		optimizer.setFromOptions();
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(1E-7, 0, 0); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(50);
		optimizer.solve();
		log.info(optimizer.getSolutionStatus().toString());
		Assert.assertEquals(ConvergenceFlags.CONVERGED_ATOL, optimizer.getSolutionStatus().reason);
		UnwrapUtility.logSolutionStatus(optimizer, log);
	}
	
}
