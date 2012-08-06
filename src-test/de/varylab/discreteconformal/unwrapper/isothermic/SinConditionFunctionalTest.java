package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.PICKABLE;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import junit.framework.Assert;

import org.junit.Before;
import org.junit.Test;

import de.jreality.geometry.Primitives;
import de.jreality.geometry.QuadMeshFactory;
import de.jreality.plugin.JRViewer;
import de.jreality.scene.Appearance;
import de.jreality.scene.proxy.scene.SceneGraphComponent;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.NormType;
import de.jtem.jpetsc.SNES;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.Method;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.TestUtility;

public class SinConditionFunctionalTest {

	public SinConditionApplication<CoVertex, CoEdge, CoFace, CoHDS>
		fun = null;
	
	static {
		NativePathUtility.set("native");
		String[] args = {
				"-snes_type", "ls",
				"-snes_test_display",
				"-pc_factor_shift_nonzero", "1.0e-10"
		};
		Tao.Initialize("Sinus Condition Test", args, false);
	}
	
	@Before
	public void createFunctional() {
		CoHDS hds = new CoHDS();
		CoFace f = HalfEdgeUtils.addNGon(hds, 3);
		CoVertex v = TopologyAlgorithms.splitFace(f);

		Map<CoEdge, Double> alphaMap = new HashMap<CoEdge, Double>();
		List<CoEdge> eIn = HalfEdgeUtils.incomingEdges(v);
		alphaMap.put(eIn.get(0), -3*PI/8);
		alphaMap.put(eIn.get(0).getPreviousEdge(), -PI/4);
		alphaMap.put(eIn.get(1), -0.2);
		alphaMap.put(eIn.get(1).getPreviousEdge(), PI/4);
		alphaMap.put(eIn.get(2), 3*PI/8);
		alphaMap.put(eIn.get(2).getPreviousEdge(), PI/2);
		TestUtility.completeOpposites(alphaMap);
		
		fun = new SinConditionApplication<CoVertex, CoEdge, CoFace, CoHDS>(hds);
		fun.initialize(alphaMap, false);
	}
	
	
	@Test
	public void testSolveCG() throws Exception {
		Vec init = fun.getSolutionVec();
		Vec startVec = new Vec(init.getSize());
		init.copy(startVec);
		
		Assert.assertTrue("energy is large at the beginning", 1E-3 < fun.evaluateObjective(init));
		
		fun.solveEnergyMinimzation(200, 1E-10, Method.CG);
		Vec solution = fun.getSolutionVec();
		Assert.assertTrue("energy is small after minimizatin", 1E-8 > fun.evaluateObjective(solution));
		
		double dif2 = 0.0;
		for (int i = 0; i < startVec.getSize(); i++) {
			double start = startVec.getValue(i);
			double sol = init.getValue(i);
			dif2 += (start-sol)*(start-sol);
		}
		Assert.assertEquals("solution start proximity", 1E-3, dif2, 1E-2);
	}
	
	@Test
	public void testSolveSNES() throws Exception {
		Vec init = fun.getSolutionVec();
		Vec startVec = new Vec(init.getSize());
		init.copy(startVec);
		Vec f = new Vec(fun.getDimension());
		
		fun.evaluateFunction(init, f);
		double res = f.norm(NormType.NORM_FROBENIUS);
		Assert.assertTrue("residual is large at the beginning", 1E-3 < res);
		
		fun.solveSNES(200, 1E-10);
		Vec solution = fun.getSolutionVec();
		fun.evaluateFunction(solution, f);
		res = f.norm(NormType.NORM_FROBENIUS);
		Assert.assertTrue("residual is small after minimizatin", 1E-8 > res);
	}
	
	
	@Test
	public void testSNESJacobian() throws Exception {
		SNES snes = SNES.create();
		Vec x = fun.getSolutionVec();
		Vec f = new Vec(x.getSize());
		Mat J = new Mat(x.getSize(), x.getSize());
		J.assemble();
		snes.setFromOptions();
		snes.setFunction(fun, f);
		snes.setJacobian(fun, J, J);
		fun.evaluateJacobian(x, J, J);
		boolean check = SNES.testJacobian(fun, J, x, f);
		Assert.assertTrue("SNES Jacobian seems correct", check);
	}
	
	@Test
	public void testCGGradient() throws Exception {
		Vec init = fun.getSolutionVec();
		Vec Gfd = new Vec(init.getSize());
		Vec Ghc = new Vec(init.getSize());
		Gfd.zeroEntries();
		Ghc.zeroEntries();
		Gfd.assemble();
		Ghc.assemble();
		fun.computeGradient(init, Ghc);
		TestUtility.calculateFDGradient(fun, init, Gfd);
		for (int i = 0; i < init.getSize(); i++) {
			double fd = Gfd.getValue(i);
			double hc = Ghc.getValue(i);
			Assert.assertEquals(fd, hc, 1E-6);
		}
	}
	
	
	@Test
	public void testCGHessian() throws Exception {
		Vec init = fun.getSolutionVec();
		Mat Hfd = fun.createHessianTemplate();
		Mat Hhc = fun.createHessianTemplate();
		fun.evaluateHessian(init, Hhc, Hhc);
		TestUtility.calculateFDHessian(fun, init, Hfd);
		System.out.println("HC:\n" + Hhc);
		System.out.println("FD:\n" + Hfd);
		for (int i = 0; i < init.getSize(); i++) {
			for (int j = 0; j < init.getSize(); j++) {
				double fd = Hfd.getValue(i, j);
				double hc = Hhc.getValue(i, j);
				Assert.assertEquals(fd, hc, 1E-5);
			}
		}
	}
	
	
	public static void main(String[] args) throws Exception {
		SinConditionFunctionalTest test = new SinConditionFunctionalTest();
		test.createFunctional();
		test.testCGGradient();
//		test.testEnergy();
	}
	
	public void testEnergy() throws Exception {
		createFunctional();
		Tao tao = Tao.createFiniteDifferenceTester(true, false);
		tao.setFromOptions();
		tao.setApplication(fun);
		tao.solve();
		
		int e1Index = 0;
		int e2Index = 1;
		double scale = 1E-1;
		
		double[] xArr = {1.550848732167379, 0.7726856928589729, -0.8016039900623643, -0.1910292974241263, -1.1581496504686548, 1.19804483972369};
		Vec x = fun.getSolutionVec();
		for (int i = 0; i < xArr.length; i++) {
			x.setValue(i, xArr[i], INSERT_VALUES);
		}
		x.assemble();
		
		int uLine = 500;
		int vline = 500;
		
		double[][] verts = new double[uLine*vline][];
		Vec g = new Vec(xArr.length);
		g.assemble();
		
		for (int i = 0; i < uLine; i++) {
			for (int j = 0; j < vline; j++) {
				double x1 = 2*PI * (i/(double)(uLine-1)) - PI;
				double x2 = 2*PI * (j/(double)(vline-1)) - PI;
				x.setValue(e1Index, x1, INSERT_VALUES);
				x.setValue(e2Index, x2, INSERT_VALUES);
				double val = fun.evaluateObjective(x);
//				fun.evaluateGradient(x, g);
//				fun.defaultComputeGradient(x, g);
//				val = g.getValue(1);
				val *= scale;
				val = Math.min(val, 0.5);
				val = Math.max(val, -0.5);
				if (Double.isNaN(val)) {
					val = 0.0;
				}
				verts[i*uLine + j] = new double[] {x1, x2, val};
			}
		}
		QuadMeshFactory qmf = new QuadMeshFactory();
		qmf.setULineCount(uLine);
		qmf.setVLineCount(vline);
		qmf.setVertexCoordinates(verts);
		qmf.setGenerateEdgesFromFaces(true);
		qmf.setGenerateFaceNormals(true);
		qmf.setGenerateVertexNormals(true);
		qmf.update();
		
		x.setValue(e1Index, xArr[e1Index], INSERT_VALUES);
		x.setValue(e2Index, xArr[e2Index], INSERT_VALUES);
		double val = fun.evaluateObjective(x);
		System.out.println("value: " + val);
		double[] markerPos = {xArr[e1Index], xArr[e2Index], scale * val};
		
		SceneGraphComponent root = new SceneGraphComponent();
		SceneGraphComponent graph = new SceneGraphComponent();
		SceneGraphComponent marker = new SceneGraphComponent();
		graph.setGeometry(qmf.getIndexedFaceSet());
		Appearance graphApp = new Appearance();
		graphApp.setAttribute(VERTEX_DRAW, false);
		graphApp.setAttribute(EDGE_DRAW, false);
		graphApp.setAttribute(POLYGON_SHADER + "." + PICKABLE, false);
		graph.setAppearance(graphApp);
		marker.setGeometry(Primitives.point(markerPos));
		root.addChild(graph);
		root.addChild(marker);
		
		JRViewer viewer = new JRViewer();
		viewer.addBasicUI();
		viewer.addContentUI();
		viewer.setContent(root);
		viewer.startup();
	}
	
}
