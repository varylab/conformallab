package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static java.lang.Math.PI;

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
import de.jreality.shader.CommonAttributes;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.Method;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class SinConditionFunctionalTest {

	public SinConditionFunctional<CoVertex, CoEdge, CoFace, CoHDS>
		fun = null;
	
	
	static {
		NativePathUtility.set("native");
		String[] taoCommand = new String[] {
			"-tao_nm_lamda", "0.01", 
			"-tao_nm_mu", "1.0"
		};
		Tao.Initialize("Sinus Condition Test", taoCommand, false);
	}
	
	@Before
	public void createFunctional() {
		CoHDS hds = new CoHDS();
		CoFace f = HalfEdgeUtils.addNGon(hds, 3);
		CoVertex v = TopologyAlgorithms.splitFace(f);
		Map<Integer, Integer> edgeMap = IsothermicUtility.createUndirectedEdgeMap(hds);
		
		Vec init = new Vec(6);
		List<CoEdge> eIn = HalfEdgeUtils.incomingEdges(v);
		init.setValue(edgeMap.get(eIn.get(0).getIndex()), -3*PI/8, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(0).getPreviousEdge().getIndex()), -PI/4, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(1).getIndex()), 0.2, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(1).getPreviousEdge().getIndex()), PI/4, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(2).getIndex()), 3*PI/8, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(2).getPreviousEdge().getIndex()), PI/2, INSERT_VALUES);
		init.assemble();
		fun = new SinConditionFunctional<CoVertex, CoEdge, CoFace, CoHDS>(hds, edgeMap);
		fun.setInitialSolutionVec(init);
	}
	
	
	@Test
	public void testSinConditionFunctional() throws Exception {
		Vec init = fun.getSolutionVec();
		Vec startVec = new Vec(init.getSize());
		init.copy(startVec);
		
		Assert.assertTrue("energy is large at the beginning", 1E-3 < fun.evaluateObjective(init));
		
		Tao tao = new Tao(Method.CG);
		tao.setFromOptions();
		tao.setApplication(fun);
		tao.setMaximumIterates(20);
		tao.setGradientTolerances(1E-10, 1E-10, 1E-10);
		tao.setTolerances(1E-10, 1E-10, 1E-10, 1E-10);
		tao.solve();
		System.out.println(tao.getSolutionStatus());
		Assert.assertTrue("energy is small after minimizatin", 1E-9 > fun.evaluateObjective(init));
		
		double dif2 = 0.0;
		for (int i = 0; i < startVec.getSize(); i++) {
			double start = startVec.getValue(i);
			double solution = init.getValue(i);
			dif2 += (start-solution)*(start-solution);
		}
		Assert.assertEquals("solution start proximity", 1E-3, dif2, 1E-2);
	}
	
	
	@Test
	public void testGradient() throws Exception {
		Tao tao = Tao.createFiniteDifferenceTester(true, false);
		tao.setFromOptions();
		tao.setApplication(fun);
		tao.solve();
	}
	
	
	public static void main(String[] args) throws Exception {
		SinConditionFunctionalTest test = new SinConditionFunctionalTest();
		test.testEnergy();
	}
	
	public void testEnergy() throws Exception {
		createFunctional();
		Tao tao = Tao.createFiniteDifferenceTester(true, false);
		tao.setFromOptions();
		tao.setApplication(fun);
		tao.solve();
		
		int e1Index = 0;
		int e2Index = 1;
		double scale = 1E-2;
		
		double[] xArr = {0.5063494740335717, 0.2551584482990371, -0.24595351698479617, 0.06080650753547159, -0.38134954860685977, 0.36865045139194175};
		Vec x = fun.getSolutionVec();
		for (int i = 0; i < xArr.length; i++) {
			x.setValue(i, PI*xArr[i], INSERT_VALUES);
		}
		x.assemble();
		
		int uLine = 500;
		int vline = 500;
		
		double[][] verts = new double[uLine*vline][];
		Vec g = new Vec(xArr.length);
		g.assemble();
		
		for (int i = 0; i < uLine; i++) {
			for (int j = 0; j < vline; j++) {
				double x1 = PI * (i/(double)(uLine-1)) - PI/2;
				double x2 = PI * (j/(double)(vline-1)) - PI/2;
				x.setValue(e1Index, x1, INSERT_VALUES);
				x.setValue(e2Index, x2, INSERT_VALUES);
				double val = fun.evaluateObjective(x);
				fun.evaluateGradient(x, g);
//				fun.defaultComputeGradient(x, g);
				val = g.getValue(4);
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
		
		x.setValue(e1Index, PI*xArr[e1Index], INSERT_VALUES);
		x.setValue(e2Index, PI*xArr[e2Index], INSERT_VALUES);
		double val = fun.evaluateObjective(x);
		System.out.println("value: " + val);
		double[] markerPos = {PI*xArr[e1Index], PI*xArr[e2Index], scale * val};
		
		SceneGraphComponent root = new SceneGraphComponent();
		SceneGraphComponent graph = new SceneGraphComponent();
		SceneGraphComponent marker = new SceneGraphComponent();
		graph.setGeometry(qmf.getIndexedFaceSet());
		Appearance graphApp = new Appearance();
		graphApp.setAttribute(VERTEX_DRAW, false);
		graphApp.setAttribute(CommonAttributes.EDGE_DRAW, false);
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
