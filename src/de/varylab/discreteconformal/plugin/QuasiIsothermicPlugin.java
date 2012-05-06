package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;
import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.swing.JButton;

import de.jreality.plugin.basic.View;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.jpetsc.Vec;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.Method;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.ConesUtility;
import de.varylab.discreteconformal.unwrapper.IsothermicUtility;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicLayout;
import de.varylab.discreteconformal.unwrapper.isothermic.SinConditionFunctional;

public class QuasiIsothermicPlugin extends ShrinkPanelPlugin implements ActionListener {

	private HalfedgeInterface
		hif = null;
	private JButton
		goButton = new JButton("Go");
	
//	static {
//		String[] taoCommand = new String[] {
//			"-tao_nm_lamda", "0.01", 
//			"-tao_nm_mu", "1.0"
//		};
//		System.out.println("initing tao: " + Arrays.toString(taoCommand));
//		Tao.Initialize("Quasiisothermic Parametrization", taoCommand, false);
//	}
	
	public QuasiIsothermicPlugin() {
		shrinkPanel.setTitle("Quasiisothermic Parametrization");
		shrinkPanel.add(goButton);
		goButton.addActionListener(this);
	}
	
	@Override
	public void actionPerformed(ActionEvent ae) {
		AdapterSet a = hif.getAdapters();
		CoHDS hds = hif.get(new CoHDS());
		
		Map<Integer, Integer> edgeMap = IsothermicUtility.createUndirectedEdgeMap(hds);
		
		SinConditionFunctional<CoVertex, CoEdge, CoFace, CoHDS> 
		fun = new SinConditionFunctional<CoVertex, CoEdge, CoFace, CoHDS>(hds, edgeMap);
		fun.calculateAndSetInitionSolution(a);
		
		System.out.println("energy before optimization: " + fun.evaluateObjective(fun.getSolutionVec()));
		
		Tao tao = new Tao(Method.LMVM);
		tao.setFromOptions();
		tao.setApplication(fun);
		tao.setMaximumIterates(200);
		tao.setMaximumFunctionEvaluations(1000000);
		tao.setTolerances(1E-10, 1E-10, 1E-10, 1E-10);
		tao.setGradientTolerances(1E-10, 1E-10, 1E-10);
		tao.solve();
		System.out.println(tao.getSolutionStatus());
		System.out.println("energy after optimization: " + fun.evaluateObjective(fun.getSolutionVec()));
		
		Vec solution = fun.getSolutionVec();
		Map<CoEdge, Double> alphaMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getEdges()) {
			int solverIndex = edgeMap.get(e.getIndex());
			double alpha = solution.getValue(solverIndex);
			alphaMap.put(e, alpha);
			alphaMap.put(e.getOppositeEdge(), alpha);
		}
		
		// remove topology
		if (HalfEdgeUtils.getGenus(hds) >= 1) {
			CoVertex cutRoot = hds.getVertex(0);
			cutManifoldToDisk(hds, cutRoot, null);
		}
		
		Set<CoVertex> innerVerts = new HashSet<CoVertex>(hds.getVertices());
		innerVerts.removeAll(HalfEdgeUtils.boundaryVertices(hds));
		for (CoVertex v : innerVerts) {
			double sum = IsothermicLayout.calculateAngleSum(v, alphaMap);
			if (abs(sum - 2*PI) > Math.PI/4) {
				int index = (int)Math.round(sum / PI);
				v.setTheta(index * PI);
				System.out.println("singularity: " + v + ", " + index + " pi");
			} else {
				v.setTheta(2 * PI);
			}
		}
		ConesUtility.cutMesh(hds);
		
		IsothermicLayout.doTexLayout(hds, alphaMap, a);
		hif.update();
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
	}
	
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

}
