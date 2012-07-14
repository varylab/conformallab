package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;
import static java.awt.GridBagConstraints.HORIZONTAL;
import static java.awt.GridBagConstraints.REMAINDER;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import de.jreality.plugin.basic.View;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.generator.TestVectorFieldGenerator;
import de.varylab.discreteconformal.unwrapper.isothermic.DBFSolution;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicLayout;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility;
import de.varylab.discreteconformal.unwrapper.isothermic.SinConditionApplication;

public class QuasiIsothermicPlugin extends ShrinkPanelPlugin implements ActionListener {

	private HalfedgeInterface
		hif = null;
	private JPanel
		circlePatternPanel = new JPanel(),
		dbfPanel = new JPanel();
	private JCheckBox
		excludeBoundaryChecker = new JCheckBox("Exclude Boundary", true);
	private JButton
		goCirclePatternButton = new JButton("Calculate Circle Pattern"),
		goDBFButton = new JButton("Calculate DBF");
	
	public QuasiIsothermicPlugin() {
		shrinkPanel.setTitle("Quasiisothermic Parametrization");
		shrinkPanel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.fill = HORIZONTAL;
		c.weightx = 1.0;
		c.gridwidth = REMAINDER;
		shrinkPanel.add(dbfPanel, c);
		shrinkPanel.add(circlePatternPanel, c);
		
		dbfPanel.setLayout(new GridBagLayout());
		dbfPanel.setBorder(BorderFactory.createTitledBorder("DBF"));
		dbfPanel.add(excludeBoundaryChecker, c);
		dbfPanel.add(goDBFButton, c);
		
		circlePatternPanel.setLayout(new GridBagLayout());
		circlePatternPanel.setBorder(BorderFactory.createTitledBorder("Circle Patterns"));
		circlePatternPanel.add(goCirclePatternButton, c);
		
		goCirclePatternButton.addActionListener(this);
		goDBFButton.addActionListener(this);
	}
	
	@Override
	public void actionPerformed(ActionEvent ae) {
		Object s = ae.getSource();
		try {
			if (goDBFButton == s) {
				calculateWithSinFunctional();
			}
			if (goCirclePatternButton == s) {
				calculateWithCirclePattern();
			}
		} catch (Exception e) {
			Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
			JOptionPane.showMessageDialog(w, e.toString());
		}
	}
	
	
	protected void calculateWithCirclePattern() {
		AdapterSet a = hif.getAdapters();
		CoHDS hds = hif.get(new CoHDS());
		
		Map<CoEdge, Double> alphaMap = IsothermicUtility.calculateAlphasFromCurvature(a, hds);
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		
//		IsothermicUtility.checkTriangleAngles(hds, betaMap);
		IsothermicUtility.createDelaunayAngleSystem(hds, betaMap);
//		IsothermicUtility.checkTriangleAngles(hds, betaMap);
		
		Map<CoEdge, Double> thetaMap = IsothermicUtility.calculateThetasFromBetas(hds, betaMap);
		Map<CoFace, Double> phiMap = IsothermicUtility.calculatePhisFromBetas(hds, betaMap);
		Map<CoFace, Double> rhoMap = IsothermicUtility.calculateCirclePatternRhos(hds, thetaMap, phiMap);

		IsothermicUtility.cutConesToBoundary(hds, betaMap);
		
		IsothermicUtility.doCirclePatternLayout(hds, thetaMap, rhoMap);
		
		IsothermicUtility.alignLayout(hds, alphaMap);
		
		
		hif.update();
		hif.addLayerAdapter(new MappedWeightAdapter(thetaMap, "Quasiisothermic Thetas"), false);
		hif.addLayerAdapter(new MappedWeightAdapter(betaMap, "Quasiisothermic Betas"), false);
		hif.addLayerAdapter(new MappedWeightAdapter(alphaMap, "Quasiisothermic Alphas"), false);
	}
	
	
	protected void calculateWithSinFunctional() {
		System.out.println("using sin-functional");
		AdapterSet a = hif.getAdapters();
		CoHDS hds = hif.get(new CoHDS());
		
		boolean excludeBoundary = excludeBoundaryChecker.isSelected();
		
		SinConditionApplication<CoVertex, CoEdge, CoFace, CoHDS> 
		fun = new SinConditionApplication<CoVertex, CoEdge, CoFace, CoHDS>(hds);
		fun.initialize(a, excludeBoundary);
		fun.solveCG(1000, 1E-10);
		
		DBFSolution<CoVertex, CoEdge, CoFace, CoHDS> solution = fun.getDBFSolution();
		Map<CoEdge, Double> alphaMap = solution.solutionAlphaMap;
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		
		// remove topology
		if (HalfEdgeUtils.getGenus(hds) >= 1) {
			CoVertex cutRoot = hds.getVertex(0);
			cutManifoldToDisk(hds, cutRoot, null);
		}
		
		IsothermicUtility.cutConesToBoundary(hds, betaMap);
		
		IsothermicLayout.doTexLayout(hds, alphaMap, a);
		hif.update();
	}

	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		c.getPlugin(TestVectorFieldGenerator.class);
	}
	
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

}
