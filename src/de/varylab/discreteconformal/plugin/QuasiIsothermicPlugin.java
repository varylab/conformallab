package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;
import static java.awt.GridBagConstraints.HORIZONTAL;
import static java.awt.GridBagConstraints.RELATIVE;
import static java.awt.GridBagConstraints.REMAINDER;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;

import de.jreality.plugin.basic.View;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.jtao.Tao.Method;
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

	private Method[]
		methods = {Method.LMVM, Method.CG, Method.NTR, Method.NLS};
	private HalfedgeInterface
		hif = null;
	private JPanel
		circlePatternPanel = new JPanel(),
		dbfPanel = new JPanel();
	private SpinnerNumberModel
		maxItModel = new SpinnerNumberModel(200, 1, 10000, 1),
		tolExpModel = new SpinnerNumberModel(-6, -25, -1, -1);
	private JSpinner
		maxItSpinner = new JSpinner(maxItModel),
		tolExpSpinner = new JSpinner(tolExpModel);
	private JCheckBox
		excludeBoundaryChecker = new JCheckBox("Exclude Boundary", false);
	private JComboBox
		methodCombo = new JComboBox(methods);
	private JButton
		goCirclePatternButton = new JButton("Calculate Circle Pattern"),
		goDBFButton = new JButton("Calculate DBF");
	
	public QuasiIsothermicPlugin() {
		shrinkPanel.setTitle("Quasiisothermic Parametrization");
		shrinkPanel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.fill = HORIZONTAL;
		c.weightx = 1.0;
		c.insets = new Insets(2, 2, 2, 2);
		c.gridwidth = REMAINDER;
		shrinkPanel.add(dbfPanel, c);
		shrinkPanel.add(circlePatternPanel, c);
		
		dbfPanel.setLayout(new GridBagLayout());
		dbfPanel.setBorder(BorderFactory.createTitledBorder("DBF"));
		dbfPanel.add(excludeBoundaryChecker, c);
		c.gridwidth = RELATIVE;
		dbfPanel.add(new JLabel("Method"), c);
		c.gridwidth = REMAINDER;
		dbfPanel.add(methodCombo, c);
		c.gridwidth = RELATIVE;
		dbfPanel.add(new JLabel("Max Iterations"), c);
		c.gridwidth = REMAINDER;
		dbfPanel.add(maxItSpinner, c);
		c.gridwidth = RELATIVE;
		dbfPanel.add(new JLabel("Tolerance"), c);
		c.gridwidth = REMAINDER;
		dbfPanel.add(tolExpSpinner, c);
		
		dbfPanel.add(goDBFButton, c);
		
		circlePatternPanel.setLayout(new GridBagLayout());
		circlePatternPanel.setBorder(BorderFactory.createTitledBorder("Circle Patterns"));
		circlePatternPanel.add(goCirclePatternButton, c);
		
		goCirclePatternButton.addActionListener(this);
		goDBFButton.addActionListener(this);
	}
	
	@Override
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "excludeBoundary", excludeBoundaryChecker.isSelected());
		c.storeProperty(getClass(), "dbfMethod", methodCombo.getSelectedIndex());
		c.storeProperty(getClass(), "maxIter", maxItModel.getNumber().intValue());
		c.storeProperty(getClass(), "tolExp", tolExpModel.getNumber().intValue());
	}

	@Override
	public void restoreStates(Controller c) throws Exception {
		super.restoreStates(c);
		excludeBoundaryChecker.setSelected(c.getProperty(getClass(), "excludeBoundary", false));
		methodCombo.setSelectedIndex(c.getProperty(getClass(), "dbfMethod", 0));
		maxItModel.setValue(c.getProperty(getClass(), "maxIter", 200));
		tolExpModel.setValue(c.getProperty(getClass(), "tolExp", -6));
	}
	
	protected Method getSelectedMethod() {
		return (Method)methodCombo.getSelectedItem();
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
			e.printStackTrace();
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
		
		int maxIt = maxItModel.getNumber().intValue();
		int tolExp = tolExpModel.getNumber().intValue();
		double tol = Math.pow(10, tolExp);
		
		Method method = getSelectedMethod();
		fun.solveEnergyMinimzation(maxIt, tol, method);
		DBFSolution<CoVertex, CoEdge, CoFace, CoHDS> solution = fun.getDBFSolution();
		Map<CoEdge, Double> alphaMap = solution.solutionAlphaMap;
		Map<CoFace, Double> orientationMap = IsothermicUtility.calculateOrientationFromAlphas(hds,alphaMap);
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		
		
		// remove topology
		if (HalfEdgeUtils.getGenus(hds) >= 1) {
			CoVertex cutRoot = hds.getVertex(0);
			cutManifoldToDisk(hds, cutRoot, null);
		}
		
		IsothermicUtility.cutConesToBoundary(hds, betaMap);
		
		IsothermicLayout.doTexLayout(hds, alphaMap, orientationMap, a);
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
