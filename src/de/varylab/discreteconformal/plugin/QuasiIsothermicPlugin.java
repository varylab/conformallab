package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.unwrapper.isothermic.IsothermicLayout.doTexLayout;
import static de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility.cutConesToBoundary;
import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;
import static java.awt.GridBagConstraints.HORIZONTAL;
import static java.awt.GridBagConstraints.RELATIVE;
import static java.awt.GridBagConstraints.REMAINDER;
import static javax.swing.ListSelectionModel.SINGLE_SELECTION;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.table.DefaultTableModel;

import de.jreality.plugin.basic.View;
import de.jtem.halfedge.Edge;
import de.jtem.halfedgetools.adapter.Adapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMax;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.VectorField;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition2d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeListener;
import de.jtem.halfedgetools.plugin.data.AdapterNameComparator;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.Method;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.generator.TestVectorFieldGenerator;
import de.varylab.discreteconformal.unwrapper.circlepattern.CirclePatternLayout;
import de.varylab.discreteconformal.unwrapper.circlepattern.CirclePatternUtility;
import de.varylab.discreteconformal.unwrapper.isothermic.DBFSolution;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility;
import de.varylab.discreteconformal.unwrapper.isothermic.SinConditionApplication;

public class QuasiIsothermicPlugin extends ShrinkPanelPlugin implements ActionListener, HalfedgeListener {

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
		goDBFButton = new JButton("Calculate DBF Variational"),
		goDBFSnesButton = new JButton("Calculate DBF SNES");
	private JTable
		dataTable = new JTable();
	private JScrollPane
		dataScroller = new JScrollPane(dataTable);
	private Set<Adapter<?>>
		vecSet = new TreeSet<Adapter<?>>(new AdapterNameComparator());
	
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
		dataScroller.setPreferredSize(new Dimension(10, 75));
		dataScroller.setMinimumSize(new Dimension(10, 75));
		dataTable.getTableHeader().setPreferredSize(new Dimension(10, 0));
		dataTable.setSelectionMode(SINGLE_SELECTION);
		dbfPanel.add(dataScroller, c);
		
		dbfPanel.add(goDBFButton, c);
		dbfPanel.add(goDBFSnesButton, c);
		
		circlePatternPanel.setLayout(new GridBagLayout());
		circlePatternPanel.setBorder(BorderFactory.createTitledBorder("Circle Patterns"));
		circlePatternPanel.add(goCirclePatternButton, c);
		
		goCirclePatternButton.addActionListener(this);
		goDBFButton.addActionListener(this);
		goDBFSnesButton.addActionListener(this);
	}
	
	protected class VecTableModel extends DefaultTableModel {
		
		private static final long 
			serialVersionUID = 1L;

		@Override
		public int getColumnCount() {
			return 1;
		}
		
		@Override
		public Class<?> getColumnClass(int columnIndex) {
			switch (columnIndex) {
				default: return String.class;
			}
		}
		
		@Override
		public int getRowCount() {
			return vecSet.size();
		}
		
		@Override
		public Object getValueAt(int row, int column) {
			if (row < 0 || row >= vecSet.size()) {
				return "-";
			}
			Object[] objects = vecSet.toArray();
			Object op = objects[row];
			return op;
		}
		
		@Override
		public boolean isCellEditable(int row, int column) {
			return false;
		}
		
	}
	
	public Adapter<double[]> getSetectedData() {
		int row = dataTable.getSelectedRow();
		if (row < 0) {
			return null;
		}
		@SuppressWarnings("unchecked")
		Adapter<double[]> data = (Adapter<double[]>)dataTable.getModel().getValueAt(row, 0);
		return data;
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
		Tao.Initialize();
		// create data structure independent data
		AdapterSet a = hif.getAdapters();
		Map<Integer, Double> indexAlphaMap = new HashMap<Integer, Double>();
		Adapter<?> data = getSetectedData();
		if (data == null) {
			Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
			JOptionPane.showMessageDialog(w, "Please generate and select a vector field");
		}
		for (Edge<?,?,?> e : hif.get().getPositiveEdges()) {
			Object kData = data.get(e, hif.getAdapters());
			double[] K = null;
			if (kData instanceof double[]) {
				K = (double[])kData;
			} else 
			if (kData instanceof double[][]) {
				K = ((double[][])kData)[0];
			} else {
				throw new RuntimeException("cannot work with VectorField data of type " + kData.getClass().getName());
			}
			
			double[] N = a.getD(Normal.class, e);
			double[] E = a.getD(EdgeVector.class, e);
			if (N == null || K == null || E == null) {
				throw new RuntimeException("Could not get curvature information at edge " + e);
			}
			double alpha = IsothermicUtility.getSignedAngle(N, K, E);
			indexAlphaMap.put(e.getIndex(), alpha);
			indexAlphaMap.put(e.getOppositeEdge().getIndex(), alpha);
		}
		
		CoHDS hds = hif.get(new CoHDS());
		Map<CoEdge, Double> edgeAlphaMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getEdges()) {
			Double alpha = indexAlphaMap.get(e.getIndex()); 
			edgeAlphaMap.put(e, alpha);
			edgeAlphaMap.put(e.getOppositeEdge(), alpha);
		}
		
		Object s = ae.getSource();
		try {
			if (goDBFButton == s) {
				calculateWithSinFunctional(hds, edgeAlphaMap, a);
			}
			if (goCirclePatternButton == s) {
				calculateWithCirclePattern(hds, edgeAlphaMap, a);
			}
			if (goDBFSnesButton == s) {
				calculateWithSNES(hds, edgeAlphaMap, a);
			}
		} catch (Exception e) {
			Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
			JOptionPane.showMessageDialog(w, e.toString());
			e.printStackTrace();
		}
	}
	
	
	protected void calculateWithCirclePattern(CoHDS hds, Map<CoEdge, Double> initAlphas, AdapterSet a) {
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, initAlphas);
		IsothermicUtility.createDelaunayAngleSystem(hds, betaMap);
		
		Map<CoEdge, Double> thetaMap = IsothermicUtility.calculateThetasFromBetas(hds, betaMap);
		Map<CoFace, Double> phiMap = IsothermicUtility.calculatePhisFromBetas(hds, betaMap);
		Map<CoFace, Double> rhoMap = CirclePatternUtility.calculateCirclePatternRhos(hds, thetaMap, phiMap);

		IsothermicUtility.cutConesToBoundary(hds, betaMap);
		CirclePatternLayout.execute(hds, rhoMap, thetaMap, a, TexturePosition2d.class, TexturePosition.class);
		IsothermicUtility.alignLayout(hds, initAlphas);
		
		hif.update();
		hif.addLayerAdapter(new MappedWeightAdapter(thetaMap, "Quasiisothermic Thetas"), false);
		hif.addLayerAdapter(new MappedWeightAdapter(betaMap, "Quasiisothermic Betas"), false);
		hif.addLayerAdapter(new MappedWeightAdapter(initAlphas, "Quasiisothermic Alphas"), false);
	}
	
	
	protected void calculateWithSinFunctional(CoHDS hds, Map<CoEdge, Double> initAlphas, AdapterSet a) {
		boolean excludeBoundary = excludeBoundaryChecker.isSelected();
		
		SinConditionApplication<CoVertex, CoEdge, CoFace, CoHDS> 
		fun = new SinConditionApplication<CoVertex, CoEdge, CoFace, CoHDS>(hds);
		fun.initialize(initAlphas, excludeBoundary);
		
		int maxIt = maxItModel.getNumber().intValue();
		int tolExp = tolExpModel.getNumber().intValue();
		double tol = Math.pow(10, tolExp);
		
		Method method = getSelectedMethod();
		fun.solveEnergyMinimzation(maxIt, tol, method);
		DBFSolution<CoVertex, CoEdge, CoFace, CoHDS> solution = fun.getDBFSolution();
		
		Map<CoEdge, Double> alphaMap = solution.solutionAlphaMap;
		Map<CoFace, Double> orientationMap = IsothermicUtility.calculateOrientationFromAlphas(hds,alphaMap);
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		
		cutManifoldToDisk(hds, hds.getVertex(0), null);
		cutConesToBoundary(hds, betaMap);
		
		doTexLayout(hds, alphaMap, orientationMap, a);
		
		hif.update();
	}
	
	protected void calculateWithSNES(CoHDS hds, Map<CoEdge, Double> initAlphas, AdapterSet a) {
		boolean excludeBoundary = excludeBoundaryChecker.isSelected();
		
		SinConditionApplication<CoVertex, CoEdge, CoFace, CoHDS> 
		fun = new SinConditionApplication<CoVertex, CoEdge, CoFace, CoHDS>(hds);
		fun.initialize(initAlphas, excludeBoundary);
		
		int maxIt = maxItModel.getNumber().intValue();
		int tolExp = tolExpModel.getNumber().intValue();
		double tol = Math.pow(10, tolExp);
		
		fun.solveSNES(maxIt, tol);
		DBFSolution<CoVertex, CoEdge, CoFace, CoHDS> solution = fun.getDBFSolution();
		
		Map<CoEdge, Double> alphaMap = solution.solutionAlphaMap;
		Map<CoFace, Double> orientationMap = IsothermicUtility.calculateOrientationFromAlphas(hds,alphaMap);
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		
		cutManifoldToDisk(hds, hds.getVertex(0), null);
		cutConesToBoundary(hds, betaMap);
		
		doTexLayout(hds, alphaMap, orientationMap, a);
		
		hif.update();
	}


	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		hif.addHalfedgeListener(this);
		c.getPlugin(TestVectorFieldGenerator.class);
		dataTable.setModel(new VecTableModel());
		c.getPlugin(DiscreteConformalPlugin.class);
	}
	
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

	@Override
	public void dataChanged(HalfedgeLayer layer) {
	}
	@Override
	public void adaptersChanged(HalfedgeLayer layer) {
		vecSet.clear();
		vecSet.addAll(hif.getAdapters().queryAll(VectorField.class));
		vecSet.addAll(hif.getAdapters().queryAll(VectorField.class));
		vecSet.addAll(hif.getAdapters().queryAll(CurvatureFieldMin.class));
		vecSet.addAll(hif.getAdapters().queryAll(CurvatureFieldMax.class));
		dataTable.setModel(new VecTableModel());
		dataTable.updateUI();
		dataTable.getSelectionModel().setSelectionInterval(0, 0);
	}
	@Override
	public void activeLayerChanged(HalfedgeLayer old, HalfedgeLayer active) {
	}
	@Override
	public void layerCreated(HalfedgeLayer layer) {
	}
	@Override
	public void layerRemoved(HalfedgeLayer layer) {
	}

}
