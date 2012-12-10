package de.varylab.discreteconformal.plugin;

import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryVertices;
import static de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicLayout.doTexLayout;
import static de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility.calculateOrientationFromAlphas;
import static de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility.createAlphaField;
import static de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility.cutConesToBoundary;
import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;
import static java.awt.GridBagConstraints.HORIZONTAL;
import static java.awt.GridBagConstraints.RELATIVE;
import static java.awt.GridBagConstraints.REMAINDER;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static javax.swing.ListSelectionModel.SINGLE_SELECTION;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
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

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import de.jreality.math.Pn;
import de.jreality.plugin.basic.View;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.util.HalfEdgeUtils;
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
import de.varylab.discreteconformal.functional.EuclideanFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.generator.TestVectorFieldGenerator;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.UnwrapException;
import de.varylab.discreteconformal.unwrapper.circlepattern.CirclePatternLayout;
import de.varylab.discreteconformal.unwrapper.circlepattern.CirclePatternUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.discreteconformal.unwrapper.quasiisothermic.ConformalStructureUtility;
import de.varylab.discreteconformal.unwrapper.quasiisothermic.DBFSolution;
import de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicLayout;
import de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility;
import de.varylab.discreteconformal.unwrapper.quasiisothermic.SinConditionApplication;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.Search;
import de.varylab.discreteconformal.util.UnwrapUtility;

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
		goVectorfieldButton = new JButton("Calculate Flat Vectorfield"),
		goDBFSnesButton = new JButton("Calculate DBF SNES"),
		goConformalStructureButton = new JButton("Calculate via Conformal Structure");
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
		dbfPanel.add(goVectorfieldButton, c);
		dbfPanel.add(goConformalStructureButton, c);
//		dbfPanel.add(goDBFSnesButton, c);
		
		circlePatternPanel.setLayout(new GridBagLayout());
		circlePatternPanel.setBorder(BorderFactory.createTitledBorder("Circle Patterns"));
		circlePatternPanel.add(goCirclePatternButton, c);
		
		goCirclePatternButton.addActionListener(this);
		goDBFButton.addActionListener(this);
		goDBFSnesButton.addActionListener(this);
		goVectorfieldButton.addActionListener(this);
		goConformalStructureButton.addActionListener(this);
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
			double alpha = QuasiisothermicUtility.getSignedAngle(N, K, E);
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
				calculateWithSinFunctional(hds, edgeAlphaMap, a, false);
			}
			if (goCirclePatternButton == s) {
				calculateWithCirclePattern(hds, edgeAlphaMap, a);
			}
			if (goDBFSnesButton == s) {
				calculateWithSNES(hds, edgeAlphaMap, a);
			}
			if (goVectorfieldButton == s) {
				calculateWithSinFunctional(hds, edgeAlphaMap, a, true);
			}
			if (goConformalStructureButton == s) {
				calculateWithConformalStructure(hds, edgeAlphaMap, a);
			}
		} catch (Exception e) {
			Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
			JOptionPane.showMessageDialog(w, e.toString());
			e.printStackTrace();
		}
	}
	
	
	protected void calculateWithConformalStructure(CoHDS hds, Map<CoEdge, Double> initAlphas, AdapterSet a) {
		try {
			Map<CoVertex, Double> thetaMap = ConformalStructureUtility.boundaryAnglesFromAlphas(hds, initAlphas);
			Map<CoEdge, Double> lcrPseudoMap = ConformalStructureUtility.calculatePseudoConformalStructure(hds, initAlphas);
			Map<CoEdge, Double> lcrMap = ConformalStructureUtility.calculateConformalStructure(hds, lcrPseudoMap);
			Map<CoEdge, Double> lengthMap = ConformalStructureUtility.lengthsFromCrossRatios(hds, lcrMap);
			Map<CoVertex, double[]> pMap = ConformalStructureUtility.calculateFlatRepresentation(hds, lengthMap, thetaMap);
			for (CoVertex v : hds.getVertices()) {
				double[] T = pMap.get(v);
				Pn.dehomogenize(T, T);
				v.T = T;
			}
		} catch (UnwrapException e) {
			e.printStackTrace();
		}
		hif.update();
	}
	
	
	protected void calculateWithCirclePattern(CoHDS hds, Map<CoEdge, Double> initAlphas, AdapterSet a) {
		Map<CoEdge, Double> betaMap = QuasiisothermicUtility.calculateBetasFromAlphas(hds, initAlphas);
		QuasiisothermicUtility.createDelaunayAngleSystem(hds, betaMap);
		
		Map<CoEdge, Double> thetaMap = QuasiisothermicUtility.calculateThetasFromBetas(hds, betaMap);
		Map<CoFace, Double> phiMap = QuasiisothermicUtility.calculatePhisFromBetas(hds, betaMap);
		Map<CoFace, Double> rhoMap = CirclePatternUtility.calculateCirclePatternRhos(hds, thetaMap, phiMap);

		QuasiisothermicUtility.cutConesToBoundary(hds, betaMap);
		CirclePatternLayout.execute(hds, rhoMap, thetaMap, a, TexturePosition2d.class, TexturePosition.class);
		QuasiisothermicUtility.alignLayout(hds, initAlphas, a);
		
		hif.update();
		hif.addLayerAdapter(new MappedWeightAdapter(thetaMap, "Quasiisothermic Thetas"), false);
		hif.addLayerAdapter(new MappedWeightAdapter(betaMap, "Quasiisothermic Betas"), false);
		hif.addLayerAdapter(new MappedWeightAdapter(initAlphas, "Quasiisothermic Alphas"), false);
	}
	
	
	protected void calculateWithSinFunctional(CoHDS hds, Map<CoEdge, Double> initAlphas, AdapterSet a, boolean vectorsOnly) {
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
		if (vectorsOnly) {
			Adapter<double[]> aInit = createAlphaField(hds, initAlphas, a, "Quasiisothermic Initial Alpha");
			Adapter<double[]> aNew = createAlphaField(hds, alphaMap, a, "Quasiisothermic Flat Alpha");
			hif.addAdapter(aInit, false);
			hif.addAdapter(aNew, false);
			hif.update();
			return;
		}
		
		inspectAlphas(hds, alphaMap, a);
//		subdivideFaceSingularities(hds, alphaMap, a);
		
		Map<CoFace, Double> orientationMap = calculateOrientationFromAlphas(hds,alphaMap);

		
		cutManifoldToDisk(hds, hds.getVertex(0), null);
		cutMesh(hds, alphaMap, a);
		
		doTexLayout(hds, alphaMap, orientationMap, a);
		Adapter<double[]> aInit = QuasiisothermicUtility.createAlphaField(hds, initAlphas, a, "Quasiisothermic Initial Alpha");
		Adapter<double[]> aOpt = QuasiisothermicUtility.createAlphaField(hds, alphaMap, a, "Quasiisothermic Flat Alpha");
		hif.addAdapter(aInit, false);
		hif.addAdapter(aOpt, false);
		hif.update();
	}
	
	
	public static void cutMesh(CoHDS hds, Map<CoEdge, Double> alpha, AdapterSet a) {
		Map<CoEdge, Double> beta =  QuasiisothermicUtility.calculateBetasFromAlphas(hds, alpha);
		Set<CoVertex> cones = new HashSet<CoVertex>();
		for (CoVertex v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			double betaSum = QuasiisothermicUtility.calculateAngleSumFromBetas(v, beta);
			double index = betaSum / PI;
			if (abs(index - 2) > 0.1) {
				cones.add(v);
			}
		}
		Set<CoEdge> validSet = new HashSet<CoEdge>(hds.getEdges());
		
		CuttingInfo<CoVertex, CoEdge, CoFace> info = new CuttingInfo<CoVertex, CoEdge, CoFace>();
		for (CoVertex c : cones) {
			if (HalfEdgeUtils.isBoundaryVertex(c)) {
				System.out.println("vertex " + c + " already on boundary");
				continue;
			}
			System.out.println("cutting from cone " + c);
			List<CoEdge> path = Search.bFS(validSet, c, boundaryVertices(hds), true, null);
			CuttingUtility.cutAlongPath(path, info);
		}
	}
	
	
	
	public void inspectAlphas(CoHDS hds, Map<CoEdge, Double> alpha, AdapterSet a) {
		Map<CoEdge, Double> beta =  QuasiisothermicUtility.calculateBetasFromAlphas(hds, alpha);
		System.out.println("alpha inspection -------------");
		Map<CoFace, Double> orientationMap = calculateOrientationFromAlphas(hds, alpha);
		for (CoFace f : hds.getFaces()) {
			Double orientation = orientationMap.get(f);
			if (orientation < 0) {
				List<CoVertex> bv = HalfEdgeUtils.boundaryVertices(f);
				System.out.println(f + " -> flipped: " + bv);
			}
		}
		System.out.println("singularities -----");
		for (CoVertex v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			double index = QuasiisothermicUtility.getSingularityIndexV(v, alpha, a);
			if (Math.abs(index) > 0.25) {
				System.out.println(v + " -> index = " + index);
			}
		}
		for (CoFace f : hds.getFaces()) {
			double index = QuasiisothermicUtility.getSingularityIndexF(f, alpha, a);
			if (Math.abs(index) > 0.25) {
				System.out.println(f + " -> index = " + index);
			}
		}
		System.out.println("angle sums (beta) --");
		for (CoVertex v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			double betaSum = QuasiisothermicUtility.calculateAngleSumFromBetas(v, beta);
			double index = betaSum / PI;
			if (abs(index - 2) > 0.1) {
				System.out.println(v + " -> sum / PI = " + index);
			}
		}
		System.out.println("------------------------------");
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
		Map<CoFace, Double> orientationMap = QuasiisothermicUtility.calculateOrientationFromAlphas(hds,alphaMap);
		Map<CoEdge, Double> betaMap = QuasiisothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		
		cutManifoldToDisk(hds, hds.getVertex(0), null);
		cutConesToBoundary(hds, betaMap);
		
		QuasiisothermicLayout.doTexLayout(hds, alphaMap, orientationMap, a);
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
		c.getPlugin(ProjectiveTexturePlugin.class);
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
