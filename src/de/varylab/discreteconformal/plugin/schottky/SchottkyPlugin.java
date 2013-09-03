package de.varylab.discreteconformal.plugin.schottky;

import static javax.swing.JOptionPane.WARNING_MESSAGE;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.AbstractTableModel;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.scripting.swing.ButtonCellEditor;
import de.jreality.plugin.scripting.swing.ButtonCellRenderer;
import de.jreality.ui.LayoutFactory;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.java2d.Viewer2D;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.image.ImageHook;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class SchottkyPlugin extends ShrinkPanelPlugin implements ActionListener {

	private HalfedgeInterface
		hif = null;
	private DiscreteConformalPlugin
		dcp = null;
	private SchottkyModeller
		schottkyModeller = new SchottkyModeller();
	private Viewer2D
		viewer = schottkyModeller.getViewer();
	private JButton
		generateButton = new JButton("Generate Surface"),
		toFuchsianButton = new JButton("Uniformize");
	private SpinnerNumberModel
		randomSeedModel = new SpinnerNumberModel(0, 0, 10000000, 1),
		equalizationIterationsModel = new SpinnerNumberModel(10, 0, 100, 1),
		cirleResModel = new SpinnerNumberModel(20, 4, 1000, 1),
		extraPointsModel = new SpinnerNumberModel(400, 0, 10000, 1);
	private JSpinner
		randomSeedSpinner = new JSpinner(randomSeedModel),
		equalizationIterationsSpinner = new JSpinner(equalizationIterationsModel),
		extraPointsSpinner = new JSpinner(extraPointsModel),
		circleResSpinner = new JSpinner(cirleResModel);
	private JCheckBox
		spherialChecker = new JCheckBox("Spherical"),
		cutFromRootChecker = new JCheckBox("Use Cut Root");
	private GeneratorModel
		generatorModel = new GeneratorModel();
	private JTable 
		generatorTable = new JTable(generatorModel);
	private JScrollPane
		generatorScroller = new JScrollPane(generatorTable);
	private Icon
		removeIcon = ImageHook.getIcon("remove.png");
	
	
	public SchottkyPlugin() {
		GridBagConstraints c1 = LayoutFactory.createLeftConstraint();
		GridBagConstraints c2 = LayoutFactory.createRightConstraint();
		shrinkPanel.setTitle("Schottky Modeller");
		shrinkPanel.setFillSpace(true);
		viewer.setPreferredSize(new Dimension(240, 240));
		viewer.setMinimumSize(viewer.getPreferredSize());
		JScrollPane scrollProtector = new JScrollPane(viewer);
		scrollProtector.setPreferredSize(viewer.getPreferredSize());
		scrollProtector.setMinimumSize(scrollProtector.getPreferredSize());
		scrollProtector.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		scrollProtector.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_NEVER);
		
		shrinkPanel.setLayout(new GridBagLayout());
		c2.weighty = 1.0;
		shrinkPanel.add(scrollProtector, c2);
		c2.weighty = 0.0;
		generatorScroller.setPreferredSize(new Dimension(10, 100));
		generatorTable.setFillsViewportHeight(true);
		generatorScroller.setBorder(BorderFactory.createTitledBorder("Generators"));
		generatorTable.setDefaultRenderer(JButton.class, new ButtonCellRenderer());
		generatorTable.setDefaultEditor(JButton.class, new ButtonCellEditor());
		generatorTable.getColumnModel().getColumn(0).setMaxWidth(40);
		generatorTable.getColumnModel().getColumn(2).setMaxWidth(40);
		generatorTable.setRowHeight(22);
		shrinkPanel.add(generatorScroller, c2);
		shrinkPanel.add(new JLabel("Random Seed"), c1);
		shrinkPanel.add(randomSeedSpinner, c2);
		shrinkPanel.add(new JLabel("Num Extra Points"), c1);
		shrinkPanel.add(extraPointsSpinner, c2);
		shrinkPanel.add(new JLabel("Point Equalizer Iterations"), c1);
		shrinkPanel.add(equalizationIterationsSpinner, c2);
		shrinkPanel.add(new JLabel("Circle Resolution"), c1);
		shrinkPanel.add(circleResSpinner, c2);
		shrinkPanel.add(cutFromRootChecker, c1);
		shrinkPanel.add(generateButton, c1);
		shrinkPanel.add(toFuchsianButton, c2);
		shrinkPanel.add(spherialChecker, c2);

		generateButton.addActionListener(this);
		toFuchsianButton.addActionListener(this);
		
		viewer.setScaleToolEnabled(true);
		viewer.setTranslateToolEnabled(true);
		viewer.setMenuToolEnabled(true);
		
		schottkyModeller.getModelContainer().addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				generatorTable.updateUI();
				schottkyModeller.getModeller().updateTools();
				schottkyModeller.getViewer().updateUI();
			}
		});
	}
	
	
	@Override
	public void actionPerformed(ActionEvent event) {
		List<SchottkyGenerator> pairs = schottkyModeller.getGenerators();
		Complex root = schottkyModeller.getBasePoint();
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		Set<Set<CoEdge>> cycles = new HashSet<Set<CoEdge>>();
		CoHDS hds = new CoHDS();
		Map<CoVertex, double[]> mapCycleMap = new HashMap<CoVertex, double[]>();
		
		int randomSeed = randomSeedModel.getNumber().intValue();
		int circleRes = cirleResModel.getNumber().intValue();
		int numExtraPoints = extraPointsModel.getNumber().intValue();
		int numSpreadIterations = equalizationIterationsModel.getNumber().intValue();
		CoVertex rootVertex = SchottkyUtility.generateSurface(hds, pairs, root, lMap, cycles, mapCycleMap, randomSeed, circleRes, numExtraPoints, numSpreadIterations);
		SchottkyLengthAdapter schottkyMetric = new SchottkyLengthAdapter(lMap);
		
		if (generateButton == event.getSource()) {
			hif.set(hds);
			hif.addLayerAdapter(schottkyMetric, false);
			return;
		}
		if (toFuchsianButton == event.getSource()) {
			AdapterSet aSet = hif.getAdapters();
			aSet.add(schottkyMetric);
			
			int genus = HalfEdgeUtils.getGenus(hds);
			System.out.println("unwrapping surface of genus " + genus + "...");
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = null;
			try {
				boolean onSphere = spherialChecker.isSelected();
				boolean cuFromRoot = cutFromRootChecker.isSelected();
				cutInfo = SchottkyUtility.unwrapSchottkySurface(hds, cycles, mapCycleMap, rootVertex, aSet, onSphere, cuFromRoot);
			} catch (Exception e) {
				JOptionPane.showMessageDialog(shrinkPanel, e.getMessage(), "Optimizer error", WARNING_MESSAGE);
				e.printStackTrace();
				return;
			}
			dcp.createVisualization(hds, genus, cutInfo);
			dcp.updateSurface();
			dcp.updateDomainImage();
		}
	}
	
	private class RemoveGeneratorButton extends JButton implements ActionListener {
		
		private static final long serialVersionUID = 1L;
		private SchottkyGenerator generator = null;

		public RemoveGeneratorButton(SchottkyGenerator generator) {
			super(removeIcon);
			this.generator = generator;
			addActionListener(this);
		}
		
		@Override
		public void actionPerformed(ActionEvent e) {
			schottkyModeller.removeGenerator(generator);
			schottkyModeller.getModeller().updateTools();
		}
		
	}
	

	private class GeneratorModel extends AbstractTableModel {

		private static final long serialVersionUID = 1L;

		@Override
		public int getColumnCount() {
			return 3;
		}

		@Override
		public int getRowCount() {
			return schottkyModeller.getGenerators().size();
		}
		
		@Override
		public boolean isCellEditable(int rowIndex, int columnIndex) {
			return columnIndex == 1 || columnIndex == 2;
		}
		
		@Override
		public String getColumnName(int col) {
			switch (col) {
			case 0: return "ID";
			case 1: return "|μ|";
			default: return "";
			}
		}
		
		@Override
		public Class<?> getColumnClass(int columnIndex) {
			switch (columnIndex) {
			case 0: return Integer.class;
			case 1: return Double.class;
			case 2: return JButton.class;
			default: return String.class;
			}
		}

		@Override
		public Object getValueAt(int row, int col) {
			SchottkyGenerator g = schottkyModeller.getGenerators().get(row);
			switch (col) {
			case 0: return col;
			case 1: return g.getMu().abs();
			case 2: return new RemoveGeneratorButton(g);
			default: return 0;
			}
		}
		
		@Override
		public void setValueAt(Object value, int row, int col) {
			SchottkyGenerator g = schottkyModeller.getGenerators().get(row);
			switch (col) {
			case 1: 
				double abs = (Double)value;
				double oldAbs = g.getMu().abs();
				Complex newMu = g.getMu().times(abs / oldAbs);
				g.setMu(newMu);
				schottkyModeller.getModeller().updateTools();
			}
		}
		
	}
	

	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		dcp = c.getPlugin(DiscreteConformalPlugin.class);
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}
	
	
	public static void main(String[] args) {
		NativePathUtility.set("native");
		JRViewer v = new JRViewer();
		v.setPropertiesFile("Schottky.jrw");
		v.setPropertiesResource(SchottkyPlugin.class, "Schottky.jrw");
		v.addContentUI();
		v.addBasicUI();
		v.registerPlugin(SchottkyPlugin.class);
		v.startup();
	}
	
}
