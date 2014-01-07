package de.varylab.discreteconformal.plugin;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.lang.annotation.Annotation;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.AbstractCellEditor;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;

import de.jreality.plugin.basic.View;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.image.ImageHook;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.conformallab.data.DataFactory;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.ConformalData;
import de.varylab.conformallab.data.types.ConformalDataList;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMetric;
import de.varylab.conformallab.data.types.HyperEllipticAlgebraicCurve;
import de.varylab.conformallab.data.types.SchottkyData;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.MappedEdgeLengthAdapter;
import de.varylab.discreteconformal.plugin.schottky.SchottkyPlugin;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class ConformalDataPlugin extends ShrinkPanelPlugin implements ActionListener {

	private static Logger
		log = Logger.getLogger(ConformalDataPlugin.class.getName());
	private List<ConformalData>
		data = new ArrayList<ConformalData>();
	
	private HalfedgeInterface
		hif = null;
	private DiscreteConformalPlugin
		discreteConformalPlugin = null;
	private SchottkyPlugin
		schottkyPlugin = null;
	private HyperellipticCurvePlugin
		hyperellipticCurvePlugin = null;
	
	private DataModel
		dataModel = new DataModel();
	private JTable
		dataTable = new JTable(dataModel);
	private JScrollPane
		dataScroller = new JScrollPane(dataTable);
	private JButton
		exportButton = new JButton("Export...", ImageHook.getIcon("disk.png")),
		importButton = new JButton("Import...", ImageHook.getIcon("folder.png")),
		clearButton = new JButton("Clear");
	private JFileChooser
		fileChooser = new JFileChooser();
	
	public ConformalDataPlugin() {
		shrinkPanel.setTitle("Conformal Data");
		setInitialPosition(SHRINKER_LEFT);
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.BOTH;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 1.0;
		c.weighty = 1.0;
		shrinkPanel.setLayout(new GridBagLayout());
		shrinkPanel.add(dataScroller, c);
		dataScroller.setPreferredSize(new Dimension(10, 150));
		dataTable.setRowHeight(24);
		dataTable.setDefaultRenderer(ConformalData.class, new DataCellRenderer());
		dataTable.getColumnModel().getColumn(0).setCellEditor(new DataNameEditor());
		dataTable.getColumnModel().getColumn(1).setMaxWidth(30);
		dataTable.getColumnModel().getColumn(1).setCellRenderer(new LoadButtonEditor());
		dataTable.getColumnModel().getColumn(1).setCellEditor(new LoadButtonEditor());
		dataTable.getColumnModel().getColumn(2).setMaxWidth(30);
		dataTable.getColumnModel().getColumn(2).setCellRenderer(new DeleteButtonEditor());
		dataTable.getColumnModel().getColumn(2).setCellEditor(new DeleteButtonEditor());
		dataTable.getTableHeader().setPreferredSize(new Dimension(0, 0));
		c.weighty = 0.0;
		c.gridwidth = GridBagConstraints.RELATIVE;
		shrinkPanel.add(exportButton, c);
		exportButton.addActionListener(this);
		c.gridwidth = GridBagConstraints.REMAINDER;
		shrinkPanel.add(importButton, c);
		importButton.addActionListener(this);
		shrinkPanel.add(clearButton, c);
		clearButton.addActionListener(this);
		
		fileChooser.setAcceptAllFileFilterUsed(true);
		fileChooser.setFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "Conformal Data (*.xml)";
			}
			@Override
			public boolean accept(File f) {
				return f.isDirectory() || f.getName().toLowerCase().endsWith(".xml");
			}
		});
		fileChooser.setDialogTitle("Export Conformal Data");
		fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if (exportButton == e.getSource()) {
			Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
			int result = fileChooser.showSaveDialog(w);
			if (result != JFileChooser.APPROVE_OPTION) {
				return;
			}
			File file = fileChooser.getSelectedFile();
			if (!file.getName().toLowerCase().endsWith(".xml")) {
				file = new File(file.getAbsolutePath() + ".xml");
			}
			if (file.exists()) {
				result = JOptionPane.showConfirmDialog(w, "The file " + file.getName() + " exists. \nDo you want to overwrite this file?", "Overwrite?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
				if (result != JOptionPane.OK_OPTION) {
					return;
				}
			}
			try {
				FileOutputStream fout = new FileOutputStream(file);
				ConformalDataList list = new ConformalDataList();
				for (ConformalData d : data) {
					list.getData().add(d);
				}
				DataFactory.writeConformalDataList(list, fout);
			} catch (Exception e1) {
				log.log(Level.SEVERE, "Could not export conformal data", e1);
			}
		}
		if (importButton == e.getSource()) {
			Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
			int result = fileChooser.showOpenDialog(w);
			if (result != JFileChooser.APPROVE_OPTION) {
				return;
			}
			File file = fileChooser.getSelectedFile();
			try {
				FileInputStream fin = new FileInputStream(file);
				ConformalDataList cdl = DataFactory.readConformalDataList(fin);
				clearData();
				for (ConformalData data : cdl.getData()) {
					addData(data);
				}
			} catch (Exception e1) {
				log.log(Level.SEVERE, "Could not import conformal data", e1);
			}	
		}
		if (clearButton == e.getSource()) {
			clearData();
		}
	}
	
	private class DataModel extends DefaultTableModel {

		private static final long serialVersionUID = 1L;

		@Override
		public int getRowCount() {
			return data.size();
		}
		@Override
		public int getColumnCount() {
			return 3;
		}
		@Override
		public Class<?> getColumnClass(int columnIndex) {
			return ConformalData.class;
		}
		@Override
		public Object getValueAt(int row, int column) {
			return data.get(row);
		}
		@Override
		public boolean isCellEditable(int rol, int column) {
			return true;
		}
		@Override
		public void setValueAt(Object arg0, int arg1, int arg2) {
		}
		
	}
	
	private class DataCellRenderer extends DefaultTableCellRenderer {
		
		private static final long serialVersionUID = 1L;

		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			Component c = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
			if (value instanceof ConformalData && c instanceof JLabel) {
				ConformalData data = (ConformalData)value;
				JLabel label = (JLabel)c;
				String text = data.getName();
				if (text == null) {
					text = value.getClass().getSimpleName();
				}
				label.setText(text);
			}
			return c;
		}
		
	}
	
	private class DataNameEditor extends DefaultCellEditor {
		
		private static final long serialVersionUID = 1L;
		private ConformalData
			editedData = null;
		
		public DataNameEditor() {
			super(new JTextField());
		}

		@Override
		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
			editedData = (ConformalData)value;
			return super.getTableCellEditorComponent(table, editedData.getName(), isSelected, row, column);
		}
		
		@Override
		public boolean stopCellEditing() {
			String text = (String)getCellEditorValue();
			editedData.setName(text);
			return super.stopCellEditing();
		}
		
	}
	
	
	private class DeleteButtonEditor extends AbstractCellEditor implements TableCellRenderer, TableCellEditor {

		private static final long 
			serialVersionUID = 1L;
		private JButton
			button = new JButton(ImageHook.getIcon("remove.png"));
		private ConformalData 
			value = null;
		
		public DeleteButtonEditor() {
			button.setToolTipText("Remove Data");
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			return button;
		}

		@Override
		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
			this.value = (ConformalData)value;
			JButton deleteButton = new JButton(ImageHook.getIcon("remove.png"));
			deleteButton.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					fireEditingStopped();
					removeData(DeleteButtonEditor.this.value);
				}
			});
			return deleteButton;
		}
		
		@Override
		public Object getCellEditorValue() {
			return value;
		}
		
	}
	
	
	private class LoadButtonEditor extends AbstractCellEditor implements TableCellRenderer, TableCellEditor {

		private static final long 
			serialVersionUID = 1L;
		private JButton
			button = new JButton(ImageHook.getIcon("cog_go.png"));
		private ConformalData 
			value = null;
		
		public LoadButtonEditor() {
			button.setToolTipText("Load into UI");
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			return button;
		}

		@Override
		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
			this.value = (ConformalData)value;
			JButton deleteButton = new JButton(ImageHook.getIcon("cog_go.png"));
			deleteButton.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					fireEditingStopped();
					loadDataIntoUI(LoadButtonEditor.this.value);
				}
			});
			return deleteButton;
		}
		
		@Override
		public Object getCellEditorValue() {
			return value;
		}
		
	}
	
	private void updateUI() {
		dataModel.fireTableDataChanged();
	}
	
	private void loadDataIntoUI(ConformalData data) {
		if (data instanceof DiscreteEmbedding) {
			DiscreteEmbedding de = (DiscreteEmbedding)data;
			int genus = DataUtility.calculateGenus(de);
			System.out.println("loading embedding of genus " + genus);
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
			CoHDS hds = DataUtility.toHDS(de, cutInfo);
			for (CoVertex v : hds.getVertices()) {
				System.arraycopy(v.P, 0, v.T, 0, 4);
			}
			discreteConformalPlugin.createUniformization(hds, genus, cutInfo);
			discreteConformalPlugin.updateGeometry();
			discreteConformalPlugin.updateDomainImage();
		}
		if (data instanceof DiscreteMetric) {
			DiscreteMetric dm = (DiscreteMetric)data;
			MappedEdgeLengthAdapter lMap = new MappedEdgeLengthAdapter(1000.0);
			CoHDS hds = DataUtility.toHalfedgeAndLengths(dm, lMap);
			hif.set(hds);
			hif.addAdapter(lMap, false);
		}
		if (data instanceof SchottkyData) {
			SchottkyData sd = (SchottkyData)data;
			schottkyPlugin.setSchottkyData(sd);
		}
		if (data instanceof HyperEllipticAlgebraicCurve) {
			HyperEllipticAlgebraicCurve curve = (HyperEllipticAlgebraicCurve)data;
			hyperellipticCurvePlugin.setCurve(DataUtility.toCurve(curve));
		}
	}

	public void addDiscreteMetric(String name, CoHDS surface, AdapterSet a) {
		DiscreteMetric dm = DataUtility.toDiscreteMetric(name, surface, a);
		addData(dm);
	}

	@Position
	@Position4d
	@TexturePosition
	@TexturePosition4d
	private static class DirectDataAdapter extends AbstractAdapter<double[]> {
		private boolean
			useTextureData = false;
		public DirectDataAdapter(boolean useTextureData) {
			super(double[].class, true, false);
			this.useTextureData = useTextureData;
		}
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return CoVertex.class.isAssignableFrom(nodeClass);
		}
		@Override
		public double getPriority() {
			return Double.MAX_VALUE;
		}
		@Override
		public <
			V extends Vertex<V, E, F>, 
			E extends Edge<V, E, F>, 
			F extends Face<V, E, F>
		> double[] getV(V v, AdapterSet a) {
			CoVertex cv = (CoVertex)v;
			if (useTextureData) {
				return cv.T;
			} else {
				return cv.P;
			}
		}
	}
	public void addDiscreteTextureEmbedding(CoHDS surface, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		AdapterSet aSet = new AdapterSet(new DirectDataAdapter(true));
		addDiscreteEmbedding("Texture Embedding", surface, aSet, TexturePosition4d.class, cutInfo);
	}
	public void addDiscretePositionEmbedding(CoHDS surface, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		AdapterSet aSet = new AdapterSet(new DirectDataAdapter(false));
		addDiscreteEmbedding("Position Embedding", surface, aSet, Position4d.class, cutInfo);
	}
	public void addDiscreteEmbedding(String name, CoHDS surface, AdapterSet a, Class<? extends Annotation> type, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		DiscreteEmbedding de = DataUtility.toDiscreteEmbedding(name, surface, a, type, cutInfo);
		addData(de);
	}
	public void addData(ConformalData data) {
		this.data.add(data);
		updateUI();
	}
	public void removeData(ConformalData data) {
		this.data.remove(data);
		updateUI();
	}
	public void clearData() {
		this.data.clear();
		updateUI();
	}
	public List<ConformalData> getData() {
		return Collections.unmodifiableList(data);
	}

	@Override
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "exportDir", fileChooser.getCurrentDirectory().getAbsolutePath());
	}
	
	@Override
	public void restoreStates(Controller c) throws Exception {
		super.restoreStates(c);
		fileChooser.setCurrentDirectory(new File(c.getProperty(getClass(), "exportDir", ".")));
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		discreteConformalPlugin = c.getPlugin(DiscreteConformalPlugin.class);
		schottkyPlugin = c.getPlugin(SchottkyPlugin.class);
		hyperellipticCurvePlugin = c.getPlugin(HyperellipticCurvePlugin.class);
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

}
