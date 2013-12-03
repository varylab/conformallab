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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.plugin.image.ImageHook;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.conformallab.data.DataFactory;
import de.varylab.conformallab.data.types.ConformalData;
import de.varylab.conformallab.data.types.ConformalDataList;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMetric;
import de.varylab.conformallab.data.types.EmbeddedTriangle;
import de.varylab.conformallab.data.types.EmbeddedVertex;
import de.varylab.conformallab.data.types.MetricEdge;
import de.varylab.conformallab.data.types.MetricTriangle;
import de.varylab.conformallab.data.types.SchottkyData;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class ConformalDataPlugin extends ShrinkPanelPlugin implements ActionListener {

	private static Logger
		log = Logger.getLogger(ConformalDataPlugin.class.getName());
	private List<ConformalData>
		data = new ArrayList<ConformalData>();
	
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
		dataTable.setDefaultRenderer(ConformalData.class, new DataCellRenderer());
//		dataTable.setRowHeight(23);
		dataTable.getColumnModel().getColumn(0).setCellEditor(new DataNameEditor());
		dataTable.getColumnModel().getColumn(1).setMaxWidth(30);
		dataTable.getColumnModel().getColumn(1).setCellRenderer(new LoadButtonEditor());
		dataTable.getColumnModel().getColumn(1).setCellEditor(new LoadButtonEditor());
		dataTable.getColumnModel().getColumn(2).setMaxWidth(30);
		dataTable.getColumnModel().getColumn(2).setCellRenderer(new DeleteButtonEditor());
		dataTable.getColumnModel().getColumn(2).setCellEditor(new DeleteButtonEditor());
		dataTable.getTableHeader().setPreferredSize(new Dimension(0, 0));
		c.weighty = 0.0;
		c.gridwidth = 1;
		shrinkPanel.add(clearButton, c);
		clearButton.addActionListener(this);
		c.gridwidth = GridBagConstraints.REMAINDER;
		shrinkPanel.add(exportButton, c);
		exportButton.addActionListener(this);
		c.gridwidth = GridBagConstraints.REMAINDER;
		shrinkPanel.add(importButton, c);
		importButton.addActionListener(this);
		
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
		System.out.println("ConformalDataPlugin.updateUI()");
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
	
	private void loadDataIntoUI(ConformalData data) {
		if (data instanceof DiscreteEmbedding) {
			System.out.println("loading as embedding...");
		}
		if (data instanceof DiscreteMetric) {
			System.out.println("loading as metric...");
		}
		if (data instanceof SchottkyData) {
			System.out.println("loading into Schottky editor...");
		}
	}
	
	
	public void addDiscreteMetric(String name, CoHDS surface, AdapterSet a) {
		DiscreteMetric dm = new DiscreteMetric();
		dm.setName(name);
		int edgeIndex = 0;
		Map<Integer, Integer> edgeIndexMap = new HashMap<Integer, Integer>();
		for (CoEdge e : surface.getPositiveEdges()) {
			MetricEdge me = new MetricEdge();
			me.setIndex(edgeIndex);
			me.setLength(a.get(Length.class, e, Double.class));
			dm.getEdges().add(me);
			edgeIndexMap.put(e.getIndex(), edgeIndex);
			edgeIndex++;
		}
		for (CoFace f : surface.getFaces()) {
			MetricTriangle mt = new MetricTriangle();
			CoEdge e1 = f.getBoundaryEdge();
			CoEdge e2 = e1.getNextEdge();
			CoEdge e3 = e2.getNextEdge();
			e1 = e1.isPositive() ? e1 : e1.getOppositeEdge(); 
			e2 = e2.isPositive() ? e2 : e2.getOppositeEdge(); 
			e3 = e3.isPositive() ? e3 : e3.getOppositeEdge();
			mt.setEdge1(edgeIndexMap.get(e1.getIndex()));
			mt.setEdge2(edgeIndexMap.get(e2.getIndex()));
			mt.setEdge3(edgeIndexMap.get(e3.getIndex()));
			dm.getTriangles().add(mt);
		}
		addData(dm);
	}
	
	
	public void addDiscreteTextureEmbedding(CoHDS surface, AdapterSet a) {
		addDiscreteEmbedding("Texture Embedding", surface, a, Position4d.class);
	}
	
	public void addDiscretePositionEmbedding(CoHDS surface, AdapterSet a) {
		addDiscreteEmbedding("Position Embedding", surface, a, TexturePosition4d.class);
	}
		
	public void addDiscreteEmbedding(String name, CoHDS surface, AdapterSet a, Class<? extends Annotation> type) {
		DiscreteEmbedding de = new DiscreteEmbedding();
		de.setName(name);
		for (CoVertex v : surface.getVertices()) {
			EmbeddedVertex ev = new EmbeddedVertex();
			double[] pos = a.getD(type, v);
			ev.setX(pos[0]);
			ev.setY(pos[1]);
			ev.setZ(pos[2]);
			ev.setW(pos[3]);
			ev.setIndex(v.getIndex());
			de.getVertices().add(ev);
		}
		for (CoFace f : surface.getFaces()) {
			List<CoEdge> b = HalfEdgeUtils.boundaryEdges(f);
			assert b.size() == 3 : "only triangulations are supported";
			CoVertex v0 = b.get(0).getStartVertex();
			CoVertex v1 = b.get(1).getStartVertex();
			CoVertex v2 = b.get(2).getStartVertex();
			EmbeddedTriangle et = new EmbeddedTriangle();
			et.setVertex1(v0.getIndex());
			et.setVertex2(v1.getIndex());
			et.setVertex3(v2.getIndex());
			de.getTriangles().add(et);
		}
		addData(de);
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
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

}
