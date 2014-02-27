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
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.image.ImageHook;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.ConformalData;
import de.varylab.conformallab.data.types.ConformalDataList;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMap;
import de.varylab.conformallab.data.types.DiscreteMetric;
import de.varylab.conformallab.data.types.EmbeddedVertex;
import de.varylab.conformallab.data.types.HyperEllipticAlgebraicCurve;
import de.varylab.conformallab.data.types.SchottkyData;
import de.varylab.conformallab.data.types.UniformizationData;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoDirectDataAdapter;
import de.varylab.discreteconformal.heds.adapter.CoDirectPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoDirectTextureAdapter;
import de.varylab.discreteconformal.heds.adapter.MappedEdgeLengthAdapter;
import de.varylab.discreteconformal.plugin.schottky.SchottkyPlugin;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
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
		dataTable.getColumnModel().getColumn(0).setMaxWidth(30);
		dataTable.getColumnModel().getColumn(1).setCellEditor(new DataNameEditor());
		dataTable.getColumnModel().getColumn(2).setMaxWidth(30);
		dataTable.getColumnModel().getColumn(2).setCellRenderer(new ExportButtonEditor());
		dataTable.getColumnModel().getColumn(2).setCellEditor(new ExportButtonEditor());
		dataTable.getColumnModel().getColumn(3).setMaxWidth(30);
		dataTable.getColumnModel().getColumn(3).setCellRenderer(new LoadButtonEditor());
		dataTable.getColumnModel().getColumn(3).setCellEditor(new LoadButtonEditor());
		dataTable.getColumnModel().getColumn(4).setMaxWidth(30);
		dataTable.getColumnModel().getColumn(4).setCellRenderer(new DeleteButtonEditor());
		dataTable.getColumnModel().getColumn(4).setCellEditor(new DeleteButtonEditor());
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
			File file = showXMLSaveDialog();
			if (file == null) return;
			try {
				FileOutputStream fout = new FileOutputStream(file);
				ConformalDataList list = new ConformalDataList();
				for (ConformalData d : data) {
					list.getData().add(d);
				}
				DataIO.writeConformalDataList(list, fout);
				fout.close();
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
				ConformalDataList cdl = DataIO.readConformalDataList(fin);
				clearData();
				for (ConformalData data : cdl.getData()) {
					addData(data);
				}
				fin.close();
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
			return 5;
		}
		@Override
		public Class<?> getColumnClass(int columnIndex) {
			switch (columnIndex) {
			case 0:
				return Integer.class;
			default:
				return ConformalData.class;
			}
		}
		@Override
		public Object getValueAt(int row, int column) {
			switch (column) {
			case 0:
				return row;
			default:
				return data.get(row);
			}
		}
		@Override
		public boolean isCellEditable(int rol, int column) {
			switch (column) {
			case 0:
				return false;
			default:
				return true;
			}
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
	
	private class ExportButtonEditor extends AbstractCellEditor implements TableCellRenderer, TableCellEditor {

		private static final long 
			serialVersionUID = 1L;
		private JButton
			button = new JButton(ImageHook.getIcon("disk.png"));
		private ConformalData 
			value = null;
		
		public ExportButtonEditor() {
			button.setToolTipText("Export Data");
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			return button;
		}

		@Override
		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
			this.value = (ConformalData)value;
			JButton saveButton = new JButton(ImageHook.getIcon("disk.png"));
			saveButton.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					fireEditingStopped();
					File file = showXMLSaveDialog();
					if (file == null) return;
					try {
						FileOutputStream fout = new FileOutputStream(file);
						DataIO.writeConformalData(ExportButtonEditor.this.value, fout);
						fout.close();
					} catch (Exception e1) {
						log.log(Level.SEVERE, "Could not export conformal data", e1);
					}
				}
			});
			return saveButton;
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
		
		private void checkButtonStatus(JButton button, Object value) {
			if (value instanceof UniformizationData) {
				button.setEnabled(false);
			} else {
				button.setEnabled(true);
			}
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			checkButtonStatus(button, value);
			return button;
		}

		@Override
		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
			this.value = (ConformalData)value;
			JButton loadButton = new JButton(ImageHook.getIcon("cog_go.png"));
			checkButtonStatus(loadButton, value);
			loadButton.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					fireEditingStopped();
					loadDataIntoUI(LoadButtonEditor.this.value);
				}
			});
			return loadButton;
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
		if (data instanceof DiscreteMap) {
			DiscreteMap map = (DiscreteMap)data;
			DiscreteEmbedding image = map.getImage();
			DiscreteEmbedding domain = map.getDomain();
			int genus = DataUtility.calculateGenus(domain);
			System.out.println("loading discreet map of a genus " + genus + " surface");
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
			CoHDS hds = DataUtility.toHDS(image, cutInfo);
			for (CoVertex v : hds.getVertices()) {
				EmbeddedVertex dv = domain.getVertices().get(v.getIndex());
				v.T = new double[]{dv.getX(), dv.getY(), dv.getZ(), dv.getW()};
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

	public void addDiscreteTextureEmbedding(CoHDS surface, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		AdapterSet aSet = new AdapterSet(new CoDirectDataAdapter(true));
		addDiscreteEmbedding("Texture Embedding", surface, aSet, TexturePosition4d.class, cutInfo);
	}
	public void addDiscretePositionEmbedding(CoHDS surface, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		AdapterSet aSet = new AdapterSet(new CoDirectDataAdapter(false));
		addDiscreteEmbedding("Position Embedding", surface, aSet, Position4d.class, cutInfo);
	}
	public void addDiscreteMap(String name, CoHDS surface, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		AdapterSet a = new AdapterSet(new CoDirectTextureAdapter(), new CoDirectPositionAdapter());
		DiscreteMap map = DataUtility.toDiscreteMap(name, surface, a, TexturePosition4d.class, Position4d.class, cutInfo);
		addData(map);
	}
	public void addDiscreteEmbedding(String name, CoHDS surface, AdapterSet a, Class<? extends Annotation> type, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		DiscreteEmbedding de = DataUtility.toDiscreteEmbedding(name, surface, a, type, cutInfo);
		addData(de);
	}
	public void addUniformizationData(String name, FundamentalPolygon P) {
		UniformizationData ud = DataUtility.toUniformizationData(name, P);
		addData(ud);
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

	public File showXMLSaveDialog() {
		Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
		int result = fileChooser.showSaveDialog(w);
		if (result != JFileChooser.APPROVE_OPTION) {
			return null;
		}
		File file = fileChooser.getSelectedFile();
		if (!file.getName().toLowerCase().endsWith(".xml")) {
			file = new File(file.getAbsolutePath() + ".xml");
		}
		if (file.exists()) {
			result = JOptionPane.showConfirmDialog(w, "The file " + file.getName() + " exists. \nDo you want to overwrite this file?", "Overwrite?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
			if (result != JOptionPane.OK_OPTION) {
				return null;
			}
		}
		return file;
	}

}
