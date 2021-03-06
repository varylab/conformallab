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
import java.io.IOException;
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
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.selection.Selection;
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
import de.varylab.conformallab.data.types.HalfedgeEmbedding;
import de.varylab.conformallab.data.types.HalfedgeMap;
import de.varylab.conformallab.data.types.HalfedgeSelection;
import de.varylab.conformallab.data.types.HalfedgeVertex;
import de.varylab.conformallab.data.types.HyperEllipticAlgebraicCurve;
import de.varylab.conformallab.data.types.SchottkyData;
import de.varylab.conformallab.data.types.UniformizationData;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoDirectPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoDirectTextureAdapter;
import de.varylab.discreteconformal.heds.adapter.MappedEdgeLengthAdapter;
import de.varylab.discreteconformal.plugin.image.ImageHook;
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
	private UniformizationDomainPlugin
		domainPlugin = null;
	
	private DataModel
		dataModel = new DataModel();
	private JTable
		dataTable = new JTable(dataModel);
	private JScrollPane
		dataScroller = new JScrollPane(dataTable);
	private JButton
		exportButton = new JButton("Export...", ImageHook.getIcon("disk.png")),
		importButton = new JButton("Import...", ImageHook.getIcon("folder.png")),
		addActiveButton = new JButton("Add Active", ImageHook.getIcon("add.png")),
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
		c.gridwidth = GridBagConstraints.RELATIVE;
		shrinkPanel.add(addActiveButton, c);
		addActiveButton.addActionListener(this);
		c.gridwidth = GridBagConstraints.REMAINDER;
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
			fileChooser.setDialogTitle("Import Conformal Data");
			int result = fileChooser.showOpenDialog(w);
			if (result != JFileChooser.APPROVE_OPTION) {
				return;
			}
			File file = fileChooser.getSelectedFile();
			FileInputStream fin = null;
			try {
				fin = new FileInputStream(file);
				Object data = DataIO.read(fin);
				clearData();
				if (data instanceof ConformalData) {
					addData((ConformalData)data);
				} else if (data instanceof ConformalDataList) {
					ConformalDataList dataList = (ConformalDataList)data;
					for (ConformalData d : dataList.getData()) {
						addData(d);
					}
				}
				fin.close();
			} catch (Exception e1) {
				log.warning("Could not import conformal data" + e1);
			} finally {
				try {
					fin.close();
				} catch (IOException e1) {
					log.warning("could not close input stream: " + e1);
				}
			}
		}
		if (clearButton == e.getSource()) {
			clearData();
		}
		if (addActiveButton == e.getSource()) {
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = discreteConformalPlugin.getCurrentCutInfo();
			addHalfedgeMap("Map", hif.get(new CoHDS()), hif.getSelection(), cutInfo);
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
		
		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			return button;
		}

		@Override
		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
			this.value = (ConformalData)value;
			JButton loadButton = new JButton(ImageHook.getIcon("cog_go.png"));
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
			log.info("loading embedding of genus " + genus);
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
			CoHDS hds = new CoHDS();
			DataUtility.toHalfedge(de, hif.getAdapters(), Position.class, hds, cutInfo);
			for (CoVertex v : hds.getVertices()) {
				System.arraycopy(v.P, 0, v.T, 0, 4);
			}
			TargetGeometry target = TargetGeometry.calculateTargetGeometry(genus, 0);
			discreteConformalPlugin.createUniformization(hds, target, cutInfo);
			if (de.getSelection() != null) {
				Selection s = DataUtility.toSelection(de.getSelection(), hds);
				hif.setSelection(s);
			}
		}
		if (data instanceof DiscreteMap) {
			DiscreteMap map = (DiscreteMap)data;
			DiscreteEmbedding image = map.getImage();
			DiscreteEmbedding domain = map.getDomain();
			int genus = DataUtility.calculateGenus(domain);
			log.info("loading discrete map of a genus " + genus + " surface");
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
			CoHDS hds = new CoHDS();
			DataUtility.toHalfedge(image, hif.getAdapters(), Position.class, hds, cutInfo);
			for (CoVertex v : hds.getVertices()) {
				EmbeddedVertex dv = domain.getVertices().get(v.getIndex());
				v.T = new double[]{dv.getX(), dv.getY(), dv.getZ(), dv.getW()};
			}
			TargetGeometry target = TargetGeometry.calculateTargetGeometry(genus, 0);
			discreteConformalPlugin.createUniformization(hds, target, cutInfo);
		}
		if (data instanceof HalfedgeEmbedding) {
			HalfedgeEmbedding he = (HalfedgeEmbedding)data;
			int genus = DataUtility.calculateGenus(he);
			log.info("loading half-edge embedding of genus " + genus);
			CoHDS hds = new CoHDS();
			CuttingInfo<CoVertex, CoEdge, CoFace> hdsCutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
			DataUtility.toHalfedge(he, hif.getAdapters(), TexturePosition.class, hds, hdsCutInfo);
			for (CoVertex v : hds.getVertices()) {
				System.arraycopy(v.T, 0, v.P, 0, 4);
			}
			TargetGeometry target = TargetGeometry.calculateTargetGeometry(genus, 0);
			discreteConformalPlugin.createUniformization(hds, target, hdsCutInfo);
			if (he.getSelection() != null) {
				Selection s = DataUtility.toSelection(he.getSelection(), hds);
				hif.setSelection(s);
			}
		}
		if (data instanceof HalfedgeMap) {
			HalfedgeMap map = (HalfedgeMap)data;
			HalfedgeEmbedding image = map.getImage();
			HalfedgeEmbedding domain = map.getDomain();
			int genus = DataUtility.calculateGenus(domain);
			log.info("loading halfedge map of a genus " + genus + " surface");
			CoHDS domainHds = new CoHDS();
			CuttingInfo<CoVertex, CoEdge, CoFace> domainInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
			DataUtility.toHalfedge(domain, hif.getAdapters(), TexturePosition.class, domainHds, domainInfo);
			for (CoVertex v : domainHds.getVertices()) {
				HalfedgeVertex iv = image.getVertices().get(v.getIndex());
				v.P = new double[]{iv.getX(), iv.getY(), iv.getZ(), iv.getW()};
			}
			TargetGeometry target = TargetGeometry.calculateTargetGeometry(genus, 0);
			discreteConformalPlugin.createUniformization(domainHds, target, domainInfo);
			if (domain.getSelection() != null) {
				Selection s = DataUtility.toSelection(domain.getSelection(), domainHds);
				hif.setSelection(s);
			}
		}		
		if (data instanceof DiscreteMetric) {
			DiscreteMetric dm = (DiscreteMetric)data;
			MappedEdgeLengthAdapter lMap = new MappedEdgeLengthAdapter(1000.0);
			CoHDS hds = DataUtility.toHalfedgeAndLengths(dm, lMap);
			log.info("loading discrete metric of genus " + HalfEdgeUtils.getGenus(hds));
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
		if (data instanceof UniformizationData) {
			UniformizationData ud = (UniformizationData)data;
			FundamentalPolygon p = DataUtility.toFundamentalPolygon(ud);
			CoHDS surface = hif.get(new CoHDS());
			domainPlugin.createUniformization(surface, p, p, p, p, p);
		}
	}

	public void addDiscreteMetric(String name, CoHDS surface, AdapterSet a) {
		DiscreteMetric dm = DataUtility.toDiscreteMetric(name, surface, a);
		addData(dm);
	}

	@Deprecated
	public void addDiscreteMap(String name, CoHDS surface, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		AdapterSet a = new AdapterSet(new CoDirectTextureAdapter(), new CoDirectPositionAdapter());
		DiscreteMap map = DataUtility.toDiscreteMap(name, surface, a, TexturePosition4d.class, Position4d.class, cutInfo);
		addData(map);
	}
	@Deprecated
	public void addDiscreteEmbedding(String name, CoHDS surface, Selection selection, AdapterSet a, Class<? extends Annotation> type, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		DiscreteEmbedding de = DataUtility.toDiscreteEmbedding(name, surface, a, type, cutInfo);
		if (selection != null) {
			de.setSelection(DataUtility.toEmbeddingSelection(selection));
		}
		addData(de);
	}
	public void addHalfedgeMap(String name, CoHDS surface, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		addHalfedgeMap(name, surface, new Selection(), cutInfo);
	}
	public void addHalfedgeMap(String name, CoHDS surface, Selection selection, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		AdapterSet a = new AdapterSet(new CoDirectTextureAdapter(), new CoDirectPositionAdapter());
		HalfedgeMap map = DataUtility.toHalfedgeMap(name, surface, a, TexturePosition4d.class, Position4d.class, cutInfo);
		if (selection != null) {
			HalfedgeSelection hs = DataUtility.toHalfedgeSelection(selection);
			map.getDomain().setSelection(hs);
			map.getImage().setSelection(hs);
		}
		addData(map);
	}	
	public void addHalfedgeEmbedding(String name, CoHDS surface, Selection selection, AdapterSet a, Class<? extends Annotation> type, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		HalfedgeEmbedding he = DataUtility.toHalfedgeEmbedding(name, surface, a, type, cutInfo);
		if (selection != null) {
			he.setSelection(DataUtility.toHalfedgeSelection(selection));
		}
		addData(he);
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
		domainPlugin = c.getPlugin(UniformizationDomainPlugin.class);
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

	public File showXMLSaveDialog() {
		Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
		fileChooser.setDialogTitle("Export Conformal Data");
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
