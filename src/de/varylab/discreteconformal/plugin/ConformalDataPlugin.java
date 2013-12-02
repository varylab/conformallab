package de.varylab.discreteconformal.plugin;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;

import de.jreality.plugin.basic.View;
import de.jtem.halfedgetools.plugin.image.ImageHook;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.conformallab.data.DataFactory;
import de.varylab.conformallab.data.types.ConformalData;
import de.varylab.conformallab.data.types.ConformalDataList;

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
		exportButton = new JButton("Export...", ImageHook.getIcon("disk.png"));
	private JFileChooser
		exportFileChooser = new JFileChooser();
	
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
		dataScroller.setPreferredSize(new Dimension(10, 200));
		dataTable.setDefaultRenderer(ConformalData.class, new DataCellRenderer());
		dataTable.getTableHeader().setPreferredSize(new Dimension(0, 0));
		c.weighty = 0.0;
		shrinkPanel.add(exportButton, c);
		exportButton.addActionListener(this);
		
		exportFileChooser.setAcceptAllFileFilterUsed(true);
		exportFileChooser.setFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "Conformal Data (*.xml)";
			}
			@Override
			public boolean accept(File f) {
				return f.isDirectory() || f.getName().toLowerCase().endsWith(".xml");
			}
		});
		exportFileChooser.setDialogTitle("Export Conformal Data");
		exportFileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if (exportButton == e.getSource()) {
			Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
			int result = exportFileChooser.showSaveDialog(w);
			if (result != JFileChooser.APPROVE_OPTION) {
				return;
			}
			File file = exportFileChooser.getSelectedFile();
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
	}
	
	private class DataModel extends DefaultTableModel {

		private static final long serialVersionUID = 1L;

		@Override
		public int getRowCount() {
			return data.size();
		}
		@Override
		public int getColumnCount() {
			return 1;
		}
		@Override
		public Class<?> getColumnClass(int columnIndex) {
			return ConformalData.class;
		}
		@Override
		public Object getValueAt(int row, int column) {
			return data.get(row);
		}
		
	}
	
	private class DataCellRenderer extends DefaultTableCellRenderer {
		
		private static final long serialVersionUID = 1L;

		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			Component c = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
			if (value instanceof ConformalData && c instanceof JLabel) {
				JLabel label = (JLabel)c;
				label.setText(value.getClass().getSimpleName());
			}
			return c;
		}
		
	}
	
	private void updateUI() {
		dataModel.fireTableDataChanged();
	}
	
	
	public void addData(ConformalData data) {
		this.data.add(data);
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
		c.storeProperty(getClass(), "exportDir", exportFileChooser.getCurrentDirectory().getAbsolutePath());
	}
	
	@Override
	public void restoreStates(Controller c) throws Exception {
		super.restoreStates(c);
		exportFileChooser.setCurrentDirectory(new File(c.getProperty(getClass(), "exportDir", ".")));
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

}
