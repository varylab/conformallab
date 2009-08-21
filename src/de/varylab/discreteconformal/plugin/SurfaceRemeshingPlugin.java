package de.varylab.discreteconformal.plugin;

import static javax.swing.JOptionPane.WARNING_MESSAGE;

import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import de.jreality.plugin.basic.View;
import de.jreality.plugin.experimental.ManagedContent;
import de.jtem.halfedge.jreality.adapter.Adapter;
import de.jtem.halfedge.jreality.adapter.Adapter.AdapterType;
import de.jtem.halfedge.jreality.standard.MyCoordinateAdapter;
import de.jtem.halfedge.jreality.standard.MyTextCoordAdapter;
import de.jtem.halfedge.jreality.standard.node.MyHDS;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.varylab.jrworkspace.plugin.Controller;
import de.varylab.jrworkspace.plugin.PluginInfo;
import de.varylab.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.varylab.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;

public class SurfaceRemeshingPlugin extends ShrinkPanelPlugin implements ActionListener {

	public static abstract class MeshPattern {
		  

		
	}
	
	
	// plug-in connection
	private ManagedContent
		managedContent = null;
	private HalfedgeConnectorPlugin
		hcp = null;
	
	// ui components
	private GridBagConstraints
		gbc1 = new GridBagConstraints(),
		gbc2 = new GridBagConstraints();
	private JButton
		meshingButton = new JButton("Remesh Surface");
	
	public SurfaceRemeshingPlugin() {
		gbc1.insets = new Insets(2, 2, 2, 2);
		gbc1.gridwidth = GridBagConstraints.RELATIVE;
		gbc1.weightx = 0.0;
		gbc2.insets = new Insets(2, 2, 2, 2);
		gbc2.gridwidth = GridBagConstraints.REMAINDER;
		gbc2.weightx = 1.0;
		shrinkPanel.add(meshingButton);

		meshingButton.addActionListener(this);
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == meshingButton) {
			remeshSurface();
		}
	}
	
	
	private void remeshSurface() {
		Adapter a1 = new MyCoordinateAdapter(AdapterType.VERTEX_ADAPTER);
		Adapter a2 = new MyTextCoordAdapter(AdapterType.VERTEX_ADAPTER);
		MyHDS hds = hcp.getHalfedgeContent(a1, a2);
		if (hds == null) {
			return;
		}
		Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
		if (hds.getVertex(0).textCoord == null) {
			JOptionPane.showMessageDialog(w, "Surface has no texture coordinates.", "Error", WARNING_MESSAGE);
			return;
		}
		
	}
	
	
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hcp = c.getPlugin(HalfedgeConnectorPlugin.class);
		managedContent = c.getPlugin(ManagedContent.class);
	}
	
	
	@Override
	public void uninstall(Controller c) throws Exception {
		super.uninstall(c);
		managedContent.removeAll(getClass());
	}
	
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = new PluginInfo("Surface Remeshing", "Stefan Sechelmann");
		return info;
	}

}
