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
import de.jtem.halfedgetools.jreality.adapter.Adapter;
import de.jtem.halfedgetools.jreality.adapter.Adapter.AdapterType;
import de.jtem.halfedgetools.jreality.adapter.standard.StandardCoordinateAdapter;
import de.jtem.halfedgetools.jreality.adapter.standard.StandardTextCoordAdapter;
import de.jtem.halfedgetools.jreality.node.standard.StandardHDS;
import de.jtem.halfedgetools.plugin.HalfedgeConnectorPlugin;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class SurfaceRemeshingPlugin extends ShrinkPanelPlugin implements ActionListener {

	public static abstract class MeshPattern {
		  

		
	}
	
	// plug-in connection
	private ManagedContent
		managedContent = null;
	private HalfedgeConnectorPlugin<CoVertex, CoEdge, CoFace, CoHDS>
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
		Adapter a1 = new StandardCoordinateAdapter(AdapterType.VERTEX_ADAPTER);
		Adapter a2 = new StandardTextCoordAdapter(AdapterType.VERTEX_ADAPTER);
		StandardHDS hds = hcp.getHalfedgeContent(a1, a2);
		if (hds == null) {
			return;
		}
		Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
		if (hds.getVertex(0).textCoord == null) {
			JOptionPane.showMessageDialog(w, "Surface has no texture coordinates.", "Error", WARNING_MESSAGE);
			return;
		}
		
	}
	
	
	
	@SuppressWarnings("unchecked")
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
