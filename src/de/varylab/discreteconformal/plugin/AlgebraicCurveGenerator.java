package de.varylab.discreteconformal.plugin;

import static java.lang.Math.PI;
import static java.lang.Math.sin;

import java.awt.Window;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import de.jreality.plugin.JRViewerUtility;
import de.jreality.plugin.basic.Scene;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.basic.ViewMenuBar;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.jrworkspace.plugin.Controller;
import de.varylab.jrworkspace.plugin.Plugin;
import de.varylab.jrworkspace.plugin.PluginInfo;

public class AlgebraicCurveGenerator extends Plugin {

	private ViewMenuBar
		viewMenuBar = null;
	private Scene
		scene = null;
	private View
		view = null;
	private HalfedgeConnectorPlugin
		halfedgeConnectorPlugin = null;
	private GenerateAction
		generateAction = new GenerateAction();

	
	private void generateCurve(int genus) {
		CoHDS hds = new CoHDS();
		// first copy
		CoVertex v0a = hds.addNewVertex();
		CoVertex v1a = hds.addNewVertex();
		CoVertex v2a = hds.addNewVertex();
		CoVertex v3a = hds.addNewVertex();
		double a = 2 * sin(PI / 4);
		v0a.getPosition().set(a, 1, 0);
		v1a.getPosition().set(-a, 1, 0);
		v2a.getPosition().set(0, -1, -a);
		v3a.getPosition().set(0, -1, a);
		
		List<CoFace> faces = new ArrayList<CoFace>();
		for (int i = 0; i < 8; i++) {
			CoFace f = hds.addNewFace();
			faces.add(f);
			CoEdge e1 = hds.addNewEdge();
			CoEdge e2 = hds.addNewEdge();
			CoEdge e3 = hds.addNewEdge();
			e1.linkNextEdge(e2);
			e2.linkNextEdge(e3);
			e3.linkNextEdge(e1);
			e1.setLeftFace(f);
			e2.setLeftFace(f);
			e3.setLeftFace(f);
			switch (i % 4) {
			case 0:
				e1.setTargetVertex(v0a);
				e2.setTargetVertex(v1a);
				e3.setTargetVertex(v3a);
				break;
			case 1:
				e1.setTargetVertex(v0a);
				e2.setTargetVertex(v1a);
				e3.setTargetVertex(v2a);
				break;
			case 2:
				e1.setTargetVertex(v2a);
				e2.setTargetVertex(v3a);
				e3.setTargetVertex(v0a);
				break;
			case 3:
				e1.setTargetVertex(v2a);
				e2.setTargetVertex(v3a);
				e3.setTargetVertex(v1a);		
				break;
			}
		}

		
		if (HalfEdgeUtils.isValidSurface(hds, true)) {
			halfedgeConnectorPlugin.updateHalfedgeContent(hds, true, new PositionAdapter());
			JRViewerUtility.encompassEuclidean(scene);
		}
	}
	
	
	private class GenerateAction extends AbstractAction {
		
		private static final long 
			serialVersionUID = 1L;

		public GenerateAction() {
			putValue(NAME, "Generate Algebraic Curve");
		}
		
		@Override
		public void actionPerformed(ActionEvent e) {
			Window w = SwingUtilities.getWindowAncestor(view.getCenterComponent());
			String genusString = JOptionPane.showInputDialog(w, "Which genus?", 1);
			int genus = Integer.parseInt(genusString);
			generateCurve(genus);
		}
		
	}

	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		view = c.getPlugin(View.class);
		scene = c.getPlugin(Scene.class);
		viewMenuBar = c.getPlugin(ViewMenuBar.class);
		viewMenuBar.addMenuItem(getClass(), 0, generateAction, "Halfedge");
		halfedgeConnectorPlugin = c.getPlugin(HalfedgeConnectorPlugin.class);
	}
	
	
	@Override
	public void uninstall(Controller c) throws Exception {
		super.uninstall(c);
		viewMenuBar.removeAll(getClass());
	}
	
	@Override
	public PluginInfo getPluginInfo() {
		return new PluginInfo("Algebraic Curve Generator", "Stefan Sechelmann");
	}

}
