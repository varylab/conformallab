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

import de.jreality.math.Rn;
import de.jreality.plugin.JRViewerUtility;
import de.jreality.plugin.basic.Scene;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.basic.ViewMenuBar;
import de.jtem.halfedge.algorithm.Coord3DAdapter;
import de.jtem.halfedge.algorithm.catmullclark.CatmullClarkSubdivision;
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
	private HalfedgeConnectorPlugin<CoVertex, CoEdge, CoFace, CoHDS>
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
				e3.setTargetVertex(v2a);
				break;
			case 1:
				e1.setTargetVertex(v1a);
				e2.setTargetVertex(v3a);
				e3.setTargetVertex(v2a);
				break;
			case 2:
				e1.setTargetVertex(v2a);
				e2.setTargetVertex(v3a);
				e3.setTargetVertex(v0a);
				break;
			case 3:
				e1.setTargetVertex(v3a);
				e2.setTargetVertex(v1a);
				e3.setTargetVertex(v0a);		
				break;
			}
		}
		
		hds.getEdge(0).linkOppositeEdge(hds.getEdge(6));
		hds.getEdge(1).linkOppositeEdge(hds.getEdge(23));
		hds.getEdge(2).linkOppositeEdge(hds.getEdge(3));
		hds.getEdge(4).linkOppositeEdge(hds.getEdge(10));
		hds.getEdge(5).linkOppositeEdge(hds.getEdge(19));
		hds.getEdge(8).linkOppositeEdge(hds.getEdge(9));
		
		
		hds.getEdge(12).linkOppositeEdge(hds.getEdge(18));
		hds.getEdge(13).linkOppositeEdge(hds.getEdge(11));
		hds.getEdge(14).linkOppositeEdge(hds.getEdge(15));
		hds.getEdge(16).linkOppositeEdge(hds.getEdge(22));
		hds.getEdge(17).linkOppositeEdge(hds.getEdge(7));
		hds.getEdge(20).linkOppositeEdge(hds.getEdge(21));
		
		CatmullClarkSubdivision<CoVertex, CoEdge, CoFace>
			subD = new CatmullClarkSubdivision<CoVertex, CoEdge, CoFace>();
		CoHDS hds2 = hds;
		for (int i = 0; i < genus; i++) {
			hds2 = new CoHDS();
			subD.subdivide(hds, hds2, new MyPositionAdapter());
			hds = hds2;
		}
		
		for (CoVertex v : hds2.getVertices()) {
			double[] pos = v.getPosition().get();
			double l = Rn.euclideanNorm(pos);
			Rn.times(pos, 1 / l, pos);
		}
		
		if (HalfEdgeUtils.isValidSurface(hds2, true)) {
			halfedgeConnectorPlugin.updateHalfedgeContent(hds2, true, new PositionAdapter());
			JRViewerUtility.encompassEuclidean(scene);
		}
	}
	
	private class MyPositionAdapter implements Coord3DAdapter<CoVertex> {

		@Override
		public double[] getCoord(CoVertex v) {
			return v.getPosition().get();
		}

		@Override
		public void setCoord(CoVertex v, double[] c) {
			v.getPosition().set(c);
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
			String genusString = JOptionPane.showInputDialog(w, "Subdivision Depth?", 1);
			int genus = Integer.parseInt(genusString);
			generateCurve(genus);
		}
		
	}

	@SuppressWarnings("unchecked")
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
