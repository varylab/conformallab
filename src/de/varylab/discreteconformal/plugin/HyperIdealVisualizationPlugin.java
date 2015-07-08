package de.varylab.discreteconformal.plugin;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;

import javax.swing.JCheckBox;

import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Radius;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition2d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeListener;
import de.jtem.halfedgetools.plugin.texturespace.TextureSpacePlugin;
import de.jtem.java2d.SceneComponent;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.sidecontainer.widget.ShrinkPanel;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.HyperIdealRadiusAdapter;

public class HyperIdealVisualizationPlugin extends Plugin implements TextureSpacePlugin, HalfedgeListener, ActionListener {

	private HalfedgeInterface
		hif = null;
	private SceneComponent
		root = new SceneComponent();
	private ShrinkPanel
		ui = new ShrinkPanel("Hyperideal Uniformization");
	private JCheckBox
		vertexChecker = new JCheckBox("Vertex Circles", true),
		faceChecker = new JCheckBox("Face Circles", true);
	
	public HyperIdealVisualizationPlugin() {
		ui.setLayout(new GridLayout(2, 1));
		ui.add(vertexChecker);
		ui.add(faceChecker);
		vertexChecker.addActionListener(this);
		faceChecker.addActionListener(this);
	}

	@Override
	public SceneComponent getSceneComponent() {
		return root;
	}

	@Override
	public ShrinkPanel getOptionPanel() {
		return ui;
	}

	protected void updateCircles(CoHDS hds) {
		root.removeAllChildren();
		AdapterSet a = hif.getActiveAdapters();
		boolean hasRadii = a.isAvailable(HyperIdealRadiusAdapter.class);
		if (!hasRadii) {
			root.fireAppearanceChange();
			return;
		}
		if (vertexChecker.isSelected()) {
			Path2D vertexCircles = new Path2D.Float();
			for (CoVertex v : hds.getVertices()) {
				Double r = a.get(Radius.class, v, Double.class);
				double[] center = a.getD(TexturePosition2d.class, v);
				Ellipse2D circle = new Ellipse2D.Double(center[0] - r, center[1] - r, 2*r, 2*r);
				vertexCircles.append(circle, false);
			}
			SceneComponent vertexComponent = new SceneComponent();
			vertexComponent.setShape(vertexCircles);
			root.addChild(vertexComponent);
		}
		if (faceChecker.isSelected()) {
			Path2D faceCircles = new Path2D.Float();
			
			SceneComponent faceComponent = new SceneComponent();
			faceComponent.setShape(faceCircles);
			root.addChild(faceComponent);
		}
		root.fireAppearanceChange();
	}
	
	@Override
	public boolean getRenderOnTop() {
		return true;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		updateCircles(hif.get(new CoHDS()));
	}
	
	@Override
	public void dataChanged(HalfedgeLayer layer) {
		updateCircles(layer.get(new CoHDS()));
	}
	@Override
	public void adaptersChanged(HalfedgeLayer layer) {
	}
	@Override
	public void activeLayerChanged(HalfedgeLayer old, HalfedgeLayer active) {
		updateCircles(active.get(new CoHDS()));
	}
	@Override
	public void layerCreated(HalfedgeLayer layer) {
	}
	@Override
	public void layerRemoved(HalfedgeLayer layer) {
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		hif.addHalfedgeListener(this);
	}

}
