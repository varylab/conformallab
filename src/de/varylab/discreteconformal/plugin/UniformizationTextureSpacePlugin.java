package de.varylab.discreteconformal.plugin;

import static java.awt.BasicStroke.CAP_SQUARE;
import static java.awt.BasicStroke.JOIN_ROUND;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;

import javax.swing.JCheckBox;

import de.jreality.ui.ColorChooseJButton;
import de.jreality.ui.ColorChooseJButton.ColorChangedEvent;
import de.jreality.ui.ColorChooseJButton.ColorChangedListener;
import de.jreality.ui.LayoutFactory;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition3d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeListener;
import de.jtem.halfedgetools.plugin.texturespace.TextureSpacePlugin;
import de.jtem.halfedgetools.util.GeometryUtility;
import de.jtem.java2d.SceneComponent;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.sidecontainer.widget.ShrinkPanel;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.uniformization.VisualizationUtility;

public class UniformizationTextureSpacePlugin extends Plugin implements TextureSpacePlugin, ColorChangedListener, ActionListener, HalfedgeListener {

	private HalfedgeInterface
		hif = null;
	
	private ShrinkPanel
		options = new ShrinkPanel("Uniformization");
	private JCheckBox
		fundamentalChecker = new JCheckBox("Fundamental Domain", true),
		triangulationChecker = new JCheckBox("Triangulation", true),
		polygonChecker = new JCheckBox("Polygon", true),
		axesChecker = new JCheckBox("Axes", true),
		boundaryChecker = new JCheckBox("Boundary", false),
		faceCirclesChecker = new JCheckBox("Face Circles");
	private ColorChooseJButton
		fundamentalColorButton = new ColorChooseJButton(new Color(0, 102, 204), true),
		triangulationColorButton = new ColorChooseJButton(new Color(102, 102, 102), true),
		polygonColorButton = new ColorChooseJButton(new Color(204, 102, 0), true),
		boundaryColorButton = new ColorChooseJButton(new Color(204, 102, 0), true),
		axesColorButton = new ColorChooseJButton(new Color(0, 153, 204), true);
	private SceneComponent
		scene = new SceneComponent(),
		boundaryComponent = new SceneComponent(),
		fundamentalDomainComponent = new SceneComponent(),
		axesComponent = new SceneComponent(),
		polygonComponent = new SceneComponent(),
		boundaryEdgesComponent = new SceneComponent(),
		triangulationComponent = new SceneComponent(),
		faceCirclesComponent = new SceneComponent();
	
	public UniformizationTextureSpacePlugin() {
		scene.addChild(fundamentalDomainComponent);
		scene.addChild(triangulationComponent);
		triangulationComponent.setFilled(false);
		scene.addChild(boundaryEdgesComponent);
		boundaryEdgesComponent.setFilled(false);
		boundaryEdgesComponent.setStroke(new BasicStroke(2));
		scene.addChild(axesComponent);
		axesComponent.setFilled(false);
		axesComponent.setStroke(new BasicStroke(2, CAP_SQUARE, JOIN_ROUND, 1, new float[] {5, 7}, 1));
		scene.addChild(polygonComponent);
		polygonComponent.setFilled(false);
		polygonComponent.setStroke(new BasicStroke(2));
		scene.addChild(boundaryComponent);
		boundaryComponent.setOutlinePaint(Color.BLACK);
		boundaryComponent.setShape(new Ellipse2D.Double(-1,-1,2,2));
		boundaryComponent.setFilled(false);
		boundaryComponent.setStroke(new BasicStroke(2));
		boundaryComponent.setVisible(false);
		fundamentalDomainComponent.setFilled(true);
		fundamentalDomainComponent.setOutlined(false);
		faceCirclesComponent.setOutlinePaint(Color.DARK_GRAY);
		faceCirclesComponent.setFilled(false);
		scene.addChild(faceCirclesComponent);
		
		GridBagConstraints lc = LayoutFactory.createLeftConstraint();
		GridBagConstraints rc = LayoutFactory.createRightConstraint();
		options.setLayout(new GridBagLayout());
		options.add(triangulationChecker, lc);
		options.add(triangulationColorButton, rc);
		options.add(fundamentalChecker, lc);
		options.add(fundamentalColorButton, rc);
		options.add(polygonChecker, lc);
		options.add(polygonColorButton, rc);
		options.add(axesChecker, lc);
		options.add(axesColorButton, rc);
		options.add(boundaryChecker, lc);
		options.add(boundaryColorButton, rc);	
		options.add(faceCirclesChecker, rc);
		
		triangulationChecker.addActionListener(this);
		triangulationColorButton.addColorChangedListener(this);
		polygonChecker.addActionListener(this);
		polygonColorButton.addColorChangedListener(this);
		axesChecker.addActionListener(this);
		axesColorButton.addColorChangedListener(this);
		boundaryChecker.addActionListener(this);
		boundaryColorButton.addColorChangedListener(this);
		fundamentalChecker.addActionListener(this);
		fundamentalColorButton.addColorChangedListener(this);
		faceCirclesChecker.addActionListener(this);
		
		updateStates();
	}

	private void updateStates() {
		triangulationComponent.setVisible(triangulationChecker.isSelected());
		triangulationComponent.setOutlinePaint(triangulationColorButton.getColor());
		polygonComponent.setVisible(polygonChecker.isSelected());
		polygonComponent.setOutlinePaint(polygonColorButton.getColor());
		axesComponent.setVisible(axesChecker.isSelected());
		axesComponent.setOutlinePaint(axesColorButton.getColor());
		boundaryEdgesComponent.setVisible(boundaryChecker.isSelected());
		boundaryEdgesComponent.setOutlinePaint(boundaryColorButton.getColor());
		fundamentalDomainComponent.setVisible(fundamentalChecker.isSelected());
		faceCirclesComponent.setVisible(faceCirclesChecker.isSelected());
		Color fc = fundamentalColorButton.getColor();
		Color fcAlpha = new Color(fc.getRed(), fc.getGreen(), fc.getBlue(), 51);
		fundamentalDomainComponent.setPaint(fcAlpha);
		scene.fireAppearanceChange();
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		updateStates();
		if (faceCirclesChecker == e.getSource()) {
			updateFaceCircles(hif.get(), hif.getAdapters());
		}
	}
	
	@Override
	public void colorChanged(ColorChangedEvent cce) {
		updateStates();
	}
	
	public void createUniformization(CoHDS surface, FundamentalPolygon P, int maxElements, double maxDistance, HyperbolicModel model, TargetGeometry geometry) {
		Path2D axesPath = new Path2D.Float();
		Path2D polyPath = new Path2D.Float();
		Path2D triangulationPath = new Path2D.Float();
		Path2D boundaryPath = new Path2D.Float();
		Path2D fundamentalDomainPath = new Path2D.Float();
		VisualizationUtility.createUniversalCover(
			P, 
			model, 
			maxElements, maxDistance, 
			true, true, 
			null, null,
			axesPath, polyPath, fundamentalDomainPath
		);
		VisualizationUtility.createTriangulation(
			surface,
			P, 
			model, 
			maxElements, maxDistance,
			triangulationPath,
			boundaryPath
		);
		triangulationComponent.setShape(triangulationPath);
		boundaryEdgesComponent.setShape(boundaryPath);
		axesComponent.setShape(axesPath);
		polygonComponent.setShape(polyPath);
		fundamentalDomainComponent.setShape(fundamentalDomainPath);
		boundaryComponent.setVisible(geometry == TargetGeometry.Hyperbolic);
		
		if (geometry == TargetGeometry.Automatic) {
			geometry = TargetGeometry.calculateTargetGeometry(surface);
		}
		switch (geometry) {
		case Euclidean:
			axesComponent.setShape(null);
			axesComponent.setVisible(false);
			break;
		default:
			break;
		}
		
		scene.fireAppearanceChange();
	}
	
	
	private <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void updateFaceCircles(HDS surface, AdapterSet a) {
		if (!faceCirclesChecker.isSelected()) {
			faceCirclesComponent.setShape(null);
			return;
		}
		Path2D circles = new Path2D.Double();
		for (F f : surface.getFaces()) {
			double[] p1 = a.getD(TexturePosition3d.class, f.getBoundaryEdge().getTargetVertex()); 
			double[] p2 = a.getD(TexturePosition3d.class, f.getBoundaryEdge().getNextEdge().getTargetVertex()); 
			double[] p3 = a.getD(TexturePosition3d.class, f.getBoundaryEdge().getPreviousEdge().getTargetVertex()); 
			double[] c = GeometryUtility.circumCircle(p1, p2, p3);
			Ellipse2D e = new Ellipse2D.Double(c[0] - c[3], c[1] - c[3], 2*c[3], 2*c[3]);
			circles.append(e, false);
		}
		faceCirclesComponent.setShape(circles);
	}
	
	public void reset() {
		boundaryComponent.setVisible(false);
		triangulationComponent.setShape(null);
		axesComponent.setShape(null);
		polygonComponent.setShape(null);
		fundamentalDomainComponent.setShape(null);
		faceCirclesComponent.setShape(null);
		scene.fireAppearanceChange();
	}
	
	@Override
	public SceneComponent getSceneComponent() {
		return scene;
	}

	@Override
	public ShrinkPanel getOptionPanel() {
		return options;
	}

	@Override
	public boolean getRenderOnTop() {
		return true;
	}

	@Override
	public void dataChanged(HalfedgeLayer layer) {
		reset();
		updateFaceCircles(layer.get(), layer.getEffectiveAdapters());
	}
	@Override
	public void adaptersChanged(HalfedgeLayer layer) {
	}
	@Override
	public void activeLayerChanged(HalfedgeLayer old, HalfedgeLayer active) {
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
