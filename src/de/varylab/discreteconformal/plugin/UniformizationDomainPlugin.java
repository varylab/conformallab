package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.plugin.TargetGeometry.Hyperbolic;
import static java.awt.BasicStroke.CAP_SQUARE;
import static java.awt.BasicStroke.JOIN_ROUND;
import static java.lang.Double.isNaN;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Shape;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import de.jreality.plugin.job.AbstractJob;
import de.jreality.plugin.job.Job;
import de.jreality.plugin.job.JobQueuePlugin;
import de.jreality.ui.ColorChooseJButton;
import de.jreality.ui.ColorChooseJButton.ColorChangedEvent;
import de.jreality.ui.ColorChooseJButton.ColorChangedListener;
import de.jreality.ui.LayoutFactory;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition3d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeListener;
import de.jtem.halfedgetools.plugin.texturespace.TextureSpacePlugin;
import de.jtem.halfedgetools.util.GeometryUtility;
import de.jtem.java2d.SceneComponent;
import de.jtem.java2dx.Line2DDouble;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.sidecontainer.widget.ShrinkPanel;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.uniformization.VisualizationUtility;

public class UniformizationDomainPlugin extends Plugin implements TextureSpacePlugin, ColorChangedListener, ActionListener, ChangeListener, HalfedgeListener {

	private HalfedgeInterface
		hif = null;
	private JobQueuePlugin
		jobQueue = null;
	
	// active data section
	private CoHDS 
		surface = null; 
	private FundamentalPolygon 
		Pcut = null, 
		Pminimal = null, 
		Pcanonical = null, 
		Popposite = null; 
	private boolean
		uniformizationUpdateRequested = false;
	
	private ShrinkPanel
		options = new ShrinkPanel("Uniformization");
	private JCheckBox
		fundamentalChecker = new JCheckBox("Fundamental Domain", true),
		triangulationChecker = new JCheckBox("Triangulation", true),
		polygonChecker = new JCheckBox("Polygon", true),
		axesChecker = new JCheckBox("Axes", true),
		boundaryChecker = new JCheckBox("Boundary", false),
		faceCirclesChecker = new JCheckBox("Face Circles"),
		vertexCirclesWhiteChecker = new JCheckBox("White Vertex Circles"),
		vertexCirclesBlackChecker = new JCheckBox("Black Vertex Circles");
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
		faceCirclesComponent = new SceneComponent(),
		vertexCirclesWhiteComponent = new SceneComponent(),
		vertexCirclesBlackComponent = new SceneComponent();
	private Shape
		unitCircleShape = new Ellipse2D.Double(-1,-1,2,2);
	
	private SpinnerNumberModel
		coverMaxDistanceModel = new SpinnerNumberModel(0.8, 0.0, 100.0, 0.01),
		coverElementsModel = new SpinnerNumberModel(1, 0, 1000, 1);
	private JSpinner
		coverMaxDistanceSpinner = new JSpinner(coverMaxDistanceModel),
		coverElementsSpinner = new JSpinner(coverElementsModel);
	private JComboBox<DomainPolygon>
		domainCombo = new JComboBox<>(DomainPolygon.values());
	private JComboBox<HyperbolicModel>
		modelCombo = new JComboBox<>(HyperbolicModel.values());
	private JComboBox<InterpolationMethod>
		interpolationCombo = new JComboBox<>(InterpolationMethod.values());
	private JComboBox<TargetGeometry>
		geometryCombo = new JComboBox<>(TargetGeometry.values());
	
	public UniformizationDomainPlugin() {
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
		boundaryComponent.setFilled(false);
		boundaryComponent.setStroke(new BasicStroke(2));
		boundaryComponent.setVisible(false);
		fundamentalDomainComponent.setFilled(true);
		fundamentalDomainComponent.setOutlined(false);
		faceCirclesComponent.setOutlinePaint(Color.DARK_GRAY);
		faceCirclesComponent.setFilled(false);
		scene.addChild(faceCirclesComponent);
		vertexCirclesWhiteComponent.setFilled(false);
		scene.addChild(vertexCirclesWhiteComponent);		
		vertexCirclesBlackComponent.setFilled(false);
		scene.addChild(vertexCirclesBlackComponent);			
		
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
		options.add(vertexCirclesWhiteChecker, rc);
		options.add(vertexCirclesBlackChecker, rc);
		options.add(new JSeparator(JSeparator.HORIZONTAL), rc);
		options.add(new JLabel("Geometry"), lc);
		options.add(geometryCombo, rc);
		options.add(new JLabel("Domain"), lc);
		options.add(domainCombo, rc);
		options.add(new JLabel("Cover Elements"), lc);
		options.add(coverElementsSpinner, rc);
		options.add(new JLabel("Cover Distance"), lc);
		options.add(coverMaxDistanceSpinner, rc);
		options.add(new JLabel("Hyperbolic Model"), lc);
		options.add(modelCombo, rc);
		options.add(new JLabel("Interpolation"), lc);
		options.add(interpolationCombo, rc);
		
		interpolationCombo.addActionListener(this);
		modelCombo.addActionListener(this);
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
		vertexCirclesWhiteChecker.addActionListener(this);
		vertexCirclesBlackChecker.addActionListener(this);
		domainCombo.addActionListener(this);
		coverElementsSpinner.addChangeListener(this);
		coverMaxDistanceSpinner.addChangeListener(this);
		geometryCombo.addActionListener(this);
		
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
		vertexCirclesWhiteComponent.setVisible(vertexCirclesWhiteChecker.isSelected());
		vertexCirclesBlackComponent.setVisible(vertexCirclesBlackChecker.isSelected());
		Color fc = fundamentalColorButton.getColor();
		Color fcAlpha = new Color(fc.getRed(), fc.getGreen(), fc.getBlue(), 51);
		fundamentalDomainComponent.setPaint(fcAlpha);
		vertexCirclesWhiteComponent.setOutlinePaint(fc);
		vertexCirclesBlackComponent.setOutlinePaint(fc);
		scene.fireAppearanceChange();
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		updateStates();
		if (faceCirclesChecker == e.getSource()) {
			updateFaceCircles(hif.get(), hif.getAdapters());
		}
		if (vertexCirclesWhiteChecker == e.getSource()) {
			updateVertexCircles(hif.get(), hif.getAdapters(), true);
		}
		if (vertexCirclesBlackChecker == e.getSource()) {
			updateVertexCircles(hif.get(), hif.getAdapters(), false);
		}
		if (domainCombo == e.getSource()) {
			updateUniformization();
		}
		if (geometryCombo == e.getSource()) {
			if (getSelectedGeometry() == TargetGeometry.Euclidean) {
				setHyperbolicModel(HyperbolicModel.Klein);
			}
			updateGeometry(true);
		}
		if (modelCombo == e.getSource()) {
			updateGeometry(true);
		}
		if (interpolationCombo == e.getSource()) {
			updateGeometry(true);
		}
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (coverElementsSpinner == e.getSource()) {
			updateUniformization();
		}	
		if (coverMaxDistanceSpinner == e.getSource()) {
			updateUniformization();
		}
	}

	private void updateGeometry(boolean updateUniformization) {
		updateAdapters();
		uniformizationUpdateRequested = updateUniformization;
		hif.updateNoUndo();
	}


	public void updateAdapters() {
		CoTexturePositionAdapter texturePositionAdapter = hif.getAdapters().query(CoTexturePositionAdapter.class);
		texturePositionAdapter.setInterpolationMethod(getSelectedInterpolation());
		texturePositionAdapter.setModel(getSelectedHyperbolicModel());
	}
	
	
	public InterpolationMethod getSelectedInterpolation() {
		return (InterpolationMethod)interpolationCombo.getSelectedItem();
	}

	public HyperbolicModel getSelectedHyperbolicModel() {
		return (HyperbolicModel)modelCombo.getSelectedItem();
	}
	
	public void setHyperbolicModel(HyperbolicModel model) {
		try {
			modelCombo.removeActionListener(this);
			modelCombo.setSelectedItem(model);
		} finally {
			modelCombo.addActionListener(this);
		}
		updateAdapters();
	}
	
	public void setGeometry(TargetGeometry geometry) {
		try {
			geometryCombo.removeActionListener(this);
			geometryCombo.setSelectedItem(geometry);
		} finally {
			geometryCombo.addActionListener(this);
		}
		updateAdapters();
	}
	
	public TargetGeometry getSelectedGeometry() {
		return (TargetGeometry)geometryCombo.getSelectedItem();
	}
	
	@Override
	public void colorChanged(ColorChangedEvent cce) {
		updateStates();
	}
	
	public void updateUniformization() {
		if (surface != null) {
			Job j = new AbstractJob() {
				@Override
				public String getJobName() {
					return "Uniformization Visualization";
				}
				@Override
				protected void executeJob() throws Exception {
					createUniformization(surface, Pcut, Pminimal, Pcanonical, Popposite);		
				}
			};
			jobQueue.queueJob(j);
		}
	}
	
	public void createUniformization(
		CoHDS surface, 
		FundamentalPolygon Pcut,
		FundamentalPolygon Pminimal,
		FundamentalPolygon Pcanonical,
		FundamentalPolygon Popposite
	) {
		this.surface = surface;
		this.Pcut = Pcut;
		this.Pminimal = Pminimal;
		this.Pcanonical = Pcanonical;
		this.Popposite = Popposite;
		HyperbolicModel model = getSelectedHyperbolicModel();
		TargetGeometry geometry = getSelectedGeometry();
		int maxElements = coverElementsModel.getNumber().intValue();
		double maxDistance = coverMaxDistanceModel.getNumber().doubleValue();
		DomainPolygon p = (DomainPolygon) domainCombo.getSelectedItem();
		FundamentalPolygon P = null;
		switch (p) {
			case Minimal: P = Pminimal; break;
			case Canonical: P = Pcanonical; break;
			case Opposite: P = Popposite; break;
			default: P = Pcut; break;
		}
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
		boundaryComponent.setVisible(geometry == Hyperbolic);
		
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
		switch (model) {
		case Halfplane:
			Rectangle2D bbox = polyPath.getBounds2D();
			bbox = bbox.createUnion(axesPath.getBounds2D());
			Shape realLineShape = new Line2DDouble(1.2*bbox.getMinX(), 0, 1.2*bbox.getMaxX(), 0); 
			boundaryComponent.setShape(realLineShape);
			break;
		default:
			boundaryComponent.setShape(unitCircleShape);
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
	
	private <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void updateVertexCircles(HDS surface, AdapterSet a, boolean white) {
		if (white) {
			if (!vertexCirclesWhiteChecker.isSelected()) {
				vertexCirclesWhiteComponent.setShape(null);
				return;
			}
		} else {
			if (!vertexCirclesBlackChecker.isSelected()) {
				vertexCirclesBlackComponent.setShape(null);
				return;
			}
		}
		Path2D circles = new Path2D.Double();
		Set<V> lattice = getSublatticeVertices(surface, white);
		for (V v : lattice) {
			V v1 = v.getIncomingEdge().getStartVertex();
			V v2 = v.getIncomingEdge().getNextEdge().getTargetVertex();
			V v3 = v.getIncomingEdge().getOppositeEdge().getPreviousEdge().getStartVertex();
			double[] p1 = a.getD(TexturePosition3d.class, v1); 
			double[] p2 = a.getD(TexturePosition3d.class, v2); 
			double[] p3 = a.getD(TexturePosition3d.class, v3); 
			double[] c = GeometryUtility.circumCircle(p1, p2, p3);
			if (isNaN(c[0]) || isNaN(c[1]) || isNaN(c[2]) || isNaN(c[3])) continue;
			Ellipse2D e = new Ellipse2D.Double(c[0] - c[3], c[1] - c[3], 2*c[3], 2*c[3]);
			circles.append(e, false);
		}
		if (white) {
			vertexCirclesWhiteComponent.setShape(circles);
		} else {
			vertexCirclesBlackComponent.setShape(circles);
		}
	}
	
	
	private <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Set<V> getSublatticeVertices(HDS surface, boolean white) {
		Set<V> result = new HashSet<>();
		Set<V> visited = new HashSet<>();
		if (surface.numVertices() == 0) return visited;
		TreeSet<V> visitQueue = new TreeSet<>();
		visitQueue.add(surface.getVertex(0));
		while (!visitQueue.isEmpty()) {
			V v = visitQueue.pollFirst();
			visited.add(v);
			if (white) result.add(v);
			for (E e : HalfEdgeUtils.incomingEdges(v)) {
				V vv = e.getStartVertex();
				if (!white) result.add(vv);
				for (E ee : HalfEdgeUtils.incomingEdges(vv)) {
					V vvv = ee.getStartVertex();
					if (!visited.contains(vvv)) {
						visitQueue.add(vvv);
					}
				}
			}
		}
		return result;
	}
	
	
	public void reset() {
		surface = null;
		Pcut = null;
		Pminimal = null;
		Pcanonical = null;
		Popposite = null;
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
		if (uniformizationUpdateRequested) {
			updateUniformization();
			uniformizationUpdateRequested = false;
		} else {
			reset();
		}
		updateFaceCircles(layer.get(), layer.getEffectiveAdapters());
		updateVertexCircles(hif.get(), hif.getAdapters(), true);
		updateVertexCircles(hif.get(), hif.getAdapters(), false);
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
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
	}
	
	@Override
	public void restoreStates(Controller c) throws Exception {
		super.restoreStates(c);
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		hif.addHalfedgeListener(this);
		jobQueue = c.getPlugin(JobQueuePlugin.class);
	}
	
	public int getMaxCoverElements() {
		return coverElementsModel.getNumber().intValue();
	}
	
	public double getMaxCoverDistance() {
		return coverMaxDistanceModel.getNumber().doubleValue();
	}

}
