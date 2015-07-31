package de.varylab.discreteconformal.plugin;

import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.LINE_SHADER;
import static de.jreality.shader.CommonAttributes.POINT_RADIUS;
import static de.jreality.shader.CommonAttributes.POINT_SHADER;
import static de.jreality.shader.CommonAttributes.SPHERES_DRAW;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryVertices;
import static de.varylab.discreteconformal.adapter.HyperbolicModel.Klein;
import static de.varylab.discreteconformal.adapter.HyperbolicModel.Poincaré;
import static de.varylab.discreteconformal.plugin.TargetGeometry.Hyperbolic;
import static de.varylab.discreteconformal.uniformization.SurfaceCurveUtility.createIntersectionVertices;
import static de.varylab.discreteconformal.uniformization.SurfaceCurveUtility.createSurfaceCurves;
import static de.varylab.discreteconformal.unwrapper.BoundaryMode.ReadIsometricAngles;
import static de.varylab.discreteconformal.unwrapper.EuclideanLayout.calculateAngleSum;
import static de.varylab.discreteconformal.util.DiscreteEllipticUtility.calculateHalfPeriodRatio;
import static de.varylab.discreteconformal.util.UnwrapUtility.prepareInvariantDataEuclidean;
import static java.awt.Color.YELLOW;
import static javax.swing.JOptionPane.ERROR_MESSAGE;
import static javax.swing.JOptionPane.showMessageDialog;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import de.jreality.geometry.Primitives;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.plugin.JRViewer;
import de.jreality.plugin.basic.ViewShrinkPanelPlugin;
import de.jreality.plugin.job.Job;
import de.jreality.plugin.job.JobListener;
import de.jreality.plugin.job.JobQueuePlugin;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedLineSet;
import de.jreality.scene.PointSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.ui.ColorChooseJButton.ColorChangedEvent;
import de.jreality.ui.ColorChooseJButton.ColorChangedListener;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.SelectionInterface;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.halfedgetools.selection.SelectionListener;
import de.jtem.halfedgetools.selection.TypedSelection;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.jrworkspace.plugin.sidecontainer.widget.ShrinkPanel;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.functional.EuclideanCyclicFunctional;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Lambda;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Phi;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomEdgeInfo;
import de.varylab.discreteconformal.heds.CustomVertexInfo;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.heds.adapter.MetricErrorAdapter;
import de.varylab.discreteconformal.uniformization.CanonicalFormUtility;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility;
import de.varylab.discreteconformal.uniformization.FundamentalVertex;
import de.varylab.discreteconformal.unwrapper.BoundaryMode;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.HyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.QuantizationMode;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CPhi;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class DiscreteConformalPlugin extends ViewShrinkPanelPlugin 
	implements ListSelectionListener, ChangeListener, ActionListener, SelectionListener, ColorChangedListener, JobListener {

	public static final Integer
		CHANNEL_BROKEN_TRIANGLES = 23435634,
		CHANNEL_BOUNDARY_CONDITION = 1236644;
	private Logger
		log = Logger.getLogger(DiscreteConformalPlugin.class.getName());
	
	// plug-in section ------------------ 
	private HalfedgeInterface
		hif = null;
	private JobQueuePlugin
		jobQueue = null;
	private ConformalDataPlugin
		conformalDataPlugin = null;
	private UniformizationDomainPlugin
		domainPlugin = null;
	
	// data section ---------------------
	private CoHDS
		surface = null,
		surfaceUnwrapped = null;
	private UnwrapJob
		unwrapJob = null;
	private TargetGeometry
		activeGeometry = TargetGeometry.Automatic;
	private List<CoVertex>
		customVertices = new LinkedList<CoVertex>();
	private List<CoEdge>
		customEdges = new LinkedList<CoEdge>();	
	private FundamentalPolygon
		cuttedPolygon = null,
		minimalPolygon = null,
		oppositePolygon = null,
		canonicalPolygon = null,
		keenPolygon = null;
	private CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;

	private MetricErrorAdapter
		metricErrorAdapter = new MetricErrorAdapter();
	public CoPositionAdapter
		positionAdapter = new CoPositionAdapter();
	protected CoTexturePositionAdapter
		texturePositionAdapter = new CoTexturePositionAdapter();
	
	private Appearance
		yellowPointsAppearance = new Appearance(),
		polygonCurvesAppearance = new Appearance(),
		axesCurvesAppearance = new Appearance();
	private SceneGraphComponent
		selectedCustomNodesRoot = new SceneGraphComponent("Selected Custom Nodes"),
		polygonCurvesRoot = new SceneGraphComponent("Polygon Curves"),
		axesCurvesRoot = new SceneGraphComponent("Axes Curves");
	
	// user interface section ------------
	private JButton
		moveToCenterButton = new JButton("Center Selected Vertex"),
		checkGaussBonnetBtn = new JButton("Check Gauß-Bonnet"),
		unwrapBtn = new JButton("Unwrap"),
		doLayoutButton = new JButton("Recalculate Layout"),
		extractCutPreparedButton = new JButton("Extract Cut-Prepared"),
		resetSurfaceButton = new JButton("Reset Surface");
	private JComboBox<TargetGeometry>
		targetGeometryCombo = new JComboBox<>(TargetGeometry.values());
	private JComboBox<CutStrategy>
		cutStrategy = new JComboBox<CutStrategy>(CutStrategy.values());
	private ShrinkPanel
		customNodePanel = new ShrinkPanel("Selected Nodes"),
		boundaryPanel = new ShrinkPanel("Boundary"),
		coneConfigPanel = new ShrinkPanel("Cones"),
		toolsPanel = new ShrinkPanel("Tools");
	private SpinnerNumberModel
		customThetaModel = new SpinnerNumberModel(360.0, 0.0, 10000.0, 1.0),
		customPhiModel = new SpinnerNumberModel(180.0, 0.0, 360.0, 1.0),
		numConesModel = new SpinnerNumberModel(0, 0, 100, 1),
		toleranceExpModel = new SpinnerNumberModel(-8, -30, -1, 1),
		maxIterationsModel = new SpinnerNumberModel(150, 1, 10000, 1),
		snapToleranceExpModel = new SpinnerNumberModel(-5, -20, 1, 1);
	private JSpinner
		customThetaSpinner = new JSpinner(customThetaModel),
		customPhiSpinner = new JSpinner(customPhiModel),
		numConesSpinner = new JSpinner(numConesModel),
		toleranceExpSpinner = new JSpinner(toleranceExpModel),
		maxIterationsSpinner = new JSpinner(maxIterationsModel),
		snapToleranceExpSpinner = new JSpinner(snapToleranceExpModel);
	private JCheckBox
		uniformizationChecker = new JCheckBox("Create Uniformization"), 
		expertChecker = new JCheckBox("Expert Mode"),
		circularEdgeChecker = new JCheckBox("Circular"), 
		useCustomThetaChecker = new JCheckBox("Theta"),
		useProjectiveTexture = new JCheckBox("Projective Texture", true),
		drawCurvesOnSurface = new JCheckBox("Draw Curves On Surface"),
		stereographicChecker = new JCheckBox("Stereographic Spherical Uniformization", true);
	private JComboBox<String>
		numericsCombo = new JComboBox<String>(new String[] {"Petsc/Tao Numerics", "Java/MTJ Numerics"});
	private JComboBox<QuantizationMode>
		conesQuantizationModeCombo = new JComboBox<QuantizationMode>(QuantizationMode.values()),
		boundaryQuantizationCombo = new JComboBox<QuantizationMode>(QuantizationMode.values()),
		customQuantizationCombo = new JComboBox<QuantizationMode>(QuantizationMode.values());
	private JComboBox<BoundaryMode>
		boundaryModeCombo = new JComboBox<BoundaryMode>(BoundaryMode.values()),
		customModeCombo = new JComboBox<BoundaryMode>(BoundaryMode.values());
		
	private JList<Node<?, ?, ?>>
		customNodesList = new JList<Node<?, ?, ?>>();
	private JScrollPane
		selectionScroller = new JScrollPane(customNodesList);
		
	public DiscreteConformalPlugin() {
		polygonCurvesAppearance.setAttribute(EDGE_DRAW, true);
		polygonCurvesRoot.setAppearance(polygonCurvesAppearance);
		axesCurvesAppearance.setAttribute(EDGE_DRAW, true);
		axesCurvesRoot.setAppearance(axesCurvesAppearance);
		
		yellowPointsAppearance.setAttribute(VERTEX_DRAW, true);
		yellowPointsAppearance.setAttribute(POINT_SHADER + "." + DIFFUSE_COLOR, YELLOW);
		yellowPointsAppearance.setAttribute(POINT_SHADER + "." + SPHERES_DRAW, true);
		yellowPointsAppearance.setAttribute(POINT_SHADER + "." + POINT_RADIUS, 0.1);
		selectedCustomNodesRoot.setAppearance(yellowPointsAppearance);
	}

	
	private void connectGUIListeners() {
		customNodesList.getSelectionModel().addListSelectionListener(this);
		customModeCombo.addActionListener(this);
		customQuantizationCombo.addActionListener(this);
		customNodePanel.setShrinked(true);
		useCustomThetaChecker.addActionListener(this);
		customThetaSpinner.addChangeListener(this);
		customPhiSpinner.addChangeListener(this);
		circularEdgeChecker.addActionListener(this);
		moveToCenterButton.addActionListener(this);
		drawCurvesOnSurface.addActionListener(this);
		extractCutPreparedButton.addActionListener(this);
		resetSurfaceButton.addActionListener(this);
		doLayoutButton.addActionListener(this);
		unwrapBtn.addActionListener(this);
		checkGaussBonnetBtn.addActionListener(this);
		useProjectiveTexture.addActionListener(this);
		expertChecker.addActionListener(this);
	}
	
	private void createLayout() {
		boolean expert = expertChecker.isSelected();
		shrinkPanel.removeAll();
		shrinkPanel.setLayout(new GridBagLayout());
		GridBagConstraints c1 = new GridBagConstraints();
		c1.insets = new Insets(1,3,1,3);
		c1.fill = GridBagConstraints.BOTH;
		c1.anchor = GridBagConstraints.WEST;
		c1.weightx = 1.0;
		c1.gridwidth = 1;
		GridBagConstraints c2 = new GridBagConstraints();
		c2.insets = new Insets(1,1,1,1);
		c2.fill = GridBagConstraints.BOTH;
		c2.anchor = GridBagConstraints.WEST;
		c2.weightx = 1.0;
		c2.gridwidth = GridBagConstraints.REMAINDER;
		shrinkPanel.add(expertChecker, c2);
		if (expert) {
			numericsCombo.setLightWeightPopupEnabled(true);
			numericsCombo.setSelectedIndex(0);
			shrinkPanel.add(numericsCombo, c2);
			shrinkPanel.add(new JLabel("Tolerance Exp"), c1);
			shrinkPanel.add(toleranceExpSpinner, c2);
			shrinkPanel.add(new JLabel("Max Iterations"), c1);
			shrinkPanel.add(maxIterationsSpinner, c2);
			shrinkPanel.add(new JLabel("Cut Strategy"), c1);
			shrinkPanel.add(cutStrategy, c2);
			shrinkPanel.add(new JLabel("Target Geometry"), c1);
			shrinkPanel.add(targetGeometryCombo, c2);
			shrinkPanel.add(uniformizationChecker, c2);
			shrinkPanel.add(stereographicChecker, c2);
		}
		shrinkPanel.add(unwrapBtn, c2);
		if (expert) {
			shrinkPanel.add(checkGaussBonnetBtn, c2);
			shrinkPanel.add(doLayoutButton, c2);
			shrinkPanel.add(resetSurfaceButton, c2);
		}
		boundaryPanel.setLayout(new GridBagLayout());
		boundaryPanel.removeAll();
		boundaryPanel.add(new JLabel("Mode"), c1);
		boundaryPanel.add(boundaryModeCombo, c2);
		boundaryPanel.add(new JLabel("Quantization"), c1);
		boundaryPanel.add(boundaryQuantizationCombo, c2);
		boundaryPanel.setShrinked(true);
		shrinkPanel.add(boundaryPanel, c2);
		if (expert) {
			coneConfigPanel.setLayout(new GridBagLayout());
			coneConfigPanel.removeAll();
			coneConfigPanel.add(new JLabel("Max"), c1);
			coneConfigPanel.add(numConesSpinner, c2);
			coneConfigPanel.add(new JLabel("Quantization"), c1);
			coneConfigPanel.add(conesQuantizationModeCombo, c2);
			coneConfigPanel.setShrinked(true);
			shrinkPanel.add(coneConfigPanel, c2);
		}
		customNodePanel.setLayout(new GridBagLayout());
		customNodePanel.removeAll();
		customNodePanel.add(selectionScroller, c2);
		selectionScroller.setPreferredSize(new Dimension(10, 120));
		selectionScroller.setMinimumSize(new Dimension(10, 120));
		customNodePanel.add(useCustomThetaChecker, c1);
		customNodePanel.add(customThetaSpinner, c2);
		if (expert) {
			customNodePanel.add(new JLabel("Mode"), c1);
			customNodePanel.add(customModeCombo, c2);
			customNodePanel.add(new JLabel("Quantization"), c1);
			customNodePanel.add(customQuantizationCombo, c2);
			customNodePanel.add(circularEdgeChecker, c1);
			customNodePanel.add(customPhiSpinner, c2);
		}
		c2.weighty = 1.0;
		shrinkPanel.add(customNodePanel, c2);
		c2.weighty = 0.0;
		if (expert) {
			toolsPanel.setLayout(new GridBagLayout());
			toolsPanel.removeAll();
			toolsPanel.add(drawCurvesOnSurface, c2);
			toolsPanel.add(new JLabel("Snap Tolerance Exp"), c1);
			toolsPanel.add(snapToleranceExpSpinner, c2);
			toolsPanel.add(extractCutPreparedButton, c2);
			toolsPanel.add(moveToCenterButton, c2);
			toolsPanel.setShrinked(true);
			shrinkPanel.add(toolsPanel, c2);
		}
		shrinkPanel.revalidate();
	}
	
	public static void main(String[] args) {
		NativePathUtility.set("native");
		JRViewer v = new JRViewer();
		v.addBasicUI();
		v.addContentUI();
		v.registerPlugin(DiscreteConformalPlugin.class);
		v.setPropertiesFile("test.xml");
		v.startup();
	}
	
	
	@Override
	public void colorChanged(ColorChangedEvent cce) {
		updateUniformization(activeGeometry);
	}
	
	@Override
	public void selectionChanged(Selection s, HalfedgeInterface hif) {
		customVertices.clear();
		customEdges.clear();
		for (Vertex<?,?,?> v : s.getVertices()) {
			if (v instanceof CoVertex) {
				CoVertex cov = (CoVertex)v;
				if (cov.info == null) {
					cov.info = new CustomVertexInfo();
				}
				customVertices.add(cov);
			}
		}
		for (Edge<?,?,?> e : s.getEdges()) {
			if (e instanceof CoEdge) {
				CoEdge coe = (CoEdge)e;
				if (coe.info == null) {
					coe.info = new CustomEdgeInfo();
					coe.getOppositeEdge().info = coe.info;
				}
				customEdges.add(coe);
			}
		}
		DefaultListModel<Node<?, ?, ?>> model = new DefaultListModel<Node<?, ?, ?>>();
		for (CoVertex v : customVertices) {
			model.addElement(v);
		}
		for (CoEdge edge : customEdges) {
			if (edge.isPositive()) {
				model.addElement(edge);
			}
		}
		customNodesList.setModel(model);
	}

	
	@Override
	public void valueChanged(ListSelectionEvent e) {
		selectedCustomNodesRoot.removeAllChildren();
		if (customNodesList.getSelectedValue() == null) return;
		Object val = customNodesList.getSelectedValue();
		if (val instanceof CoVertex) {
			CoVertex v = (CoVertex)val;
			customModeCombo.setSelectedItem(v.info.boundaryMode);
			customQuantizationCombo.setSelectedItem(v.info.quantizationMode);
			useCustomThetaChecker.setSelected(v.info.useCustomTheta);
			customThetaModel.setValue(Math.toDegrees(v.info.theta));
			PointSet ps = Primitives.point(v.P);
			SceneGraphComponent psgc = new SceneGraphComponent("Vertex " + v.getIndex());
			psgc.setGeometry(ps);
			selectedCustomNodesRoot.addChild(psgc);
		}
		if (val instanceof CoEdge) {
			CoEdge edge = (CoEdge)val;
			circularEdgeChecker.setSelected(edge.info.circularHoleEdge);
			customPhiModel.setValue(Math.toDegrees(edge.info.phi));
		}
		hif.removeTemporaryGeometry(selectedCustomNodesRoot);
		hif.addTemporaryGeometry(selectedCustomNodesRoot);
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (customThetaSpinner == e.getSource()) {
			if (customNodesList.getSelectedValue() == null) return;
			for (Node<?,?,?> s : customNodesList.getSelectedValuesList()) {
				if (!(s instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)s;
				double thetaDeg = customThetaModel.getNumber().doubleValue();
				v.info.theta = Math.toRadians(thetaDeg);
			}
		}
		if (customPhiSpinner == e.getSource()) {
			if (customNodesList.getSelectedValue() == null) return;
			for (Node<?,?,?> s : customNodesList.getSelectedValuesList()) {
				if (!(s instanceof CoEdge)) continue;
				CoEdge edge = (CoEdge)s;
				double phiDeg = customPhiModel.getNumber().doubleValue();
				edge.info.phi = Math.toRadians(phiDeg);
				edge.getOppositeEdge().info.phi = Math.toRadians(phiDeg);
			}
		}
	}
	
	protected void updateGeometry(TargetGeometry target) {
		if (surfaceUnwrapped == null) return;
		if (target != Hyperbolic) {
			domainPlugin.setHyperbolicModel(Klein);
		} else {
			domainPlugin.setHyperbolicModel(Poincaré);
		}
		hif.addLayerAdapter(metricErrorAdapter, false);
		hif.set(surfaceUnwrapped);
	}
	
	@Override
	public void jobFinished(Job job) {
		unwrapJob = (UnwrapJob)job;
		surfaceUnwrapped = unwrapJob.getSurface();
		surfaceUnwrapped.revertNormalization();
		cutInfo = unwrapJob.getCutInfo();
		activeGeometry = unwrapJob.getTargetGeometry();
		metricErrorAdapter.setLengthMap(unwrapJob.getLengthMap());
		metricErrorAdapter.setSignature(Pn.EUCLIDEAN);
		conformalDataPlugin.addHalfedgeMap("Uniformizing Map", surfaceUnwrapped, cutInfo);
		if (uniformizationChecker.isSelected()) {
			createUniformization(surfaceUnwrapped, activeGeometry, cutInfo);
		} else {
			updateGeometry(activeGeometry);
		}
		if (unwrapJob.getBoundaryMode() == ReadIsometricAngles) {
			extractBoundaryAngles();
		}
	}

	protected void extractBoundaryAngles() {
		Selection boudaryVertices = new Selection();
		for (CoVertex v : boundaryVertices(surfaceUnwrapped)) {
			double angle = calculateAngleSum(v);
			v.info = new CustomVertexInfo();
			v.info.theta = angle;
			v.info.useCustomTheta = true;
			boudaryVertices.add(v, CHANNEL_BOUNDARY_CONDITION);
		}
		hif.addSelection(boudaryVertices);
	}
	
	@Override
	public void jobFailed(Job job, Exception e) {
		UnwrapJob unwrapper = (UnwrapJob)job;
		unwrapper.getSurface().revertNormalization();
		surfaceUnwrapped = null;
		Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
		JOptionPane.showMessageDialog(w, e, "Optimization Error", JOptionPane.WARNING_MESSAGE);
		e.printStackTrace();
	}
	@Override
	public void jobCancelled(Job job) {
	}
	@Override
	public void jobProgress(Job job, double progress) {
	}
	@Override
	public void jobStarted(Job job) {
	}
	
	public void createUniformization(CoHDS surfaceUnwrapped, TargetGeometry targetGeometry, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		this.surfaceUnwrapped = surfaceUnwrapped;
		this.activeGeometry = targetGeometry;
		this.cutInfo = cutInfo;
		if (
			(targetGeometry == TargetGeometry.Euclidean || 
			targetGeometry == TargetGeometry.Hyperbolic) &&
			!cutInfo.getBranchSet().isEmpty()
		) {
			int signature = targetGeometry.getSignature();
			try {
				System.out.println("Constructing fundamental cut polygon...");
				cuttedPolygon = FundamentalPolygonUtility.constructFundamentalPolygon(cutInfo, signature);
				conformalDataPlugin.addUniformizationData("Direct Uniformization", cuttedPolygon);
				System.out.println(cuttedPolygon);
				FundamentalVertex root = cuttedPolygon.getEdges().get(0).start;
				System.out.println("Constructing minimal polygon...");
				minimalPolygon = FundamentalPolygonUtility.minimize(cuttedPolygon, root);
				conformalDataPlugin.addUniformizationData("Minimal Uniformization", minimalPolygon);
				System.out.println(minimalPolygon);
				minimalPolygon.checkRelation();
				if (targetGeometry == TargetGeometry.Hyperbolic) {
					System.out.println("Constructing opposites sides polygon...");
					oppositePolygon = CanonicalFormUtility.canonicalizeOpposite(minimalPolygon);
					System.out.println(oppositePolygon);
					conformalDataPlugin.addUniformizationData("Opposite Edges Uniformization", oppositePolygon);
					oppositePolygon.checkRelation();	
					System.out.println("Constructing fast canonical polygon...");
					canonicalPolygon = FundamentalPolygonUtility.canonicalize(minimalPolygon, false);
					System.out.println(canonicalPolygon);
					conformalDataPlugin.addUniformizationData("Canonical Uniformization", canonicalPolygon);
					canonicalPolygon.checkRelation();
					System.out.println("Constructing Keen polygon...");
					keenPolygon = CanonicalFormUtility.createKeenPolygon(canonicalPolygon);
					System.out.println(keenPolygon);
					conformalDataPlugin.addUniformizationData("Keen Uniformization", keenPolygon);
					keenPolygon.checkRelation();
				} else {
					oppositePolygon = minimalPolygon;
					canonicalPolygon = minimalPolygon;
					keenPolygon = minimalPolygon;
				}
				metricErrorAdapter.setSignature(signature);
			} catch (Exception e) {
				log.log(Level.WARNING, e.getMessage(), e);
			}
			if (targetGeometry == TargetGeometry.Euclidean) {
				Complex tau = calculateHalfPeriodRatio(cutInfo);
				log.info("torus modulus: " + tau);
			}
		} else {
			cutInfo = null;
			cuttedPolygon = null;
			minimalPolygon = null;
			oppositePolygon = null;
			canonicalPolygon = null;
			keenPolygon = null;
		}
		updateGeometry(targetGeometry);
		updateUniformization(targetGeometry);
	}

	protected void updateUniformization(final TargetGeometry target) {
		Runnable r = new Runnable() {
			@Override
			public void run() {
				updateUniformizationDirect(target);
			}
		};
		EventQueue.invokeLater(r);
	}
	
	protected void updateUniformizationDirect(TargetGeometry target) {
		int signature = target.getSignature();
		AdapterSet aSet = hif.getActiveAdapters();
		int maxGroupElements = domainPlugin.getMaxCoverElements();
		double maxDrawDistance = domainPlugin.getMaxCoverDistance();
		boolean drawCurves = drawCurvesOnSurface.isSelected();
		hif.removeTemporaryGeometry(polygonCurvesRoot);
		hif.removeTemporaryGeometry(axesCurvesRoot);
		FundamentalPolygon P = getActiveFundamentalPoygon();
		if (P == null) return;
		if (target == TargetGeometry.Hyperbolic && drawCurves) {
			IndexedLineSet curves = createSurfaceCurves(P, surfaceUnwrapped, aSet, maxGroupElements, maxDrawDistance, true, false, signature);
			polygonCurvesAppearance.setAttribute(LINE_SHADER + "." + DIFFUSE_COLOR, Color.RED);
			polygonCurvesRoot.setGeometry(curves);
			hif.addTemporaryGeometry(polygonCurvesRoot);
			IndexedLineSet axes = createSurfaceCurves(P, surfaceUnwrapped, aSet, maxGroupElements, maxDrawDistance, false, true, signature);
			axesCurvesAppearance.setAttribute(LINE_SHADER + "." + DIFFUSE_COLOR, Color.BLUE);
			axesCurvesRoot.setGeometry(axes);
			hif.addTemporaryGeometry(axesCurvesRoot);
		}
		if (domainPlugin != null) {
			domainPlugin.setGeometry(target);
			domainPlugin.createUniformization(surfaceUnwrapped, cuttedPolygon, minimalPolygon, canonicalPolygon, oppositePolygon, keenPolygon);
		}
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Window w = SwingUtilities.getWindowAncestor(this.shrinkPanel);
		Object s = e.getSource();
		if (unwrapBtn == s) {
			CoHDS unwrapped = getLoaderGeometry();
			if (unwrapped == null) return;
			surface = copySurface(unwrapped);
			Selection selection = hif.getSelection();
			AdapterSet aSet = hif.getAdapters();
			conformalDataPlugin.addDiscreteMetric("Input Discrete Metric", unwrapped, aSet);
			conformalDataPlugin.addHalfedgeEmbedding("Input Discrete Position Embedding", unwrapped, selection, aSet, Position4d.class, null);
			unwrapped.normalizeCoordinates();
			UnwrapJob uw = new UnwrapJob(unwrapped, aSet);
			uw.setTargetGeometry((TargetGeometry)targetGeometryCombo.getSelectedItem());
			uw.setToleranceExponent(toleranceExpModel.getNumber().intValue());
			uw.setMaxIterations(maxIterationsModel.getNumber().intValue());
			uw.setNumCones(numConesModel.getNumber().intValue());
			uw.setQuantizationMode((QuantizationMode)conesQuantizationModeCombo.getSelectedItem());
			uw.setBoundaryQuantizationMode((QuantizationMode)boundaryQuantizationCombo.getSelectedItem());
			uw.setBoundaryMode((BoundaryMode)boundaryModeCombo.getSelectedItem());
			uw.setUsePetsc(numericsCombo.getSelectedIndex() == 0);
			uw.setSelectedVertices(selection.getVertices(unwrapped));
			uw.setSelectedEdges(selection.getEdges(unwrapped));
			uw.setCutStrategy((CutStrategy)cutStrategy.getSelectedItem());
			uw.setUseStereographicSphericalUniformization(stereographicChecker.isSelected());
			uw.addJobListener(this);
			jobQueue.queueJob(uw);
		}
		if (doLayoutButton == s) {
			if (surfaceUnwrapped == null) return;
			for (CoEdge ee : surfaceUnwrapped.getEdges()) {
				// fill missing edge lengths
				if (unwrapJob.getLengthMap().containsKey(ee)) {
					unwrapJob.getLengthMap().put(ee.getOppositeEdge(), unwrapJob.getLengthMap().get(ee));
				}
			}
			switch (activeGeometry) {
			case Euclidean:
				EuclideanLayout.doLayout(surfaceUnwrapped, unwrapJob.getLengthMap(), unwrapJob.getAngleMap());
				break;
			case Hyperbolic:
				HyperbolicLayout.doLayout(surfaceUnwrapped, surfaceUnwrapped.getVertex(0), unwrapJob.getLengthMap());
				break;
			default:
				log.warning("recalculation of layout not implemented for " + activeGeometry + " geometry");
				break;
			}
			hif.update();
		}
		if (customModeCombo == s) {
			for (Node<?,?,?> sel : customNodesList.getSelectedValuesList()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.boundaryMode = (BoundaryMode)customModeCombo.getSelectedItem();
			}
		}
		if (customQuantizationCombo == s) {
			for (Node<?,?,?> sel : customNodesList.getSelectedValuesList()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.quantizationMode = (QuantizationMode)customQuantizationCombo.getSelectedItem();
			}
		}
		if (useCustomThetaChecker == s) {
			for (Node<?,?,?> sel : customNodesList.getSelectedValuesList()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.useCustomTheta = useCustomThetaChecker.isSelected(); 
			}
		}
		if (circularEdgeChecker == s) {
			for (Node<?,?,?> sel : customNodesList.getSelectedValuesList()) {
				if (!(sel instanceof CoEdge)) continue;
				CoEdge edge = (CoEdge)sel;
				edge.info.circularHoleEdge = circularEdgeChecker.isSelected();
			}
		}
		if (checkGaussBonnetBtn == s) {
			CoHDS hds = hif.get(new CoHDS());
			BoundaryMode boundaryMode = (BoundaryMode)boundaryModeCombo.getSelectedItem();
			QuantizationMode boundaryQuantMode = (QuantizationMode)boundaryQuantizationCombo.getSelectedItem();
			try {
				Theta<CoVertex, CoEdge> theta = new CTheta();
				Phi<CoEdge> phi = new CPhi();
				Variable<CoVertex, CoEdge> variable = new CVariable();
				Lambda<CoEdge> lambda = new CLambda();
				Alpha<CoEdge> alpha = new CAlpha();
				InitialEnergy<CoFace> initE = new CInitialEnergy();
				EuclideanCyclicFunctional<CoVertex, CoEdge, CoFace> fun = new EuclideanCyclicFunctional<CoVertex, CoEdge, CoFace>(variable, theta, phi, lambda, alpha, initE);
				prepareInvariantDataEuclidean(fun, hds, boundaryMode, boundaryQuantMode, hif.getAdapters());
			} catch (Exception e1) {
				showMessageDialog(w, e1.getMessage(), "Error", ERROR_MESSAGE);
			}
		}
		if (moveToCenterButton == s) {
			if (surfaceUnwrapped == null) return;
			TypedSelection<CoVertex> sel = hif.getSelection().getVertices(surfaceUnwrapped);
			if (sel.isEmpty()) return;
			CoVertex v = sel.iterator().next();
			double[] pos = v.T;
			MatrixBuilder mb = MatrixBuilder.hyperbolic();
			mb.translateFromTo(pos, new double[] {0,0,0,1});
			Matrix T = mb.getMatrix();
			for (CoVertex vv : surfaceUnwrapped.getVertices()) {
				T.transformVector(vv.T);
			}
			hif.update();
			createUniformization(surfaceUnwrapped, activeGeometry, cutInfo);
			conformalDataPlugin.addHalfedgeMap("Uniformizing Map", surfaceUnwrapped, cutInfo);
		}
		if (drawCurvesOnSurface == s) {
			updateGeometry(activeGeometry);
			updateUniformization(activeGeometry);				
		}
		if (expertChecker == s) {
			createLayout();
		}
		if (extractCutPreparedButton == s) {
			if (getActiveFundamentalPoygon() == null) {
				showMessageDialog(w, "No uniformization has been calculated.");
				return;
			}
			int signature = activeGeometry.getSignature();
			AdapterSet a = hif.getAdapters();
			CoHDS intersected = copySurface(surface);
			double snapTolerance = Math.pow(10, snapToleranceExpModel.getNumber().intValue());
			Set<Set<CoVertex>> intersectingVertices = createIntersectionVertices(getActiveFundamentalPoygon(), intersected, surfaceUnwrapped, cutInfo, a, snapTolerance, signature);
			Selection pathSelection = new Selection();
			for (Set<CoVertex> segmentSet : intersectingVertices) {
				for (CoVertex v : segmentSet) {
					for (CoFace f : HalfEdgeUtilsExtra.getFaceStar(v)) {
						for (CoVertex nextOnSegment : HalfEdgeUtils.boundaryVertices(f)) {
							List<CoVertex> vStar = HalfEdgeUtilsExtra.getVertexStar(nextOnSegment);
							if (nextOnSegment == v || vStar.contains(v)) {
								if (vStar.contains(v) && segmentSet.contains(nextOnSegment)) {
									CoEdge pathedge = HalfEdgeUtils.findEdgeBetweenVertices(v, nextOnSegment);
									if (pathedge != null) {
										pathSelection.add(pathedge);
										pathSelection.add(pathedge.getOppositeEdge());
									}
								}
								continue;
							}
							if (segmentSet.contains(nextOnSegment)) {
								CoEdge segmentEdge = TopologyAlgorithms.splitFaceAt(f, v, nextOnSegment);
								pathSelection.add(segmentEdge);
								pathSelection.add(segmentEdge.getOppositeEdge());
							}
						}
					}
				}
			}
			Triangulator.triangulateByCuttingCorners(intersected, a);
			HalfedgeLayer l = new HalfedgeLayer(hif);
			hif.addLayer(l);
			l.set(intersected);
			l.setSelection(pathSelection);
		}
		if (resetSurfaceButton == s && surface != null) {
			hif.set(surface);
		}
	}

	public static LinkedList<CoEdge> selectCutPath(CoHDS intersected,	Set<CoVertex> newVertices, Selection pathSelection) {
		LinkedList<CoEdge> newCut = new LinkedList<CoEdge>();
		for (CoEdge edge : intersected.getEdges()) {
			if (newVertices.contains(edge.getTargetVertex()) && newVertices.contains(edge.getStartVertex())) {
				pathSelection.add(edge);	
				newCut.add(edge);
			} 
		}
		return newCut;
	}


	private FundamentalPolygon getActiveFundamentalPoygon() {
		switch (domainPlugin.getSelectedPoygonType()) {
			default:
			case Cut:
				return cuttedPolygon;
			case Minimal:
				return minimalPolygon;
			case Canonical:
				return canonicalPolygon;
			case Opposite:
				return oppositePolygon;
			case Keen:
				return keenPolygon;
		}
	}
	
	public TargetGeometry getActiveGeometry() {
		return activeGeometry;
	}
	
	protected CoHDS getLoaderGeometry() {
		CoHDS surface = new CoHDS();
		surface.setTexCoordinatesValid(false);
		surface = hif.get(surface);
		if (surface.numVertices() == 0) {
			return null;
		}
		// clear non-custom vertices
		for (CoVertex v : surface.getVertices()) {
			if (!customVertices.contains(v)) {
				v.info = null;
			}
		}
		for (CoEdge e : surface.getEdges()) {
			if (!customEdges.contains(e) && !customEdges.contains(e.getOppositeEdge())) {
				e.info = null;
			}
		}
		Triangulator.triangulateSingleSource(surface);
		return surface;
	}
	
	protected void setSurface(CoHDS hds) {
		this.surface = hds;
	}
	
	
	private CoHDS copySurface(CoHDS hds) {
		CoHDS copy = new CoHDS();
		HalfEdgeUtils.copy(hds, copy);
		for (CoVertex v : hds.getVertices()) {
			CoVertex cv = copy.getVertex(v.getIndex());
			cv.copyData(v);
		}
		for (CoEdge e : hds.getEdges()) {
			CoEdge ce = copy.getEdge(e.getIndex());
			ce.copyData(e);
		}
		for (CoFace f : hds.getFaces()){
			CoFace cf = copy.getFace(f.getIndex());
			cf.copyData(f);
		}
		return copy;
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		hif.addAdapter(positionAdapter, true);
		hif.addAdapter(texturePositionAdapter, true);
		hif.addSelectionListener(this);
		jobQueue = c.getPlugin(JobQueuePlugin.class);
		conformalDataPlugin = c.getPlugin(ConformalDataPlugin.class);
		domainPlugin = c.getPlugin(UniformizationDomainPlugin.class);
		SelectionInterface sif = c.getPlugin(SelectionInterface.class);
		sif.registerChannelName(CHANNEL_BROKEN_TRIANGLES, "Broken Triangles");
		sif.registerChannelName(CHANNEL_BOUNDARY_CONDITION, "Boundary Conditions");
		createLayout();
		connectGUIListeners();
	}

	
	@Override
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "expertMode", expertChecker.isSelected());
		c.storeProperty(getClass(), "numCones", numConesModel.getNumber().intValue());
		c.storeProperty(getClass(), "conesQuantizationModeIndex", conesQuantizationModeCombo.getSelectedIndex());
		c.storeProperty(getClass(), "boundaryModeIndex", boundaryModeCombo.getSelectedIndex());
		c.storeProperty(getClass(), "boundaryQuantModeIndex", boundaryQuantizationCombo.getSelectedIndex());
		c.storeProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex());
		c.storeProperty(getClass(), "useProjectiveTexture", useProjectiveTexture.isSelected());
		c.storeProperty(getClass(), "toleranceExponent", toleranceExpModel.getNumber());
		c.storeProperty(getClass(), "maxIterations", maxIterationsModel.getNumber());
		c.storeProperty(getClass(), "boundaryPanelShrinked", boundaryPanel.isShrinked());	
		c.storeProperty(getClass(), "conesPanelShrinked", coneConfigPanel.isShrinked());	
		c.storeProperty(getClass(), "modelPanelShrinked", toolsPanel.isShrinked());	
		c.storeProperty(getClass(), "customVertexPanelShrinked", customNodePanel.isShrinked());
		c.storeProperty(getClass(), "doUniformization", uniformizationChecker.isSelected());
	} 
	
 
	@Override
	public void restoreStates(Controller c) throws Exception { 
		super.restoreStates(c);
		expertChecker.setSelected(c.getProperty(getClass(), "expertMode", expertChecker.isSelected()));
		numConesModel.setValue(c.getProperty(getClass(), "numCones", numConesModel.getNumber().intValue()));
		numericsCombo.setSelectedIndex(c.getProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex()));
		useProjectiveTexture.setSelected(c.getProperty(getClass(), "useProjectiveTexture", useProjectiveTexture.isSelected()));
		toleranceExpModel.setValue(c.getProperty(getClass(), "toleranceExponent", toleranceExpModel.getNumber()));
		maxIterationsModel.setValue(c.getProperty(getClass(), "maxIterations", maxIterationsModel.getNumber()));
		conesQuantizationModeCombo.setSelectedIndex(c.getProperty(getClass(), "conesQuantizationModeIndex", conesQuantizationModeCombo.getSelectedIndex()));
		boundaryModeCombo.setSelectedIndex(c.getProperty(getClass(), "boundaryModeIndex", boundaryModeCombo.getSelectedIndex()));
		boundaryQuantizationCombo.setSelectedIndex(c.getProperty(getClass(), "boundaryQuantModeIndex", boundaryQuantizationCombo.getSelectedIndex()));
		boundaryPanel.setShrinked(c.getProperty(getClass(), "boundaryPanelShrinked", true));
		coneConfigPanel.setShrinked(c.getProperty(getClass(), "conesPanelShrinked", true));
		toolsPanel.setShrinked(c.getProperty(getClass(), "modelPanelShrinked", true));
		customNodePanel.setShrinked(c.getProperty(getClass(), "customVertexPanelShrinked", customNodePanel.isShrinked()));
		uniformizationChecker.setSelected(c.getProperty(getClass(), "doUniformization", false));
	}
	
	
	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = new PluginInfo();
		info.name = "Discrete Conformal Parametrization";
		info.vendorName = "Stefan Sechelmann";
		info.email = "sechel@math.tu-berlin.de";
		return info;
	}
	
	public CuttingInfo<CoVertex, CoEdge, CoFace> getCurrentCutInfo() {
		return cutInfo;
	}
	public void setCutCurrentInfo(CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		this.cutInfo = cutInfo;
	}

}