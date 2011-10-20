package de.varylab.discreteconformal.plugin;

import static de.jreality.math.MatrixBuilder.euclidean;
import static de.jreality.scene.Appearance.DEFAULT;
import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.LIGHTING_ENABLED;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.TEXTURE_2D;
import static de.jreality.shader.CommonAttributes.TRANSPARENCY;
import static de.jreality.shader.CommonAttributes.TRANSPARENCY_ENABLED;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.varylab.discreteconformal.util.UnwrapUtility.prepareInvariantDataEuclidean;
import static java.awt.Color.BLACK;
import static java.awt.Color.WHITE;
import static java.awt.GridBagConstraints.RELATIVE;
import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static javax.swing.JOptionPane.ERROR_MESSAGE;
import static javax.swing.JOptionPane.OK_CANCEL_OPTION;
import static javax.swing.JOptionPane.showMessageDialog;
import static javax.swing.SwingUtilities.getWindowAncestor;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;

import javax.imageio.ImageIO;
import javax.swing.ButtonGroup;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import de.jreality.geometry.IndexedFaceSetFactory;
import de.jreality.geometry.Primitives;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.content.ContentAppearance;
import de.jreality.scene.Appearance;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.shader.ImageData;
import de.jreality.shader.Texture2D;
import de.jreality.shader.TextureUtility;
import de.jreality.ui.AppearanceInspector;
import de.jreality.ui.TextureInspector;
import de.jreality.ui.viewerapp.FileFilter;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeSelection;
import de.jtem.halfedgetools.plugin.SelectionListener;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.jrworkspace.plugin.sidecontainer.widget.ShrinkPanel;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.functional.EuclideanFunctional;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Lambda;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomEdgeInfo;
import de.varylab.discreteconformal.heds.CustomVertexInfo;
import de.varylab.discreteconformal.heds.adapter.BranchPointColorAdapter;
import de.varylab.discreteconformal.heds.adapter.BranchPointRadiusAdapter;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.MarkedEdgesColorAdapter;
import de.varylab.discreteconformal.heds.adapter.MarkedEdgesRadiusAdapter;
import de.varylab.discreteconformal.heds.adapter.MetricErrorAdapter;
import de.varylab.discreteconformal.plugin.tasks.Unwrap;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility;
import de.varylab.discreteconformal.uniformization.FundamentalVertex;
import de.varylab.discreteconformal.uniformization.VisualizationUtility;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.util.UnwrapUtility.QuantizationMode;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin implements ListSelectionListener, ChangeListener, ActionListener, PropertyChangeListener, SelectionListener {

	private static int
		coverResolution = 1024;
	
	private enum Domain {
		Cut,
		Minimal,
		Canonical
	}
	
	// plug-in section ------------------ 
	private HalfedgeInterface
		hif = null;
	private ContentAppearance
		contentAppearance = null;
	
	// data section ---------------------
	private CoHDS
		surface = null;
	private List<CoVertex>
		customVertices = new LinkedList<CoVertex>();
	private List<CoEdge>
		customEdges = new LinkedList<CoEdge>();	
	private FundamentalPolygon 
		cuttedPolygon = null,
		minimalPolygon = null,
		canonicalPolygon = null;
	private Matrix 
		polygonTextureMatrix = euclidean().translate(-0.5, -0.5, 0).scale(0.5).scale(1, -1, 1).getMatrix();
	private BufferedImage
		cutCoverImage = null,
		minimalCoverImage = null,
		canonicalCoverImage = null;
	private int
		genus = -1;
	private CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;

	private MetricErrorAdapter
		metricErrorAdapter = new MetricErrorAdapter();
	private CoTexturePositionAdapter
		texturePositionAdapter = new CoTexturePositionAdapter();
	private CoTexturePositionPositionAdapter
		texCoordPositionAdapter = new CoTexturePositionPositionAdapter();
	private MarkedEdgesColorAdapter
		cutColorAdapter = new MarkedEdgesColorAdapter();
	private MarkedEdgesRadiusAdapter
		cutRadiusAdapter = new MarkedEdgesRadiusAdapter();
	private BranchPointColorAdapter
		pointColorAdapter = new BranchPointColorAdapter();
	private BranchPointRadiusAdapter
		pointRadiusAdapter = new BranchPointRadiusAdapter();
	
//	private Tool
//		hyperbolicCopyTool = new HyperbolicCopyTool(this);
	private Appearance
		universalCoverAppearance = new Appearance();
	private SceneGraphComponent
		unitCircle = new SceneGraphComponent("Hyperbolic Boundary"),
		cutCoverRoot = new SceneGraphComponent("Cut Cover"),
		minimalCoverRoot = new SceneGraphComponent("Minimal Cover"),
		canonicalCoverRoot = new SceneGraphComponent("Canonical Cover");

	// user interface section ------------
	private JButton
		saveTextureButton = new JButton("Save Texture"),
		resetXFormButton = new JButton("Reset Transformation"),
		moveToCenterButton = new JButton("Center Selected Vertex"),
		coverToTextureButton = new JButton("Create Texture"),
		checkGaussBonnetBtn = new JButton("Check Gauß-Bonnet"),
		unwrapBtn = new JButton("Unwrap"),
		quantizeToQuads = new JButton("Quads");
	private JComboBox
		domainCombo = new JComboBox(Domain.values());
	private ShrinkPanel
		customNodePanel = new ShrinkPanel("Custom Vertices"),
		boundaryPanel = new ShrinkPanel("Boundary"),
		coneConfigPanel = new ShrinkPanel("Automatic Cones"),
		modelPanel = new ShrinkPanel("Hyperbolic Model"),
		visualizationPanel = new ShrinkPanel("Visualization"),
		texQuantizationPanel = new ShrinkPanel("Cone Texture Quantization");
	private SpinnerNumberModel
		coverRecursionModel = new SpinnerNumberModel(2, 0, 10, 1),
		customThetaModel = new SpinnerNumberModel(360.0, 0.0, 1000.0, 1.0),
		numConesModel = new SpinnerNumberModel(0, 0, 100, 1),
		toleranceExpModel = new SpinnerNumberModel(-8, -30, -1, 1),
		maxIterationsModel = new SpinnerNumberModel(150, 1, 10000, 1);
	private JSpinner
		coverRecursionSpinner = new JSpinner(coverRecursionModel),
		customThetaSpinner = new JSpinner(customThetaModel),
		numConesSpinner = new JSpinner(numConesModel),
		toleranceExpSpinner = new JSpinner(toleranceExpModel),
		maxIterationsSpinner = new JSpinner(maxIterationsModel);
	private JCheckBox
		circularEdgeChecker = new JCheckBox("Is Circular Edge"), 
		useDistanceToCanonicalize = new JCheckBox("Use Isometry Distances"),
		useCustomThetaChecker = new JCheckBox("Custom Theta"),
		useProjectiveTexture = new JCheckBox("Projective Texture", true),
		showUnwrapped = new JCheckBox("Unwrapped");
	private JComboBox
		numericsCombo = new JComboBox(new String[] {"Java/MTJ Numerics", "Petsc/Tao Numerics"}),
		quantizationModeCombo = new JComboBox(QuantizationMode.values()),
		boundaryModeCombo = new JComboBox(BoundaryMode.values()),
		boundaryQuantizationCombo = new JComboBox(QuantizationMode.values()),
		customModeCombo = new JComboBox(BoundaryMode.values()),
		customQuantizationCombo = new JComboBox(QuantizationMode.values());
	private JRadioButton
		kleinButton = new JRadioButton("Klein"),
		poincareButton = new JRadioButton("Poincaré", true),
		halfplaneButton = new JRadioButton("Half-Plane"); 
	private JList
		selectedNodesList = new JList();
	private JScrollPane
		selectionScroller = new JScrollPane(selectedNodesList);
		
	public DiscreteConformalPlugin() {
		createLayout();
		unwrapBtn.addActionListener(this);
		checkGaussBonnetBtn.addActionListener(this);
		showUnwrapped.addActionListener(this);
		kleinButton.addActionListener(this);
		poincareButton.addActionListener(this);
		halfplaneButton.addActionListener(this);
		domainCombo.addActionListener(this);
		useProjectiveTexture.addActionListener(this);
		coverToTextureButton.addActionListener(this);
		quantizeToQuads.addActionListener(this);
		
		ButtonGroup modelGroup = new ButtonGroup();
		modelGroup.add(kleinButton);
		modelGroup.add(poincareButton);
		modelGroup.add(halfplaneButton);
		
		IndexedFaceSetFactory ifsf = new IndexedFaceSetFactory();
		ifsf.setVertexCount(4);
		ifsf.setFaceCount(1);
		ifsf.setFaceIndices(new int[][] {{0,1,2,3}});
		ifsf.setVertexTextureCoordinates(new double[] {-1,-1,1,-1,1,1,-1,1});
		ifsf.setVertexCoordinates(new double[]{-1,-1,0.001, 1,-1,0.001, 1,1,0.001, -1,1,0.001});
		ifsf.setGenerateFaceNormals(true);
		ifsf.update();
		cutCoverRoot.setGeometry(ifsf.getGeometry());
		cutCoverRoot.setAppearance(universalCoverAppearance);
		minimalCoverRoot.setGeometry(ifsf.getGeometry());
		minimalCoverRoot.setAppearance(universalCoverAppearance);
		canonicalCoverRoot.setGeometry(ifsf.getGeometry());
		canonicalCoverRoot.setAppearance(universalCoverAppearance);
		universalCoverAppearance.setAttribute(VERTEX_DRAW, false);
		universalCoverAppearance.setAttribute(EDGE_DRAW, false);
		universalCoverAppearance.setAttribute(FACE_DRAW, true);
		universalCoverAppearance.setAttribute(LIGHTING_ENABLED, false);
		universalCoverAppearance.setAttribute(TRANSPARENCY_ENABLED, true);
		universalCoverAppearance.setAttribute(POLYGON_SHADER + "." + DIFFUSE_COLOR, WHITE);
		universalCoverAppearance.setAttribute(POLYGON_SHADER + "." + TRANSPARENCY, 0.0);
		
		Appearance circleApp = new Appearance();
		circleApp.setAttribute(EDGE_DRAW, false);
		circleApp.setAttribute(VERTEX_DRAW, false); 
		circleApp.setAttribute(FACE_DRAW, true);
		circleApp.setAttribute(POLYGON_SHADER + "." + DIFFUSE_COLOR, BLACK);
		circleApp.setAttribute(POLYGON_SHADER + "." + TEXTURE_2D, DEFAULT); 
		unitCircle.setAppearance(circleApp);
		euclidean().rotate(PI / 2, 1, 0, 0).assignTo(unitCircle);
		unitCircle.setGeometry(Primitives.torus(1.0025, 0.005, 200, 5));
		cutCoverRoot.addChild(unitCircle);
		minimalCoverRoot.addChild(unitCircle);
		canonicalCoverRoot.addChild(unitCircle);
	}

	
	private void createLayout() {
		shrinkPanel.setLayout(new GridBagLayout());
		GridBagConstraints c1 = new GridBagConstraints();
		c1.insets = new Insets(1,1,1,1);
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
		
		numericsCombo.setLightWeightPopupEnabled(true);
		numericsCombo.setSelectedIndex(0);
		shrinkPanel.add(numericsCombo, c2);
		shrinkPanel.add(new JLabel("Tolerance Exp"), c1);
		shrinkPanel.add(toleranceExpSpinner, c2);
		shrinkPanel.add(new JLabel("Max Iterations"), c1);
		shrinkPanel.add(maxIterationsSpinner, c2);
		shrinkPanel.add(checkGaussBonnetBtn, c1);
		shrinkPanel.add(unwrapBtn, c2);
		
		boundaryPanel.setLayout(new GridBagLayout());
		boundaryPanel.add(new JLabel("Mode"), c1);
		boundaryPanel.add(boundaryModeCombo, c2);
		boundaryPanel.add(new JLabel("Quantization"), c1);
		boundaryPanel.add(boundaryQuantizationCombo, c2);
		boundaryPanel.setShrinked(true);
		shrinkPanel.add(boundaryPanel, c2);
		
		coneConfigPanel.setLayout(new GridBagLayout());
		coneConfigPanel.add(new JLabel("Cones"), c1);
		coneConfigPanel.add(numConesSpinner, c2);
		coneConfigPanel.add(new JLabel("Quantization"), c1);
		coneConfigPanel.add(quantizationModeCombo, c2);
		coneConfigPanel.setShrinked(true);
		shrinkPanel.add(coneConfigPanel, c2);
		
		customNodePanel.add(selectionScroller, c2);
		selectionScroller.setPreferredSize(new Dimension(10, 70));
		selectionScroller.setMinimumSize(new Dimension(10, 70));
		customNodePanel.add(useCustomThetaChecker, c1);
		customNodePanel.add(customThetaSpinner, c2);
		customNodePanel.add(new JLabel("Mode"), c1);
		customNodePanel.add(customModeCombo, c2);
		customNodePanel.add(new JLabel("Quantization"), c1);
		customNodePanel.add(customQuantizationCombo, c2);
		customNodePanel.add(circularEdgeChecker, c2);
		shrinkPanel.add(customNodePanel, c2);
		
		texQuantizationPanel.add(quantizeToQuads, c2);
		shrinkPanel.add(texQuantizationPanel, c2);
		
		visualizationPanel.setLayout(new GridBagLayout());
		visualizationPanel.add(showUnwrapped, c1);
		visualizationPanel.add(resetXFormButton, c2);
		visualizationPanel.add(new JLabel("Domain"), c1);
		visualizationPanel.add(domainCombo, c2);
		visualizationPanel.add(new JLabel("Cover Recursion"), c1);
		visualizationPanel.add(coverRecursionSpinner, c2);
		visualizationPanel.add(useDistanceToCanonicalize, c2);
		visualizationPanel.add(useProjectiveTexture, c1);
		visualizationPanel.add(coverToTextureButton, c2);
		visualizationPanel.add(saveTextureButton, c2);
		visualizationPanel.setShrinked(true);
		shrinkPanel.add(visualizationPanel, c2);
		
		modelPanel.setLayout(new GridBagLayout()); 	c1.gridwidth = 1;
		modelPanel.add(kleinButton, c1); 			c1.gridwidth = RELATIVE;
		modelPanel.add(poincareButton, c1);
		modelPanel.add(halfplaneButton, c2);
		modelPanel.add(moveToCenterButton, c2);
		modelPanel.setShrinked(true);
		shrinkPanel.add(modelPanel, c2);
		
		selectedNodesList.getSelectionModel().addListSelectionListener(this);
		customModeCombo.addActionListener(this);
		customQuantizationCombo.addActionListener(this);
		customNodePanel.setShrinked(true);
		useCustomThetaChecker.addActionListener(this);
		customThetaSpinner.addChangeListener(this);
		circularEdgeChecker.addActionListener(this);
		moveToCenterButton.addActionListener(this);
		resetXFormButton.addActionListener(this);
		saveTextureButton.addActionListener(this);
	}
	
	
	@Override
	public void selectionChanged(HalfedgeSelection s, HalfedgeInterface hif) {
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
		DefaultListModel model = new DefaultListModel();
		for (CoVertex v : customVertices) {
			model.addElement(v);
		}
		for (CoEdge edge : customEdges) {
			if (edge.isPositive()) {
				model.addElement(edge);
			}
		}
		selectedNodesList.setModel(model);
	}

	
	@Override
	public void valueChanged(ListSelectionEvent e) {
		if (selectedNodesList.getSelectedValue() == null) return;
		Object val = selectedNodesList.getSelectedValue();
		if (val instanceof CoVertex) {
			CoVertex v = (CoVertex)val;
			customModeCombo.setSelectedItem(v.info.boundaryMode);
			customQuantizationCombo.setSelectedItem(v.info.quantizationMode);
			useCustomThetaChecker.setSelected(v.info.useCustomTheta);
			customThetaModel.setValue(Math.toDegrees(v.info.theta));
		}
		if (val instanceof CoEdge) {
			CoEdge edge = (CoEdge)val;
			circularEdgeChecker.setSelected(edge.info.circularHoleEdge);
		}
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (selectedNodesList.getSelectedValue() == null) return;
		for (Object s : selectedNodesList.getSelectedValues()) {
			if (!(s instanceof CoVertex)) continue;
			CoVertex v = (CoVertex)s;
			double thetaDeg = customThetaModel.getNumber().doubleValue();
			v.info.theta = Math.toRadians(thetaDeg);
		}
	}
	
	
	@Override
	public void propertyChange(PropertyChangeEvent evt) {
		Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
		if (SwingWorker.StateValue.DONE == evt.getNewValue()) {
			Unwrap unwrapper = (Unwrap)evt.getSource();
			if (unwrapper.isCancelled()) {
				System.out.println("Unwrap jop cancelled: " + unwrapper.getState());
				return;
			}
			try {
				surface = unwrapper.get();
			} catch (InterruptedException e) {
				return;
			} catch (ExecutionException e) {
				String name = e.getCause().getClass().getSimpleName();
				String msg = e.getCause().getLocalizedMessage();
				StringBuffer stackBuffer = new StringBuffer(name + ": " + msg + "\n");
				for (StackTraceElement ste : e.getCause().getStackTrace()) {
					stackBuffer.append(ste.toString() + "\n");
				}
				int result = JOptionPane.showConfirmDialog(w, msg + "\nShow stack trace?", name, OK_CANCEL_OPTION, ERROR_MESSAGE);
				if (result == JOptionPane.OK_OPTION) {
					JOptionPane.showMessageDialog(w, stackBuffer.toString(), name, ERROR_MESSAGE);
				}
				return;
			} finally {
				unwrapper.getSurface().revertNormalization();
			}
			genus = unwrapper.genus;
			cutInfo = unwrapper.cutInfo;
			metricErrorAdapter.setLengthMap(unwrapper.lengthMap);
			metricErrorAdapter.setSignature(Pn.EUCLIDEAN);

			createVisualization(surface, genus, cutInfo);
			updateSurface();
			updateStates();
		}
	}
	
	
	public void createVisualization(CoHDS surface, int genus, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		this.genus = genus;
		this.surface = surface;
		if (genus > 0) {
			cutColorAdapter.setContext(cutInfo);
			cutRadiusAdapter.setContext(cutInfo);
			pointRadiusAdapter.setContext(cutInfo);
			pointColorAdapter.setContext(cutInfo);
		}
		if (genus > 1) {
			try {
			System.out.println("Constructing fundamental cut polygon...");
			cuttedPolygon = FundamentalPolygonUtility.constructFundamentalPolygon(cutInfo);
			System.out.println(cuttedPolygon);
			FundamentalVertex root = cuttedPolygon.getMaxValenceVertex();
			System.out.println("Constructing minimal polygon...");
			minimalPolygon = FundamentalPolygonUtility.minimize(cuttedPolygon, root);
			System.out.println(minimalPolygon);
			minimalPolygon.checkRelation();
			System.out.println("Constructing fast canonical polygon...");
			canonicalPolygon = FundamentalPolygonUtility.canonicalize(minimalPolygon, useDistanceToCanonicalize.isSelected());
			System.out.println(canonicalPolygon);
			canonicalPolygon.checkRelation();
			updatePolygonTexture(coverResolution);
			metricErrorAdapter.setSignature(Pn.HYPERBOLIC);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		if (genus == 0) {
			kleinButton.setSelected(true);
			cutInfo = null;
			cuttedPolygon = null;
		}
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object s = e.getSource();
		if (domainCombo == s) {
			updateStates();
			return;
		}
		if (coverToTextureButton == s) {
			AppearanceInspector ai = contentAppearance.getAppearanceInspector();
			TextureInspector ti = ai.getTextureInspector();
			if (cutCoverImage != null) {
				ti.addTexture("Cut Polygon", cutCoverImage);
			}
			if (minimalCoverImage != null) {
				ti.addTexture("Minimal Polygon", minimalCoverImage);
			}
			if (canonicalCoverImage != null) {
				ti.addTexture("Canonical Polygon", canonicalCoverImage);
			}
		}
		if (unwrapBtn == s) {
			CoHDS surface = getLoaderGeometry();
			if (surface == null) return;
			surface.normalizeCoordinates();
			AdapterSet aSet = hif.getAdapters();
			Unwrap uw = new Unwrap(surface, aSet);
			uw.setToleranceExponent(toleranceExpModel.getNumber().intValue());
			uw.setMaxIterations(maxIterationsModel.getNumber().intValue());
			uw.setNumCones(numConesModel.getNumber().intValue());
			uw.setQuantizationMode((QuantizationMode)quantizationModeCombo.getSelectedItem());
			uw.setBoundaryQuantMode((QuantizationMode)boundaryQuantizationCombo.getSelectedItem());
			uw.setBoundaryMode((BoundaryMode)boundaryModeCombo.getSelectedItem());
			uw.setUsePetsc(numericsCombo.getSelectedIndex() == 1);
			uw.setSelectedVertices(hif.getSelection().getVertices(surface));
			uw.addPropertyChangeListener(this);
			uw.execute();
		}
		if (showUnwrapped == s || useProjectiveTexture == s) {
			updateSurface();
			updateStates();
		}
		if (kleinButton == s || poincareButton == s || halfplaneButton == s) {
			updateSurface();
			if (genus > 1) {
				updatePolygonTexture(coverResolution);
			}
			updateStates();
		}
		if (customModeCombo == s) {
			for (Object sel : selectedNodesList.getSelectedValues()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.boundaryMode = (BoundaryMode)customModeCombo.getSelectedItem();
			}
		}
		if (customQuantizationCombo == s) {
			for (Object sel : selectedNodesList.getSelectedValues()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.quantizationMode = (QuantizationMode)customQuantizationCombo.getSelectedItem();
			}
		}
		if (useCustomThetaChecker == s) {
			for (Object sel : selectedNodesList.getSelectedValues()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.useCustomTheta = useCustomThetaChecker.isSelected(); 
			}
		}
		if (circularEdgeChecker == s) {
			for (Object sel : selectedNodesList.getSelectedValues()) {
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
				Theta<CoVertex> theta = new CTheta();
				Variable<CoVertex, CoEdge> variable = new CVariable();
				Lambda<CoEdge> lambda = new CLambda();
				Alpha<CoEdge> alpha = new CAlpha();
				InitialEnergy<CoFace> initE = new CInitialEnergy();
				EuclideanFunctional<CoVertex, CoEdge, CoFace> fun = new EuclideanFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, initE);
				prepareInvariantDataEuclidean(fun, hds, boundaryMode, boundaryQuantMode, hif.getAdapters());
			} catch (Exception e1) {
				Window w = getWindowAncestor(shrinkPanel);
				showMessageDialog(w, e1.getMessage(), "Error", ERROR_MESSAGE);
			}
		}
		if (quantizeToQuads == s) {
			AdapterSet a = hif.getAdapters();
			HalfedgeSelection sel = hif.getSelection();
			List<CoVertex> cones = new LinkedList<CoVertex>();
			CoHDS hds = hif.get(new CoHDS());
			for (CoVertex v : hds.getVertices()) {
				if (v.info != null && v.info.useCustomTheta) {
					cones.add(v);
				}
			}
			if (cones.size() > 2) throw new RuntimeException("More than two cones not supported");
			Matrix T = contentAppearance.getAppearanceInspector().getTextureMatrix();
			Matrix Ti = new Matrix(T); Ti.invert();
			if (cones.size() == 0) return;
			if (cones.size() == 1) { // translation only
				CoVertex cv = cones.get(0);
				double[] texpos = a.getD(TexturePosition4d.class, cv);
				double[] quantPos = texpos.clone();
				T.transformVector(quantPos);
				Pn.dehomogenize(quantPos, quantPos);
				double offset = sel.isSelected(cv) ? 0.25 : 0;
				double difX = (quantPos[0] + offset) % 0.5;
				double difY = (quantPos[1] + offset) % 0.5;
				double[] difVec = {difX, difY, 0, 0};
				Ti.transformVector(difVec);
				Matrix QT = MatrixBuilder.euclidean().translate(-difVec[0], -difVec[1], 0).getMatrix();
				for (CoVertex v : hds.getVertices()) {
					double[] tp = hif.getAdapters().getD(TexturePosition4d.class, v);
					double[] qtp = tp.clone();
					QT.transformVector(qtp);
					a.set(TexturePosition.class, v, qtp);
				}
			}
			if (cones.size() == 2) { // affine transform
				// TODO fix this!
				CoVertex cv1 = cones.get(0);
				CoVertex cv2 = cones.get(1);
				double[] texpos1 = a.getD(TexturePosition4d.class, cv1);
				double[] texpos2 = a.getD(TexturePosition4d.class, cv2);
				double[] quantPos1 = texpos1.clone();
				double[] quantPos2 = texpos2.clone();
				T.transformVector(quantPos1);
				T.transformVector(quantPos2);
				Pn.dehomogenize(quantPos1, quantPos1);
				Pn.dehomogenize(quantPos2, quantPos2);
				double dist = Rn.euclideanDistance(quantPos1, quantPos2);
				double sDist = dist;
				double distOffset = (!sel.isSelected(cv2) && sel.isSelected(cv1)) || (sel.isSelected(cv2) && !sel.isSelected(cv1)) ? 0.25 : 0;
				dist -= (dist + distOffset) % 0.5;
				dist = Math.max(dist, 0.5);
				double angle = atan2(quantPos1[1] - quantPos2[1], quantPos1[0] - quantPos2[0]) % PI/2;
				double offset = sel.isSelected(cv1) ? 0.25 : 0;
				double difX1 = (quantPos1[0] + offset) % 0.5;
				double difY1 = (quantPos1[1] + offset) % 0.5;
				
				MatrixBuilder PivotB = MatrixBuilder.euclidean();
				PivotB.translate(quantPos1[0], quantPos1[1], 0);
				
				MatrixBuilder QB = MatrixBuilder.euclidean();
				QB.scale(dist / sDist);
				QB.rotate(-angle, new double[]{0,0,1});
				QB.conjugateBy(PivotB.getMatrix().getArray());
				QB.translate(-difX1, -difY1, 0);				
				Matrix QT = QB.getMatrix();
				QT.transformVector(quantPos1);
				QT.transformVector(quantPos2);
				System.out.println("new dist: " + Rn.euclideanDistance(quantPos1, quantPos2));
				System.out.println("new angle: " + atan2(quantPos1[1] - quantPos2[1], quantPos1[0] - quantPos2[0]) % 2*PI);
				QT = Matrix.conjugate(QT, Ti);
				for (CoVertex v : hds.getVertices()) {
					double[] tp = hif.getAdapters().getD(TexturePosition4d.class, v);
					double[] qtp = tp.clone();
					QT.transformVector(qtp);
					Pn.dehomogenize(qtp, qtp);
					a.set(TexturePosition.class, v, qtp);
				}
			}
			System.out.println("Quantizing texture cones: " + cones);
			hif.update();
		}
		if (moveToCenterButton == s) {
			Set<CoVertex> sel = hif.getSelection().getVertices(surface);
			if (sel.isEmpty()) return;
			CoVertex v = sel.iterator().next();
			double[] pos = v.T;
			MatrixBuilder mb = MatrixBuilder.hyperbolic();
			mb.translateFromTo(pos, new double[] {0,0,0,1});
			Matrix T = mb.getMatrix();
			for (CoVertex vv : surface.getVertices()) {
				T.transformVector(vv.T);
			}
			hif.update();
			createVisualization(surface, genus, cutInfo);
			updateSurface();
			updateStates();
		}
		if (resetXFormButton == s) {
			
		}
		if (saveTextureButton == s) {
			JFileChooser imageChooser = new JFileChooser();
			imageChooser.addChoosableFileFilter(new FileFilter() {
				
				@Override
				public String getDescription() {
					return "PNG Files (*.png)";
				}
				
				@Override
				public boolean accept(File f) {
					return f.isDirectory() || f.getName().toLowerCase().endsWith(".png");
				}
			});
			Window w = SwingUtilities.getWindowAncestor(this.shrinkPanel);
			int result = imageChooser.showSaveDialog(w);
			if (result != JFileChooser.APPROVE_OPTION) {
				return;
			}
			File file = imageChooser.getSelectedFile();
			if (!file.getName().toLowerCase().endsWith(".png")) {
				file = new File(file.getAbsolutePath() + ".png");
			}
			try {
				switch ((Domain)domainCombo.getSelectedItem()) {
				case Cut:
					ImageIO.write(cutCoverImage, "png", imageChooser.getSelectedFile());
					break;
				case Canonical:
					ImageIO.write(canonicalCoverImage, "png", imageChooser.getSelectedFile());
					break;
				case Minimal:
					ImageIO.write(minimalCoverImage, "png", imageChooser.getSelectedFile());
					break;
				}
			} catch (Exception e2) {
				JOptionPane.showMessageDialog(w, e2.getMessage(), "Error", ERROR_MESSAGE);
			}
		}
	}
	
	/*
	 * Set up scene
	 */
	private void updateStates() {
		HalfedgeLayer l = hif.getActiveLayer();
		l.removeTemporaryGeometry(cutCoverRoot);
		l.removeTemporaryGeometry(minimalCoverRoot);
		l.removeTemporaryGeometry(canonicalCoverRoot);
		if (genus > 1 && showUnwrapped.isSelected()) {
			switch ((Domain)domainCombo.getSelectedItem()) {
			case Cut:
				l.addTemporaryGeometry(cutCoverRoot);
				break;
			case Minimal:
				l.addTemporaryGeometry(minimalCoverRoot);
				break;
			case Canonical:
				l.addTemporaryGeometry(canonicalCoverRoot);
				break;
			}
			boolean showCircle = false;
			showCircle |= getSelectedModel() == HyperbolicModel.Poincaré;
			showCircle |= getSelectedModel() == HyperbolicModel.Klein;
			unitCircle.setVisible(showCircle);
		}
		ImageData imgData = null;
		switch ((Domain)domainCombo.getSelectedItem()) {
		case Cut:
			if (cutCoverImage != null) {
				imgData = new ImageData(cutCoverImage);
			}
			break;
		case Minimal:
			if (minimalCoverImage != null) {
				imgData = new ImageData(minimalCoverImage);
			}
			break;
		case Canonical:
			if (canonicalCoverImage != null) {
				imgData = new ImageData(canonicalCoverImage);
			}
			break;
		}
		if (imgData != null) {
			Texture2D tex2d = TextureUtility.createTexture(universalCoverAppearance, POLYGON_SHADER, imgData);
			tex2d.setTextureMatrix(polygonTextureMatrix);
		} else {
			TextureUtility.removeTexture(universalCoverAppearance, POLYGON_SHADER);
		}
	}
	
	
	private CoHDS getLoaderGeometry() {
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
		Triangulator.triangulate(surface);
		return surface;
	}
	
	
	private void updateSurface() {
		if (surface == null) {
			return;
		}
		hif.clearSelection();
		texturePositionAdapter.setProjective(useProjectiveTexture.isSelected());
		texturePositionAdapter.setModel(getSelectedModel());
		texCoordPositionAdapter.setModel(getSelectedModel());
		hif.addAdapter(metricErrorAdapter, false);
		if (showUnwrapped.isSelected()) {
			hif.addAdapter(texCoordPositionAdapter, false);	
		}
		hif.set(surface);
	}
	


	public void updatePolygonTexture(int resolution) {
		int coverRecusionDepth = coverRecursionModel.getNumber().intValue();
		HyperbolicModel model = getSelectedModel();
		cutCoverImage = VisualizationUtility.drawUniversalCoverImage(
			cuttedPolygon, 
			coverRecusionDepth, 
			model, 
			resolution,
			Color.BLUE
		);
		minimalCoverImage = VisualizationUtility.drawUniversalCoverImage(
			minimalPolygon, 
			coverRecusionDepth,
			model, 
			resolution,
			Color.GREEN
		);
		canonicalCoverImage = VisualizationUtility.drawUniversalCoverImage(
			canonicalPolygon, 
			coverRecusionDepth, 
			model, 
			resolution,
			Color.RED
		);
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		hif.addAdapter(new CoPositionAdapter(), true);
		hif.addAdapter(texturePositionAdapter, true);
		hif.addSelectionListener(this);
		contentAppearance = c.getPlugin(ContentAppearance.class);
	}
	
	@Override
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "numCones", numConesModel.getNumber().intValue());
		c.storeProperty(getClass(), "quantizeMode", quantizationModeCombo.getSelectedItem());
		c.storeProperty(getClass(), "boundaryMode", boundaryModeCombo.getSelectedItem());
		c.storeProperty(getClass(), "boundaryQuantMode", boundaryQuantizationCombo.getSelectedItem());
		c.storeProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex());
		c.storeProperty(getClass(), "klein", kleinButton.isSelected()); 
		c.storeProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected());
		c.storeProperty(getClass(), "coverRecursion", coverRecursionModel.getNumber());
		c.storeProperty(getClass(), "quantizationMode", quantizationModeCombo.getSelectedIndex());
		c.storeProperty(getClass(), "useProjectiveTexture", useProjectiveTexture.isSelected());
		c.storeProperty(getClass(), "toleranceExponent", toleranceExpModel.getNumber());
		c.storeProperty(getClass(), "maxIterations", maxIterationsModel.getNumber());
		c.storeProperty(getClass(), "boundaryPanelShrinked", boundaryPanel.isShrinked());	
		c.storeProperty(getClass(), "conesPanelShrinked", coneConfigPanel.isShrinked());	
		c.storeProperty(getClass(), "visualizationPanelShrinked", visualizationPanel.isShrinked());	
		c.storeProperty(getClass(), "modelPanelShrinked", modelPanel.isShrinked());	
		c.storeProperty(getClass(), "coneTexQuantPanelShrinked", texQuantizationPanel.isShrinked());
		c.storeProperty(getClass(), "customVertexPanelShrinked", customNodePanel.isShrinked());
	} 
	
 
	@Override
	public void restoreStates(Controller c) throws Exception { 
		super.restoreStates(c);
		numConesModel.setValue(c.getProperty(getClass(), "numCones", numConesModel.getNumber().intValue()));
		numericsCombo.setSelectedIndex(c.getProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex()));
		kleinButton.setSelected(c.getProperty(getClass(), "klein", kleinButton.isSelected()));
		poincareButton.setSelected(!kleinButton.isSelected());
		showUnwrapped.setSelected(c.getProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected()));
		coverRecursionModel.setValue(c.getProperty(getClass(), "coverRecursion", coverRecursionModel.getNumber()));
		quantizationModeCombo.setSelectedIndex(c.getProperty(getClass(), "quantizationMode", quantizationModeCombo.getSelectedIndex()));
		useProjectiveTexture.setSelected(c.getProperty(getClass(), "useProjectiveTexture", useProjectiveTexture.isSelected()));
		toleranceExpModel.setValue(c.getProperty(getClass(), "toleranceExponent", toleranceExpModel.getNumber()));
		maxIterationsModel.setValue(c.getProperty(getClass(), "maxIterations", maxIterationsModel.getNumber()));
		quantizationModeCombo.setSelectedItem(c.getProperty(getClass(), "quantizeMode", quantizationModeCombo.getSelectedItem()));
		boundaryModeCombo.setSelectedItem(c.getProperty(getClass(), "boundaryMode", boundaryModeCombo.getSelectedItem()));
		boundaryQuantizationCombo.setSelectedItem(c.getProperty(getClass(), "boundaryQuantMode", boundaryQuantizationCombo.getSelectedItem()));
		boundaryPanel.setShrinked(c.getProperty(getClass(), "boundaryPanelShrinked", true));
		coneConfigPanel.setShrinked(c.getProperty(getClass(), "conesPanelShrinked", true));
		visualizationPanel.setShrinked(c.getProperty(getClass(), "visualizationPanelShrinked", true));
		modelPanel.setShrinked(c.getProperty(getClass(), "modelPanelShrinked", true));
		texQuantizationPanel.setShrinked(c.getProperty(getClass(), "coneTexQuantPanelShrinked", true));
		customNodePanel.setShrinked(c.getProperty(getClass(), "customVertexPanelShrinked", customNodePanel.isShrinked()));
	}
	
	
	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = new PluginInfo();
		info.name = "Discrete Conformal Parametrization";
		info.vendorName = "Stefan Sechelmann";
		info.email = "sechel@math.tu-berlin.de";
		return info;
	}

	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}
	

	protected HyperbolicModel getSelectedModel() {
		if (kleinButton.isSelected()) {
			return HyperbolicModel.Klein;
		}
		if (poincareButton.isSelected()) {
			return HyperbolicModel.Poincaré;
		}
		if (halfplaneButton.isSelected()) {
			return HyperbolicModel.Halfplane;
		}
		return HyperbolicModel.Klein;
	}
	
}
