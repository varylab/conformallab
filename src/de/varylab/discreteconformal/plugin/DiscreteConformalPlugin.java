package de.varylab.discreteconformal.plugin;

import static de.jreality.math.MatrixBuilder.euclidean;
import static de.jreality.scene.Appearance.DEFAULT;
import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.LIGHTING_ENABLED;
import static de.jreality.shader.CommonAttributes.LINE_SHADER;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.TEXTURE_2D;
import static de.jreality.shader.CommonAttributes.TUBE_RADIUS;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.varylab.discreteconformal.uniformization.FundamentalDomainUtility.createFundamentalPolygon;
import static de.varylab.discreteconformal.uniformization.UniformizationUtility.constructFundamentalPolygon;
import static java.awt.Color.BLACK;
import static java.awt.Color.RED;
import static java.awt.Color.WHITE;
import static java.lang.Math.PI;
import static javax.swing.JOptionPane.ERROR_MESSAGE;
import static javax.swing.JOptionPane.OK_CANCEL_OPTION;
import static javax.swing.JOptionPane.showMessageDialog;
import static javax.swing.SwingUtilities.getWindowAncestor;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;

import javax.swing.ButtonGroup;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
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
import de.jreality.math.Pn;
import de.jreality.plugin.basic.View;
import de.jreality.scene.Appearance;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.shader.ImageData;
import de.jreality.shader.Texture2D;
import de.jreality.shader.TextureUtility;
import de.jtem.halfedge.Vertex;
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
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
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
import de.varylab.discreteconformal.uniformization.FundamentalDomainUtility;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.util.UnwrapUtility.QuantizationMode;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin implements ListSelectionListener, ChangeListener, ActionListener, PropertyChangeListener, SelectionListener {

	private static int
		polyResolution = 1000,
		coverRecursion = 3,
		coverResolution = 2048;
	
	// plug-in section ------------------ 
	private HalfedgeInterface
		hif = null;
	
	// data section ---------------------
	private CoHDS
		surface = null;
	private CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	private List<CoVertex>
		customVertices = new LinkedList<CoVertex>();
	private FundamentalPolygon 
		fundamentalPolygon = null;
	private Matrix 
		polygonTextureMatrix = euclidean().translate(-0.5, -0.5, 0).scale(0.5).scale(1, -1, 1).getMatrix();
	private Image
		polygonImage = null;
	private int
		genus = -1;

	private MetricErrorAdapter
		edgeLengthAdapter = new MetricErrorAdapter();
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
		fundamentalPolygonRoot = new SceneGraphComponent("Fundamental Polygon"), 
		unitCircle = new SceneGraphComponent("Hyperbolic Boundary"),
		universalCoverRoot = new SceneGraphComponent("Universal Cover");

	// user interface section ------------
	private JButton
		checkGaussBonnetBtn = new JButton("Check Gauß-Bonnet"),
		unwrapBtn = new JButton("Unwrap"),
		normalizePolygonBtn = new JButton("Normalize Polygon"),
		normalizePolygonFastBtn = new JButton("Normalize Polygon Fast");
	private ShrinkPanel
		customVertexPanel = new ShrinkPanel("Custom Vertices"),
		boundaryPanel = new ShrinkPanel("Boundary"),
		coneConfigPanel = new ShrinkPanel("Automatic Cones"),
		modelPanel = new ShrinkPanel("Hyperbolic Model"),
		visualizationPanel = new ShrinkPanel("Visualization");
	private SpinnerNumberModel
		customThetaModel = new SpinnerNumberModel(360.0, 0.0, 1000.0, 1.0),
		numConesModel = new SpinnerNumberModel(0, 0, 100, 1),
		toleranceExpModel = new SpinnerNumberModel(-8, -30, -1, 1),
		maxIterationsModel = new SpinnerNumberModel(150, 1, 10000, 1);
	private JSpinner
		customThetaSpinner = new JSpinner(customThetaModel),
		numConesSpinner = new JSpinner(numConesModel),
		toleranceExpSpinner = new JSpinner(toleranceExpModel),
		maxIterationsSpinner = new JSpinner(maxIterationsModel);
	private JCheckBox
		useCustomThetaChecker = new JCheckBox("Custom Theta"),
		showFundamentalPolygon = new JCheckBox("Fundamental Polygon"),
		useProjectiveTexture = new JCheckBox("Projective Texture", true),
		showUnwrapped = new JCheckBox("Show Unwrapped"),
		showUniversalCover = new JCheckBox("Universal Cover");
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
		selectedVertexList = new JList();
	private JScrollPane
		selectionScroller = new JScrollPane(selectedVertexList);
		
	public DiscreteConformalPlugin() {
		createLayout();
		unwrapBtn.addActionListener(this);
		checkGaussBonnetBtn.addActionListener(this);
		showUnwrapped.addActionListener(this);
		kleinButton.addActionListener(this);
		poincareButton.addActionListener(this);
		halfplaneButton.addActionListener(this);
		showUniversalCover.addActionListener(this);
		showFundamentalPolygon.addActionListener(this);
		useProjectiveTexture.addActionListener(this);
		normalizePolygonBtn.addActionListener(this);
		normalizePolygonFastBtn.addActionListener(this);
		
		ButtonGroup modelGroup = new ButtonGroup();
		modelGroup.add(kleinButton);
		modelGroup.add(poincareButton);
		modelGroup.add(halfplaneButton);
		
		Appearance funPolyApp = new Appearance();
		funPolyApp.setAttribute(VERTEX_DRAW, false);
		funPolyApp.setAttribute(EDGE_DRAW, true);
		funPolyApp.setAttribute(LINE_SHADER + "." + TUBE_RADIUS, 0.03);
		funPolyApp.setAttribute(LINE_SHADER + "." + DIFFUSE_COLOR, RED);
		fundamentalPolygonRoot.setAppearance(funPolyApp);
		
		IndexedFaceSetFactory ifsf = new IndexedFaceSetFactory();
		ifsf.setVertexCount(4);
		ifsf.setFaceCount(1);
		ifsf.setFaceIndices(new int[][] {{0,1,2,3}});
		ifsf.setVertexTextureCoordinates(new double[] {-1,-1,1,-1,1,1,-1,1});
		ifsf.setVertexCoordinates(new double[]{-1,-1,-0.001, 1,-1,-0.001,  1,1,-0.001,  -1,1,-0.001});
		ifsf.setGenerateFaceNormals(true);
		ifsf.update();
		universalCoverRoot.setGeometry(ifsf.getGeometry());
		universalCoverRoot.setAppearance(universalCoverAppearance);
		universalCoverAppearance.setAttribute(VERTEX_DRAW, false);
		universalCoverAppearance.setAttribute(EDGE_DRAW, false);
		universalCoverAppearance.setAttribute(FACE_DRAW, true);
		universalCoverAppearance.setAttribute(LIGHTING_ENABLED, false);
		universalCoverAppearance.setAttribute(POLYGON_SHADER + "." + DIFFUSE_COLOR, WHITE);
		
		Appearance circleApp = new Appearance();
		circleApp.setAttribute(EDGE_DRAW, false);
		circleApp.setAttribute(VERTEX_DRAW, false); 
		circleApp.setAttribute(FACE_DRAW, true);
		circleApp.setAttribute(POLYGON_SHADER + "." + DIFFUSE_COLOR, BLACK);
		circleApp.setAttribute(POLYGON_SHADER + "." + TEXTURE_2D, DEFAULT); 
		unitCircle.setAppearance(circleApp);
		euclidean().rotate(PI / 2, 1, 0, 0).assignTo(unitCircle);
		unitCircle.setGeometry(Primitives.torus(1.0025, 0.005, 200, 5));
		universalCoverRoot.addChild(unitCircle);
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
		
		customVertexPanel.add(selectionScroller, c2);
		selectionScroller.setPreferredSize(new Dimension(10, 70));
		selectionScroller.setMinimumSize(new Dimension(10, 70));
		customVertexPanel.add(useCustomThetaChecker, c1);
		customVertexPanel.add(customThetaSpinner, c2);
		customVertexPanel.add(new JLabel("Mode"), c1);
		customVertexPanel.add(customModeCombo, c2);
		customVertexPanel.add(new JLabel("Quantization"), c1);
		customVertexPanel.add(customQuantizationCombo, c2);
		selectedVertexList.getSelectionModel().addListSelectionListener(this);
		customModeCombo.addActionListener(this);
		customQuantizationCombo.addActionListener(this);
		customVertexPanel.setShrinked(true);
		useCustomThetaChecker.addActionListener(this);
		customThetaSpinner.addChangeListener(this);
		shrinkPanel.add(customVertexPanel, c2);
		
		visualizationPanel.setLayout(new GridBagLayout());
		visualizationPanel.add(showUnwrapped, c2);
		visualizationPanel.add(showFundamentalPolygon, c2);
		visualizationPanel.add(normalizePolygonBtn, c2);
		visualizationPanel.add(normalizePolygonFastBtn, c2);
		visualizationPanel.add(showUniversalCover, c2);
		visualizationPanel.add(useProjectiveTexture, c2);
		visualizationPanel.setShrinked(true);
		shrinkPanel.add(visualizationPanel, c2);
		
		modelPanel.setLayout(new GridLayout(1, 3, 2, 2));
		modelPanel.add(kleinButton);
		modelPanel.add(poincareButton);
		modelPanel.add(halfplaneButton);
		modelPanel.setShrinked(true);
		shrinkPanel.add(modelPanel, c2);
		
		normalizePolygonBtn.setEnabled(false);
		normalizePolygonFastBtn.setEnabled(false);
	}
	
	
	@Override
	public void selectionChanged(HalfedgeSelection s, HalfedgeInterface hif) {
		Set<CoVertex> oldSelection = new HashSet<CoVertex>(customVertices);
		customVertices.clear();
		for (Vertex<?,?,?> v : s.getVertices()) {
			if (v instanceof CoVertex) {
				CoVertex cov = (CoVertex)v;
				if (cov.info == null) {
					cov.info = new CustomVertexInfo();
				}
				customVertices.add(cov);
			}
		}
		Set<CoVertex> newSelection = new HashSet<CoVertex>(customVertices);
		newSelection.removeAll(oldSelection);
		DefaultListModel model = new DefaultListModel();
		for (CoVertex v : customVertices) {
			model.addElement(v);
		}
		selectedVertexList.setModel(model);
		if (!newSelection.isEmpty()) {
			CoVertex v = newSelection.iterator().next();
			selectedVertexList.setSelectedValue(v, true);
		}
	}

	
	@Override
	public void valueChanged(ListSelectionEvent e) {
		if (selectedVertexList.getSelectedValue() == null) return;
		CoVertex v = (CoVertex)selectedVertexList.getSelectedValue();
		customModeCombo.setSelectedItem(v.info.boundaryMode);
		customQuantizationCombo.setSelectedItem(v.info.quantizationMode);
		useCustomThetaChecker.setSelected(v.info.useCustomTheta);
		customThetaModel.setValue(Math.toDegrees(v.info.theta));
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (selectedVertexList.getSelectedValue() == null) return;
		for (Object s : selectedVertexList.getSelectedValues()) {
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
			edgeLengthAdapter.setLengthMap(unwrapper.lengthMap);
			edgeLengthAdapter.setSignature(Pn.EUCLIDEAN);
			if (genus > 0) {
				cutInfo = unwrapper.cutInfo;
				cutColorAdapter.setContext(cutInfo);
				cutRadiusAdapter.setContext(cutInfo);
				pointRadiusAdapter.setContext(cutInfo);
				pointColorAdapter.setContext(cutInfo);
			}
			if (genus > 1) {
				fundamentalPolygon = constructFundamentalPolygon(cutInfo);
				updateFundamentalPolygon(polyResolution);
				updatePolygonTexture(coverRecursion, coverResolution);
				normalizePolygonBtn.setEnabled(true);
				normalizePolygonFastBtn.setEnabled(true);
				edgeLengthAdapter.setSignature(Pn.HYPERBOLIC);
			} else {
				normalizePolygonBtn.setEnabled(false);
				normalizePolygonFastBtn.setEnabled(false);
			}
			if (genus == 0) {
				kleinButton.setSelected(true);
				cutInfo = null;
				fundamentalPolygon = null;
			}
			updateSurface();
			updateStates();
		}
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object s = e.getSource();
		if (showFundamentalPolygon == s ||
			showUniversalCover == s
		) {
			updateStates();
			return;
		}
		if (unwrapBtn == s) {
			CoHDS surface = getLoaderGeometry();
			surface.normalizeCoordinates();
			Unwrap uw = new Unwrap(surface, hif.getAdapters());
			uw.setToleranceExponent(toleranceExpModel.getNumber().intValue());
			uw.setMaxIterations(maxIterationsModel.getNumber().intValue());
			uw.setNumCones(numConesModel.getNumber().intValue());
			uw.setQuantizationMode((QuantizationMode)quantizationModeCombo.getSelectedItem());
			uw.setBoundaryQuantMode((QuantizationMode)boundaryQuantizationCombo.getSelectedItem());
			uw.setBoundaryMode((BoundaryMode)boundaryModeCombo.getSelectedItem());
			uw.setUsePetsc(numericsCombo.getSelectedIndex() == 1);
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
				updateFundamentalPolygon(polyResolution);
				updatePolygonTexture(coverRecursion, coverResolution);
			}
			updateStates();
		}
		if (customModeCombo == s) {
			for (Object sel : selectedVertexList.getSelectedValues()) {
				CoVertex v = (CoVertex)sel;
				v.info.boundaryMode = (BoundaryMode)customModeCombo.getSelectedItem();
			}
		}
		if (customQuantizationCombo == s) {
			for (Object sel : selectedVertexList.getSelectedValues()) {
				CoVertex v = (CoVertex)sel;
				v.info.quantizationMode = (QuantizationMode)customQuantizationCombo.getSelectedItem();
			}
		}
		if (useCustomThetaChecker == s) {
			for (Object sel : selectedVertexList.getSelectedValues()) {
				CoVertex v = (CoVertex)sel;
				v.info.useCustomTheta = useCustomThetaChecker.isSelected(); 
			}
		}
		if (checkGaussBonnetBtn == s) {
			CoHDS hds = hif.get(new CoHDS());
			BoundaryMode boundaryMode = (BoundaryMode)boundaryModeCombo.getSelectedItem();
			QuantizationMode boundaryQuantMode = (QuantizationMode)boundaryQuantizationCombo.getSelectedItem();
			try {
				UnwrapUtility.prepareInvariantDataEuclidean(
					hds, 
					boundaryMode, 
					boundaryQuantMode, 
					hif.getAdapters()
				);
			} catch (Exception e1) {
				Window w = getWindowAncestor(shrinkPanel);
				showMessageDialog(w, e1.getMessage(), "Error", ERROR_MESSAGE);
			}
		}
		if (normalizePolygonBtn == s) {
			fundamentalPolygon = fundamentalPolygon.getNaiveCanonical();
			updateFundamentalPolygon(polyResolution);
			updatePolygonTexture(coverRecursion, coverResolution);
		}
		if (normalizePolygonFastBtn == s) {
			fundamentalPolygon = fundamentalPolygon.getFastCanonical();
			updateFundamentalPolygon(polyResolution);
			updatePolygonTexture(coverRecursion, coverResolution);
		}
	}
	
	/*
	 * Set up scene
	 */
	private void updateStates() {
		HalfedgeLayer l = hif.getActiveLayer();
		l.removeTemporaryGeometry(fundamentalPolygonRoot);
		l.removeTemporaryGeometry(universalCoverRoot);
		if (genus > 1 && showUnwrapped.isSelected()) {
			if (showFundamentalPolygon.isSelected()) {
				l.addTemporaryGeometry(fundamentalPolygonRoot);
			} 
			if (showUniversalCover.isSelected()) {
				l.addTemporaryGeometry(universalCoverRoot);
				unitCircle.setVisible(getSelectedModel() == HyperbolicModel.Poincaré);
			}
		}
		ImageData imgData = null;
		if (polygonImage != null && (showUniversalCover.isSelected())) {
			imgData = new ImageData(polygonImage);
		}
		if (polygonImage != null && showUniversalCover.isSelected()) {
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
		hif.addGlobalAdapter(edgeLengthAdapter, false);
		if (showUnwrapped.isSelected()) {
			hif.addGlobalAdapter(texCoordPositionAdapter, false);	
		}
//		if (genus >= 1) {
//			hif.addLayerAdapter(cutRadiusAdapter,false);
//			hif.addLayerAdapter(cutColorAdapter,false);
//			hif.addLayerAdapter(pointRadiusAdapter,false);
//			hif.addLayerAdapter(pointColorAdapter,false);
//		} else {
//			hif.removeAdapter(cutRadiusAdapter);
//			hif.removeAdapter(cutColorAdapter);
//			hif.removeAdapter(pointRadiusAdapter);
//			hif.removeAdapter(pointColorAdapter);
//		}
		hif.set(surface);
	}
	

	public void updateFundamentalPolygon(int resolution) {
		double[] rootPos = cutInfo.cutRoot.T;
		createFundamentalPolygon(
			fundamentalPolygon, 
			rootPos, 
			fundamentalPolygonRoot, 
			resolution, 
			getSelectedModel()
		);
	}
	
	
	
	public void updatePolygonTexture(int depth, int resolution) {
		double[] pRoot = cutInfo.cutRoot.T;
		double[] root = new double[] {pRoot[0], pRoot[1], 0.0, pRoot[3]};
		HyperbolicModel model = getSelectedModel();
		polygonImage = FundamentalDomainUtility.createCoverTexture(root, fundamentalPolygon, depth, model, resolution);
		updateStates();
	}
	
	
	public void copyAtEdge(int edgeIndex) {
//		if (surface == null) {
//			return;
//		}
//		final boolean klein = kleinButton.isSelected(); 
//		CoEdge edge = null;
//		int i = 0;
//		for (CoEdge e : surface.getEdges()) {
//			if (e.isPositive()) {
//				if (i == edgeIndex) {
//					if (e.getLeftFace() != null) {
//						edge = e;
//					} else {
//						edge = e.getOppositeEdge();
//					}
//					break;
//				}
//				i++;
//			}
//		}
//		
//		if (edge == null) {
//			System.err.println("Edge corresponding to index " + edgeIndex + " not found");
//			return;
//		}
//		CoEdge coEdge = cutInfo.edgeCutMap.get(edge);
//		if (coEdge == null) {
//			System.err.println("CoEdge not found");
//			return;
//		}
//		if (edge.getOppositeEdge().getLeftFace() != null) {
//			System.err.println("Picked no boundary edge!");
//			return;
//		}
//		System.out.println("hyperbolic motion: " + edge + " -> " + coEdge);
		
//		final Matrix A = UniformizationUtility.makeHyperbolicMotion(coEdge, edge.getOppositeEdge());

//		ConverterHeds2JR converter = new ConverterHeds2JR();
//		Adapter texAdapter = new TextCoordsAdapter2Ifs<CoVertex> () {
//
//			@Override
//			public double[] getTextCoordinate(CoVertex v) {
//				Point t = v.getTextureCoord();
//				double[] raw = new double[] {t.x(), t.y(), 0.0, t.z()};
//				A.transformVector(raw);
//				if (klein) {
//					return new double[] {raw[0], raw[1], 0.0, raw[3]};
//				} else {
//					return new double[] {raw[0], raw[1], 0.0, raw[3] + 1};
//				}
//			}
//
//			@Override
//			public AdapterType getAdapterType() {
//				return AdapterType.VERTEX_ADAPTER;
//			}
//			
//		};
//		Adapter texPosAdapter = new CoordinateAdapter2Ifs<CoVertex> () {
//
//			@Override
//			public AdapterType getAdapterType() {
//				return AdapterType.VERTEX_ADAPTER;
//			}
//
//			@Override
//			public double[] getCoordinate(CoVertex v) {
//				Point t = v.getTextureCoord();
//				double[] raw = new double[] {t.x(), t.y(), 0.0, t.z()};
//				A.transformVector(raw);
//				if (klein) {
//					return new double[] {raw[0], raw[1], 0.0, raw[3]};
//				} else {
//					return new double[] {raw[0], raw[1], 0.0, raw[3] + 1};
//				}
//			}
//			
//		};
//		IndexedFaceSet ifs = converter.heds2ifs(
//			surface, 
//			texPosAdapter, 
//			texAdapter, 
//			cutColorAdapter, 
//			pointAdapter
//		);
//		IndexedFaceSetUtility.calculateAndSetNormals(ifs);
//		SceneGraphComponent copy = new SceneGraphComponent();
//		copy.setGeometry(ifs);
////		copiedGeometry.addChild(copy);
////		copiedGeometry.setGeometry(ifs);
////		copiedGeometry.setVisible(true);
//		surfaceRoot.addChild(copy);
	}
	
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		hif.addGlobalAdapter(new CoPositionAdapter(), true);
		hif.addGlobalAdapter(texturePositionAdapter, true);
		hif.addSelectionListener(this);
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
		c.storeProperty(getClass(), "showUniversalCover", showUniversalCover.isSelected());
		c.storeProperty(getClass(), "showFundamentalPolygon", showFundamentalPolygon.isSelected());
		c.storeProperty(getClass(), "quantizationMode", quantizationModeCombo.getSelectedIndex());
		c.storeProperty(getClass(), "useProjectiveTexture", useProjectiveTexture.isSelected());
		c.storeProperty(getClass(), "toleranceExponent", toleranceExpModel.getNumber());
		c.storeProperty(getClass(), "maxIterations", maxIterationsModel.getNumber());
		c.storeProperty(getClass(), "boundaryPanelShrinked", boundaryPanel.isShrinked());	
		c.storeProperty(getClass(), "conesPanelShrinked", coneConfigPanel.isShrinked());	
		c.storeProperty(getClass(), "visualizationPanelShrinked", visualizationPanel.isShrinked());	
		c.storeProperty(getClass(), "modelPanelShrinked", modelPanel.isShrinked());	
		c.storeProperty(getClass(), "customVertexPanelShrinked", modelPanel.isShrinked());
	} 
	
 
	@Override
	public void restoreStates(Controller c) throws Exception { 
		super.restoreStates(c);
		numConesModel.setValue(c.getProperty(getClass(), "numCones", numConesModel.getNumber().intValue()));
		numericsCombo.setSelectedIndex(c.getProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex()));
		kleinButton.setSelected(c.getProperty(getClass(), "klein", kleinButton.isSelected()));
		poincareButton.setSelected(!kleinButton.isSelected());
		showUnwrapped.setSelected(c.getProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected()));
		showUniversalCover.setSelected(c.getProperty(getClass(), "showUniversalCover", showUniversalCover.isSelected()));
		showFundamentalPolygon.setSelected(c.getProperty(getClass(), "showFundamentalPolygon", showFundamentalPolygon.isSelected()));
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
		customVertexPanel.setShrinked(c.getProperty(getClass(), "customVertexPanelShrinked", customVertexPanel.isShrinked()));
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
