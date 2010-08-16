package de.varylab.discreteconformal.plugin;

import static de.jreality.math.MatrixBuilder.euclidean;
import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.LIGHTING_ENABLED;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.VERTEX_COLORS_ENABLED;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.varylab.discreteconformal.util.FundamentalDomainUtility.createFundamentalPolygon;
import static de.varylab.discreteconformal.util.UniformizationUtility.constructFundamentalPolygon;
import static java.awt.Color.WHITE;
import static java.lang.Math.PI;
import static javax.swing.JOptionPane.ERROR_MESSAGE;
import static javax.swing.JOptionPane.OK_CANCEL_OPTION;
import geom3d.Point;

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
import de.jreality.plugin.basic.View;
import de.jreality.scene.Appearance;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.shader.ImageData;
import de.jreality.shader.Texture2D;
import de.jreality.shader.TextureUtility;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
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
import de.varylab.discreteconformal.heds.adapter.MarkedEdgesColorAdapter;
import de.varylab.discreteconformal.heds.adapter.MarkedEdgesRadiusAdapter;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.adapter.TexCoordAdapter;
import de.varylab.discreteconformal.heds.adapter.TexCoordPositionAdapter;
import de.varylab.discreteconformal.heds.calculator.SubdivisionCalculator;
import de.varylab.discreteconformal.plugin.tasks.Unwrap;
import de.varylab.discreteconformal.unwrapper.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.unwrapper.UnwrapUtility.QuantizationMode;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.FundamentalDomainUtility;
import de.varylab.discreteconformal.util.UniformizationUtility.FundamentalPolygon;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin implements ListSelectionListener, ChangeListener, ActionListener, PropertyChangeListener, SelectionListener {

	private static int
		polyResolution = 1000,
		coverRecursion = 2,
		coverResolution = 1000;
	
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
		universalCoverAppearance = new Appearance(),
		surfaceAppearance = new Appearance();
	private SceneGraphComponent
		surfaceRoot = new SceneGraphComponent("Surface Geometry"),
		auxGeometry = new SceneGraphComponent("Aux Geometry"),
		copiedGeometry = new SceneGraphComponent("Click Copy Geometry"),
		fundamentalPolygonRoot = new SceneGraphComponent("Fundamental Polygon"), 
		unitCircle = new SceneGraphComponent("Hyperbolic Boundary"),
		universalCoverRoot = new SceneGraphComponent("Universal Cover");

	// user interface section ------------
	private JButton
		unwrapBtn = new JButton("Unwrap");
	private ShrinkPanel
		customVertexPanel = new ShrinkPanel("Custom Vertices"),
		boundaryPanel = new ShrinkPanel("Boundary"),
		coneConfigPanel = new ShrinkPanel("Automatic Cones"),
		modelPanel = new ShrinkPanel("Hyperbolic Model"),
		visualizationPanel = new ShrinkPanel("Visualization");
	private SpinnerNumberModel
		customThetaModel = new SpinnerNumberModel(360, 0, 1000, 1),
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
		showGeometry = new JCheckBox("Geometry", true),
		showFundamentalPolygon = new JCheckBox("Fundamental Polygon"),
		showPoygonTexture = new JCheckBox("Polygon Texture"),
		useProjectiveTexture = new JCheckBox("Projective Texture", true),
		showUnwrapped = new JCheckBox("Show Unwrapped"),
		showUniversalCover = new JCheckBox("Universal Cover"),
		showUnitCircle = new JCheckBox("Show Unit Cirlce");
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
		showGeometry.addActionListener(this);
		showPoygonTexture.addActionListener(this);
		unwrapBtn.addActionListener(this);
		showUnwrapped.addActionListener(this);
		showUnitCircle.addActionListener(this);
		kleinButton.addActionListener(this);
		poincareButton.addActionListener(this);
		halfplaneButton.addActionListener(this);
		showUniversalCover.addActionListener(this);
		showFundamentalPolygon.addActionListener(this);
		useProjectiveTexture.addActionListener(this);
		
		ButtonGroup modelGroup = new ButtonGroup();
		modelGroup.add(kleinButton);
		modelGroup.add(poincareButton);
		modelGroup.add(halfplaneButton);
		
		Appearance circleApp = new Appearance();
		circleApp.setAttribute(EDGE_DRAW, false);
		circleApp.setAttribute(VERTEX_DRAW, false); 
		circleApp.setAttribute(FACE_DRAW, true);
		unitCircle.setAppearance(circleApp);
		unitCircle.setVisible(false);
		euclidean().rotate(PI / 2, 1, 0, 0).assignTo(unitCircle);
		unitCircle.setGeometry(Primitives.torus(1.0, 0.005, 200, 5));
		auxGeometry.addChild(unitCircle);
		
		surfaceAppearance.setAttribute(POLYGON_SHADER + "." + VERTEX_COLORS_ENABLED, false);
		surfaceRoot.setAppearance(surfaceAppearance);
		
		Appearance funPolyApp = new Appearance();
		funPolyApp.setAttribute(VERTEX_DRAW, false);
		funPolyApp.setAttribute(EDGE_DRAW, true);
		fundamentalPolygonRoot.setAppearance(funPolyApp);
		
		IndexedFaceSetFactory ifsf = new IndexedFaceSetFactory();
		ifsf.setVertexCount(4);
		ifsf.setFaceCount(1);
		ifsf.setFaceIndices(new int[][] {{0,1,2,3}});
		ifsf.setVertexTextureCoordinates(new double[] {-1,-1,1,-1,1,1,-1,1});
		ifsf.setVertexCoordinates(new double[]{-1,-1,-0.01, 1,-1,-0.01,  1,1,-0.001,  -1,1,-0.01});
		ifsf.setGenerateFaceNormals(true);
		ifsf.update();
		universalCoverRoot.setGeometry(ifsf.getGeometry());
		universalCoverRoot.setAppearance(surfaceAppearance);
		universalCoverRoot.setAppearance(universalCoverAppearance);
		universalCoverAppearance.setAttribute(VERTEX_DRAW, false);
		universalCoverAppearance.setAttribute(EDGE_DRAW, false);
		universalCoverAppearance.setAttribute(FACE_DRAW, true);
		universalCoverAppearance.setAttribute(LIGHTING_ENABLED, false);
		universalCoverAppearance.setAttribute(DIFFUSE_COLOR, WHITE);
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
		visualizationPanel.add(showGeometry, c2);
		visualizationPanel.add(showUnitCircle, c2);
		visualizationPanel.add(showUnwrapped, c2);
		visualizationPanel.add(showFundamentalPolygon, c2);
		visualizationPanel.add(showPoygonTexture, c2);
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
	}
	
	
	@Override
	public void selectionChanged(HalfedgeSelection s, HalfedgeInterface hif) {
		Set<CoVertex> oldSelection = new HashSet<CoVertex>(customVertices);
		customVertices.clear();
		for (Vertex<?,?,?> v : s.getVertices()) {
			if (v instanceof CoVertex) {
				CoVertex cov = (CoVertex)v;
				if (cov.getCustomInfo() == null) {
					cov.setCustomInfo(new CustomVertexInfo());
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
		customModeCombo.setSelectedItem(v.getCustomInfo().boundaryMode);
		customQuantizationCombo.setSelectedItem(v.getCustomInfo().quantizationMode);
		useCustomThetaChecker.setSelected(v.getCustomInfo().useCustomTheta);
		customThetaModel.setValue(Math.toDegrees(v.getCustomInfo().theta));
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (selectedVertexList.getSelectedValue() == null) return;
		for (Object s : selectedVertexList.getSelectedValues()) {
			CoVertex v = (CoVertex)s;
			double thetaDeg = customThetaModel.getNumber().doubleValue();
			v.getCustomInfo().theta = Math.toRadians(thetaDeg);
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
			}
			genus = unwrapper.genus;
			if (genus > 0) {
				cutInfo = unwrapper.cutInfo;
				cutColorAdapter.setContext(cutInfo);
				cutRadiusAdapter.setContext(cutInfo);
				pointRadiusAdapter.setContext(cutInfo);
				pointColorAdapter.setContext(cutInfo);
			}
			if (unwrapper.genus > 1) {
				fundamentalPolygon = constructFundamentalPolygon(cutInfo);
				updateFundamentalPolygon(polyResolution);
				updatePolygonTexture(coverRecursion, coverResolution);
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
		if (
			useProjectiveTexture == s ||
			showFundamentalPolygon == s ||
			showGeometry == s ||
			showPoygonTexture == s ||
			showUnitCircle == s ||
			showUniversalCover == s
		) {
			updateStates();
			return;
		}
		if (unwrapBtn == s) {
			CoHDS surface = getLoaderGeometry();
			Unwrap uw = new Unwrap(surface);
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
		if (showUnwrapped == s || kleinButton == s || poincareButton == s || halfplaneButton == s) {
			copiedGeometry.setGeometry(null);
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
				v.getCustomInfo().boundaryMode = (BoundaryMode)customModeCombo.getSelectedItem();
			}
		}
		if (customQuantizationCombo == s) {
			for (Object sel : selectedVertexList.getSelectedValues()) {
				CoVertex v = (CoVertex)sel;
				v.getCustomInfo().quantizationMode = (QuantizationMode)customQuantizationCombo.getSelectedItem();
			}
		}
		if (useCustomThetaChecker == s) {
			for (Object sel : selectedVertexList.getSelectedValues()) {
				CoVertex v = (CoVertex)sel;
				v.getCustomInfo().useCustomTheta = useCustomThetaChecker.isSelected(); 
			}
		}
	}
	
	/*
	 * Set up scene
	 */
	private void updateStates() {
//		managedContent.addContentUnique(getClass(), auxGeometry);
//		managedContent.addContentUnique(getClass(), surfaceRoot);
//		managedContent.addContentUnique(getClass(), copiedGeometry);
//		managedContent.addContentUnique(getClass(), fundamentalPolygonRoot);
//		managedContent.addContentUnique(getClass(), universalCoverRoot);
////		managedContent.addToolUnique(getClass(), hyperbolicCopyTool);
//		managedContent.update();
		surfaceRoot.setVisible(showGeometry.isSelected());
		if (genus > 1) {
			fundamentalPolygonRoot.setVisible(showUnwrapped.isSelected() && showFundamentalPolygon.isSelected());
			universalCoverRoot.setVisible(showUniversalCover.isSelected() && showUnwrapped.isSelected());
			unitCircle.setVisible(showUnitCircle.isSelected() && showUnwrapped.isSelected() && (getSelectedModel() != HyperbolicModel.Halfplane));
		} else {
			unitCircle.setVisible(false);
			fundamentalPolygonRoot.setVisible(false);
			universalCoverRoot.setVisible(false);
		}
		ImageData imgData = null;
		if (polygonImage != null && (showPoygonTexture.isSelected() || showUniversalCover.isSelected())) {
			imgData = new ImageData(polygonImage);
		}
		if (polygonImage != null && showPoygonTexture.isSelected()) {
			Texture2D tex2d = TextureUtility.createTexture(surfaceAppearance, POLYGON_SHADER, imgData);
			tex2d.setTextureMatrix(polygonTextureMatrix);
		} else {
			TextureUtility.removeTexture(surfaceAppearance, POLYGON_SHADER);
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
				v.setCustomInfo(null);
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
//		ConverterHeds2JR converter = new ConverterHeds2JR();
		boolean projective = useProjectiveTexture.isSelected();
		AdapterSet adapters = new AdapterSet();
		adapters.add(new TexCoordAdapter(getSelectedModel(), projective));
		if (showUnwrapped.isSelected()) {
			adapters.add(new TexCoordPositionAdapter(getSelectedModel(), projective));
		} else {
			adapters.add(new PositionAdapter());
		}
//		IndexedFaceSet ifs = null;
		if (genus >= 1) {
			adapters.add(cutRadiusAdapter);
			adapters.add(cutColorAdapter);
			adapters.add(pointRadiusAdapter);
			adapters.add(pointColorAdapter);
//			ifs = converter.heds2ifs(surface, adapters, null);
		} else {
//			ifs = converter.heds2ifs(surface, adapters, null);
		}
		hif.set(surface, adapters);
		hif.getActiveLayer().addTemporaryGeometry(auxGeometry);
//		hif.getActiveLayer().addTemporaryGeometry(surfaceRoot);
		hif.getActiveLayer().addTemporaryGeometry(copiedGeometry);
		hif.getActiveLayer().addTemporaryGeometry(fundamentalPolygonRoot);
		hif.getActiveLayer().addTemporaryGeometry(universalCoverRoot);
		hif.encompassAll();
//		IndexedFaceSetUtility.calculateAndSetNormals(ifs);
//		surfaceRoot.setGeometry(ifs);
//		surfaceRoot.setVisible(true);
	}
	

	public void updateFundamentalPolygon(int resolution) {
		createFundamentalPolygon(fundamentalPolygon, cutInfo, fundamentalPolygonRoot, resolution, getSelectedModel());
	}
	
	
	
	public void updatePolygonTexture(int depth, int resolution) {
		Point pRoot = cutInfo.cutRoot.getTextureCoord();
		double[] root = new double[] {pRoot.x(), pRoot.y(), 0.0, pRoot.z()};
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
		hif.addAdapter(new PositionAdapter());
		hif.addAdapter(new TexCoordAdapter(0));
		hif.addCalculator(new SubdivisionCalculator());
		hif.addSelectionListener(this);
//		managedContent = c.getPlugin(ManagedContent.class);
	}
	
	@Override
	public void uninstall(Controller c) throws Exception {
		super.uninstall(c);
//		managedContent.removeAll(getClass());
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
		c.storeProperty(getClass(), "showGeometry", showGeometry.isSelected());
		c.storeProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected());
		c.storeProperty(getClass(), "showUnitCircle", showUnitCircle.isSelected());
		c.storeProperty(getClass(), "showPoygonTexture", showPoygonTexture.isSelected());
		c.storeProperty(getClass(), "showUniversalCover", showUniversalCover.isSelected());
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
		showGeometry.setSelected(c.getProperty(getClass(), "showGeometry", showGeometry.isSelected()));
		showUnwrapped.setSelected(c.getProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected()));
		showUnitCircle.setSelected(c.getProperty(getClass(), "showUnitCircle", showUnitCircle.isSelected()));
		showPoygonTexture.setSelected(c.getProperty(getClass(), "showPoygonTexture", showPoygonTexture.isSelected()));
		showUniversalCover.setSelected(c.getProperty(getClass(), "showUniversalCover", showUniversalCover.isSelected()));
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
