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
import static java.awt.GridBagConstraints.RELATIVE;
import static java.awt.GridBagConstraints.REMAINDER;
import static java.lang.Math.PI;
import static javax.swing.JOptionPane.ERROR_MESSAGE;
import geom3d.Point;

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
import java.util.concurrent.ExecutionException;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;

import de.jreality.geometry.IndexedFaceSetFactory;
import de.jreality.geometry.IndexedFaceSetUtility;
import de.jreality.geometry.Primitives;
import de.jreality.math.Matrix;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.experimental.ManagedContent;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.tool.Tool;
import de.jreality.shader.ImageData;
import de.jreality.shader.Texture2D;
import de.jreality.shader.TextureUtility;
import de.jtem.halfedge.algorithm.triangulation.Triangulator;
import de.jtem.halfedge.jreality.ConverterHeds2JR;
import de.jtem.halfedge.jreality.adapter.Adapter;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.jtem.halfedge.jreality.adapter.TextCoordsAdapter2Ifs;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.adapter.MarkedEdgesAdapter;
import de.varylab.discreteconformal.adapter.PointAdapter;
import de.varylab.discreteconformal.adapter.PositionAdapter;
import de.varylab.discreteconformal.adapter.PositionTexCoordAdapter;
import de.varylab.discreteconformal.adapter.TexCoordAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.tasks.Unwrap;
import de.varylab.discreteconformal.plugin.tasks.Unwrap.QuantizationMode;
import de.varylab.discreteconformal.util.FundamentalDomainUtility;
import de.varylab.discreteconformal.util.UniformizationUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UniformizationUtility.FundamentalPolygon;
import de.varylab.jrworkspace.plugin.Controller;
import de.varylab.jrworkspace.plugin.PluginInfo;
import de.varylab.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.varylab.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin implements ActionListener, PropertyChangeListener {

	private static int
		polyResolution = 1000,
		coverRecursion = 2,
		coverResolution = 1000;
	
	// plug-in section ------------------ 
	private ManagedContent
		managedContent = null;
	private HalfedgeConnectorPlugin
		hcp = null;
	
	// data section ---------------------
	private CoHDS
		surface = null;
	private CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	private FundamentalPolygon 
		fundamentalPolygon = null;
	private Matrix 
		polygonTextureMatrix = euclidean().translate(-0.5, -0.5, 0).scale(0.5).scale(1, -1, 1).getMatrix();
	private Image
		polygonImage = null;
	private int
		genus = -1;
	
	private Triangulator<CoVertex, CoEdge, CoFace>
		triangulator = new Triangulator<CoVertex, CoEdge, CoFace>();

	private MarkedEdgesAdapter
		cutColorAdapter = new MarkedEdgesAdapter();
	private PointAdapter
		pointAdapter = new PointAdapter();
	
	private Tool
		hyperbolicCopyTool = new HyperbolicCopyTool(this);
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
	private JPanel
		coneConfigPanel = new JPanel(),
		modelPanel = new JPanel(),
		visualizationPanel = new JPanel();
	private SpinnerNumberModel
		numConesModel = new SpinnerNumberModel(0, 0, 100, 1),
		toleranceExpModel = new SpinnerNumberModel(-8, -30, -1, 1),
		maxIterationsModel = new SpinnerNumberModel(150, 1, 10000, 1);
	private JSpinner
		numConesSpinner = new JSpinner(numConesModel),
		toleranceExpSpinner = new JSpinner(toleranceExpModel),
		maxIterationsSpinner = new JSpinner(maxIterationsModel);
	private JCheckBox
		showGeometry = new JCheckBox("Geometry", true),
		showFundamentalPolygon = new JCheckBox("Fundamental Polygon"),
		showPoygonTexture = new JCheckBox("Polygon Texture"),
		quantizeChecker = new JCheckBox("Quantize Cone Angles"),
		useProjectiveTexture = new JCheckBox("Projective Texture", true),
		showUnwrapped = new JCheckBox("Show Unwrapped"),
		showUniversalCover = new JCheckBox("Universal Cover"),
		showUnitCircle = new JCheckBox("Show Unit Cirlce");
	private JComboBox
		quantizationModeCombo = new JComboBox(QuantizationMode.values());
	private JRadioButton
		kleinButton = new JRadioButton("Klein"),
		poincareButton = new JRadioButton("Poincaré", true),
		halfplaneButton = new JRadioButton("Half-Plane"); 
	private JComboBox
		numericsCombo = new JComboBox(new String[] {"Java/MTJ Numerics", "Petsc/Tao Numerics"});
	
	
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
		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(2,2,2,2);
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 1.0;
		c.gridwidth = GridBagConstraints.REMAINDER;
		
		coneConfigPanel.setBorder(BorderFactory.createTitledBorder("Cone Singularities"));
		coneConfigPanel.setLayout(new GridBagLayout());
		shrinkPanel.add(coneConfigPanel, c);
		
		JLabel numConesLabel = new JLabel("Cones");
		c.weightx = 0.0;
		c.gridwidth = RELATIVE;
		coneConfigPanel.add(numConesLabel, c);
		c.weightx = 1.0;
		c.gridwidth = REMAINDER;
		coneConfigPanel.add(numConesSpinner, c);
		quantizeChecker.setSelected(true);
		coneConfigPanel.add(quantizeChecker, c);
		quantizationModeCombo.setSelectedIndex(0);
		coneConfigPanel.add(quantizationModeCombo, c);
		
		numericsCombo.setLightWeightPopupEnabled(true);
		numericsCombo.setSelectedIndex(0);
		shrinkPanel.add(numericsCombo, c);
		c.gridwidth = RELATIVE;
		c.weightx = 0.0;
		shrinkPanel.add(new JLabel("Tolerance Exp"), c);
		c.gridwidth = REMAINDER;
		c.weightx = 1.0;
		shrinkPanel.add(toleranceExpSpinner, c);
		c.gridwidth = RELATIVE;
		c.weightx = 0.0;
		shrinkPanel.add(new JLabel("Max Iterations"), c);
		c.gridwidth = REMAINDER;
		c.weightx = 1.0;
		shrinkPanel.add(maxIterationsSpinner, c);
		shrinkPanel.add(unwrapBtn, c);
		
		visualizationPanel.setBorder(BorderFactory.createTitledBorder("Visualization"));
		visualizationPanel.setLayout(new GridBagLayout());
		visualizationPanel.add(showGeometry, c);
		visualizationPanel.add(showUnitCircle, c);
		visualizationPanel.add(showUnwrapped, c);
		visualizationPanel.add(showFundamentalPolygon, c);
		visualizationPanel.add(showPoygonTexture, c);
		visualizationPanel.add(showUniversalCover, c);
		visualizationPanel.add(useProjectiveTexture, c);
		shrinkPanel.add(visualizationPanel, c);
		
		modelPanel.setBorder(BorderFactory.createTitledBorder("Model"));
		modelPanel.setLayout(new GridLayout(3, 1, 2, 2));
		modelPanel.add(kleinButton);
		modelPanel.add(poincareButton);
		modelPanel.add(halfplaneButton);
		shrinkPanel.add(modelPanel, c);
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
				JOptionPane.showMessageDialog(w, msg, name, ERROR_MESSAGE);
				return;
			}
			genus = unwrapper.genus;
			if (genus > 0) {
				cutInfo = unwrapper.cutInfo;
				cutColorAdapter.setContext(cutInfo);
				pointAdapter.setContext(cutInfo);
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
			updateSurface(true);
			updateStates();
		}
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object s = e.getSource();
		if (useProjectiveTexture == s) {
			
		}
		if (unwrapBtn == s) {
			CoHDS surface = getLoaderGeometry();
			Unwrap unwrapper = new Unwrap(surface);
			unwrapper.setToleranceExponent(toleranceExpModel.getNumber().intValue());
			unwrapper.setMaxIterations(maxIterationsModel.getNumber().intValue());
			unwrapper.setNumCones(numConesModel.getNumber().intValue());
			unwrapper.setQuantizeCones(quantizeChecker.isSelected());
			unwrapper.setQuantizationMode((QuantizationMode)quantizationModeCombo.getSelectedItem());
			unwrapper.setUsePetsc(numericsCombo.getSelectedIndex() == 1);
			unwrapper.addPropertyChangeListener(this);
			unwrapper.execute();
		}
		if (showUnwrapped == s || kleinButton == s || poincareButton == s || halfplaneButton == s) {
			copiedGeometry.setGeometry(null);
			updateSurface(showUnwrapped == s);
			if (genus > 1) {
				updateFundamentalPolygon(polyResolution);
				updatePolygonTexture(coverRecursion, coverResolution);
			}
		}
		if (unwrapBtn != s) {
			updateStates();
		}
	}
	
	/*
	 * Set up scene
	 */
	private void updateStates() {
		managedContent.addContentUnique(getClass(), auxGeometry);
		managedContent.addContentUnique(getClass(), surfaceRoot);
		managedContent.addContentUnique(getClass(), copiedGeometry);
		managedContent.addContentUnique(getClass(), fundamentalPolygonRoot);
		managedContent.addContentUnique(getClass(), universalCoverRoot);
		managedContent.addToolUnique(getClass(), hyperbolicCopyTool);
		managedContent.update();
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
		hcp.getHalfedgeContent(surface, new PositionAdapter());
		if (surface.numVertices() == 0) {
			return null;
		}
		triangulator.triangulate(surface);
		surface.normalizeCoordinates();
		return surface;
	}
	
	
	private void updateSurface(final boolean align) {
		if (surface == null) {
			return;
		}
		boolean projective = useProjectiveTexture.isSelected();
		Adapter texAdapter = new TexCoordAdapter(getSelectedModel(), projective);
		Adapter posAdapter = null;
		if (showUnwrapped.isSelected()) {
			posAdapter = new PositionTexCoordAdapter(getSelectedModel(), projective);
		} else {
			posAdapter = new PositionAdapter();
		}
		IndexedFaceSet ifs = null;
		if (genus >= 1) {
			ifs = hcp.toIndexedFaceSet(surface, true, posAdapter, texAdapter, cutColorAdapter, pointAdapter);
		} else {
			ifs = hcp.toIndexedFaceSet(surface, true, posAdapter, texAdapter);
		}
		surfaceRoot.setGeometry(ifs);
		surfaceRoot.setVisible(true);
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
		if (surface == null) {
			return;
		}
		final boolean klein = kleinButton.isSelected(); 
		CoEdge edge = null;
		int i = 0;
		for (CoEdge e : surface.getEdges()) {
			if (e.isPositive()) {
				if (i == edgeIndex) {
					if (e.getLeftFace() != null) {
						edge = e;
					} else {
						edge = e.getOppositeEdge();
					}
					break;
				}
				i++;
			}
		}
		
		if (edge == null) {
			System.err.println("Edge corresponding to index " + edgeIndex + " not found");
			return;
		}
		CoEdge coEdge = cutInfo.edgeCutMap.get(edge);
		if (coEdge == null) {
			System.err.println("CoEdge not found");
			return;
		}
		if (edge.getOppositeEdge().getLeftFace() != null) {
			System.err.println("Picked no boundary edge!");
			return;
		}
		System.out.println("hyperbolic motion: " + edge + " -> " + coEdge);
		
		final Matrix A = UniformizationUtility.makeHyperbolicMotion(coEdge, edge.getOppositeEdge());

		ConverterHeds2JR<CoVertex, CoEdge, CoFace>
			converter = new ConverterHeds2JR<CoVertex, CoEdge, CoFace>();
		Adapter texAdapter = new TextCoordsAdapter2Ifs<CoVertex> () {

			@Override
			public double[] getTextCoordinate(CoVertex v) {
				Point t = v.getTextureCoord();
				double[] raw = new double[] {t.x(), t.y(), 0.0, t.z()};
				A.transformVector(raw);
				if (klein) {
					return new double[] {raw[0], raw[1], 0.0, raw[3]};
				} else {
					return new double[] {raw[0], raw[1], 0.0, raw[3] + 1};
				}
			}

			@Override
			public AdapterType getAdapterType() {
				return AdapterType.VERTEX_ADAPTER;
			}
			
		};
		Adapter texPosAdapter = new CoordinateAdapter2Ifs<CoVertex> () {

			@Override
			public AdapterType getAdapterType() {
				return AdapterType.VERTEX_ADAPTER;
			}

			@Override
			public double[] getCoordinate(CoVertex v) {
				Point t = v.getTextureCoord();
				double[] raw = new double[] {t.x(), t.y(), 0.0, t.z()};
				A.transformVector(raw);
				if (klein) {
					return new double[] {raw[0], raw[1], 0.0, raw[3]};
				} else {
					return new double[] {raw[0], raw[1], 0.0, raw[3] + 1};
				}
			}
			
		};
		IndexedFaceSet ifs = converter.heds2ifs(
			surface, 
			texPosAdapter, 
			texAdapter, 
			cutColorAdapter, 
			pointAdapter
		);
		IndexedFaceSetUtility.calculateAndSetNormals(ifs);
		copiedGeometry.setGeometry(ifs);
		copiedGeometry.setVisible(true);
	}
	
	
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
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "numCones", numConesModel.getNumber().intValue());
		c.storeProperty(getClass(), "quantizeCones", quantizeChecker.isSelected());
		c.storeProperty(getClass(), "quantizeMode", quantizationModeCombo.getSelectedItem());
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
	} 
	
 
	@Override
	public void restoreStates(Controller c) throws Exception { 
		super.restoreStates(c);
		numConesModel.setValue(c.getProperty(getClass(), "numCones", numConesModel.getNumber().intValue()));
		quantizeChecker.setSelected(c.getProperty(getClass(), "quantizeCones", quantizeChecker.isSelected()));
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
