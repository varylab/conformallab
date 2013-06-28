package de.varylab.discreteconformal.plugin;

import static de.jreality.math.MatrixBuilder.euclidean;
import static de.jreality.scene.Appearance.DEFAULT;
import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.LIGHTING_ENABLED;
import static de.jreality.shader.CommonAttributes.LINE_SHADER;
import static de.jreality.shader.CommonAttributes.POINT_SHADER;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.SPHERES_DRAW;
import static de.jreality.shader.CommonAttributes.TEXTURE_2D;
import static de.jreality.shader.CommonAttributes.TRANSPARENCY;
import static de.jreality.shader.CommonAttributes.TRANSPARENCY_ENABLED;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.varylab.discreteconformal.adapter.HyperbolicModel.Klein;
import static de.varylab.discreteconformal.adapter.HyperbolicModel.Poincaré;
import static de.varylab.discreteconformal.uniformization.VisualizationUtility.drawTriangulation;
import static de.varylab.discreteconformal.uniformization.VisualizationUtility.drawUniversalCoverImage;
import static de.varylab.discreteconformal.util.UnwrapUtility.prepareInvariantDataEuclidean;
import static java.awt.Color.BLACK;
import static java.awt.Color.WHITE;
import static java.awt.Color.YELLOW;
import static java.awt.image.BufferedImage.TYPE_INT_ARGB;
import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static javax.swing.JOptionPane.ERROR_MESSAGE;
import static javax.swing.JOptionPane.OK_CANCEL_OPTION;
import static javax.swing.JOptionPane.WARNING_MESSAGE;
import static javax.swing.JOptionPane.showMessageDialog;
import static javax.swing.SwingUtilities.getWindowAncestor;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.apache.batik.svggen.SVGGraphics2D;

import com.itextpdf.text.Document;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;

import de.jreality.geometry.IndexedFaceSetFactory;
import de.jreality.geometry.Primitives;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.content.ContentAppearance;
import de.jreality.plugin.job.Job;
import de.jreality.plugin.job.JobListener;
import de.jreality.plugin.job.JobQueuePlugin;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.PointSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.shader.ImageData;
import de.jreality.shader.Texture2D;
import de.jreality.shader.TextureUtility;
import de.jreality.ui.AppearanceInspector;
import de.jreality.ui.ColorChooseJButton;
import de.jreality.ui.ColorChooseJButton.ColorChangedEvent;
import de.jreality.ui.ColorChooseJButton.ColorChangedListener;
import de.jreality.ui.TextureInspector;
import de.jreality.ui.viewerapp.FileFilter;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.halfedgetools.jreality.ConverterHeds2JR;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeSelection;
import de.jtem.halfedgetools.plugin.SelectionListener;
import de.jtem.halfedgetools.plugin.image.ImageHook;
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
import de.varylab.discreteconformal.heds.adapter.MarkedEdgesColorAdapter;
import de.varylab.discreteconformal.heds.adapter.MarkedEdgesRadiusAdapter;
import de.varylab.discreteconformal.heds.adapter.MetricErrorAdapter;
import de.varylab.discreteconformal.plugin.tasks.Unwrap;
import de.varylab.discreteconformal.uniformization.CanonicalFormUtility;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility;
import de.varylab.discreteconformal.uniformization.FundamentalVertex;
import de.varylab.discreteconformal.uniformization.SurfaceCurveUtility;
import de.varylab.discreteconformal.uniformization.VisualizationUtility;
import de.varylab.discreteconformal.unwrapper.BoundaryMode;
import de.varylab.discreteconformal.unwrapper.QuantizationMode;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin 
	implements ListSelectionListener, ChangeListener, ActionListener, SelectionListener, ColorChangedListener, JobListener {

	private static int
		coverResolution = 1024;
	
	private enum Domain {
		Cut,
		Minimal,
		Opposite,
		Canonical
	}
	
	// plug-in section ------------------ 
	private HalfedgeInterface
		hif = null;
	private ConformalVisualizationPlugin
		vis = null;
	private ContentAppearance
		contentAppearance = null;
	private DomainVisualisationPlugin
		domainVisualisationPlugin = null;
	private JobQueuePlugin
		jobQueue = null;
	
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
		oppositePolygon = null,
		canonicalPolygon = null;
	private Matrix 
		polygonTextureMatrix = euclidean().translate(-0.5, -0.5, 0).scale(0.5).scale(1, -1, 1).getMatrix();
	private int
		genus = -1;
	private CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;

	private MetricErrorAdapter
		metricErrorAdapter = new MetricErrorAdapter();
	public CoPositionAdapter
		positionAdapter = new CoPositionAdapter();
	public CoTexturePositionAdapter
		texturePositionAdapter = new CoTexturePositionAdapter();
	private MarkedEdgesColorAdapter
		cutColorAdapter = new MarkedEdgesColorAdapter();
	private MarkedEdgesRadiusAdapter
		cutRadiusAdapter = new MarkedEdgesRadiusAdapter();
	private BranchPointColorAdapter
		pointColorAdapter = new BranchPointColorAdapter();
	private BranchPointRadiusAdapter
		pointRadiusAdapter = new BranchPointRadiusAdapter();
	
	private Appearance
		yellowPointsAppearance = new Appearance(),
		polygonCurvesAppearance = new Appearance(),
		axesCurvesAppearance = new Appearance(),
		universalCoverAppearance = new Appearance();
	private SceneGraphComponent
		selectedCustomNodesRoot = new SceneGraphComponent("Selected Custom Nodes"),
		polygonCurvesRoot = new SceneGraphComponent("Polygon Curves"),
		axesCurvesRoot = new SceneGraphComponent("Axes Curves"),
		unitCircle = new SceneGraphComponent("Hyperbolic Boundary"),
		domainRoot = new SceneGraphComponent("Domain");
	
	private BufferedImage
		activeDomainImage = null;
		

	// user interface section ------------
	private JPanel
		visButtonsPanel = new JPanel();
	private JButton
		saveTextureButton = new JButton("Save Texture"),
		exportHyperbolicButton = new JButton(ImageHook.getIcon("disk.png")),
		exportHyperbolicSVGButton = new JButton(ImageHook.getIcon("disk.png")),
		moveToCenterButton = new JButton("Center Selected Vertex"),
		coverToTextureButton = new JButton("Create Texture"),
		checkGaussBonnetBtn = new JButton("Check Gauß-Bonnet"),
		unwrapBtn = new JButton("Unwrap"),
		spherizeButton = new JButton("Spherize"),
		quantizeToQuads = new JButton("Quads");
	private ColorChooseJButton
		triangulationColorButton = new ColorChooseJButton(Color.GRAY, true),
		polygonColorButton = new ColorChooseJButton(Color.RED, true),
		axesColorButton = new ColorChooseJButton(Color.BLUE, true);
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
		coverRecursionModel = new SpinnerNumberModel(0, 0, 10, 1),
		customThetaModel = new SpinnerNumberModel(360.0, 0.0, 10000.0, 1.0),
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
		rescaleChecker = new JCheckBox("Rescale Geometry", true),
		circularEdgeChecker = new JCheckBox("Is Circular Edge"), 
		useDistanceToCanonicalize = new JCheckBox("Use Isometry Distances"),
		useCustomThetaChecker = new JCheckBox("Custom Theta"),
		useProjectiveTexture = new JCheckBox("Projective Texture", true),
		drawTriangulationChecker = new JCheckBox("Draw Triangulation"),
		drawAxesChecker = new JCheckBox("Draw Axes"),
		drawPolygonChecker = new JCheckBox("Draw Polygon", true),
		drawCurvesOnSurface = new JCheckBox("Draw Curves On Surface");
	private JComboBox
		numericsCombo = new JComboBox(new String[] {"Java/MTJ Numerics", "Petsc/Tao Numerics"}),
		conesQuantizationModeCombo = new JComboBox(QuantizationMode.values()),
		boundaryModeCombo = new JComboBox(BoundaryMode.values()),
		boundaryQuantizationCombo = new JComboBox(QuantizationMode.values()),
		customModeCombo = new JComboBox(BoundaryMode.values()),
		customQuantizationCombo = new JComboBox(QuantizationMode.values());
	private JList
		customNodesList = new JList();
	private JScrollPane
		selectionScroller = new JScrollPane(customNodesList);
	private JFileChooser
		pngChooser = new JFileChooser(),
		pdfChooser = new JFileChooser(),
		svgChooser = new JFileChooser();
		
	public DiscreteConformalPlugin() {
		createLayout();
		unwrapBtn.addActionListener(this);
		checkGaussBonnetBtn.addActionListener(this);
		domainCombo.addActionListener(this);
		useProjectiveTexture.addActionListener(this);
		coverToTextureButton.addActionListener(this);
		quantizeToQuads.addActionListener(this);
		
		IndexedFaceSetFactory ifsf = new IndexedFaceSetFactory();
		ifsf.setVertexCount(4);
		ifsf.setFaceCount(1);
		ifsf.setFaceIndices(new int[][] {{0,1,2,3}});
		ifsf.setVertexTextureCoordinates(new double[] {-1,-1,1,-1,1,1,-1,1});
		ifsf.setVertexCoordinates(new double[]{-1,-1,0.001, 1,-1,0.001, 1,1,0.001, -1,1,0.001});
		ifsf.setGenerateFaceNormals(true);
		ifsf.update();
		domainRoot.setGeometry(ifsf.getGeometry());
		domainRoot.setAppearance(universalCoverAppearance);
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
		domainRoot.addChild(unitCircle);
		
		polygonCurvesAppearance.setAttribute(EDGE_DRAW, true);
		polygonCurvesRoot.setAppearance(polygonCurvesAppearance);
		axesCurvesAppearance.setAttribute(EDGE_DRAW, true);
		axesCurvesRoot.setAppearance(axesCurvesAppearance);
		
		yellowPointsAppearance.setAttribute(VERTEX_DRAW, true);
		yellowPointsAppearance.setAttribute(POINT_SHADER + "." + DIFFUSE_COLOR, YELLOW);
		yellowPointsAppearance.setAttribute(POINT_SHADER + "." + SPHERES_DRAW, true);
		selectedCustomNodesRoot.setAppearance(yellowPointsAppearance);
		
		pngChooser = new JFileChooser();
		pngChooser.addChoosableFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "PNG Files (*.png)";
			}
			@Override
			public boolean accept(File f) {
				return f.isDirectory() || f.getName().toLowerCase().endsWith(".png");
			}
		});
		pdfChooser = new JFileChooser();
		pdfChooser.addChoosableFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "PDF Files (*.pdf)";
			}
			@Override
			public boolean accept(File f) {
				return f.isDirectory() || f.getName().toLowerCase().endsWith(".pdf");
			}
		});
		svgChooser = new JFileChooser();
		svgChooser.addChoosableFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "Scalable Vector Graphics (*.svg)";
			}
			@Override
			public boolean accept(File f) {
				return f.isDirectory() || f.getName().toLowerCase().endsWith(".svg");
			}
		});		
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
		shrinkPanel.add(rescaleChecker, c2);
		shrinkPanel.add(checkGaussBonnetBtn, c2);
		shrinkPanel.add(unwrapBtn, c2);
		shrinkPanel.add(spherizeButton, c2);
		
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
		coneConfigPanel.add(conesQuantizationModeCombo, c2);
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
		visualizationPanel.add(new JLabel("Domain"), c1);
		visualizationPanel.add(domainCombo, c2);
		visualizationPanel.add(new JLabel("Cover Recursion"), c1);
		visualizationPanel.add(coverRecursionSpinner, c2);
		visualizationPanel.add(drawTriangulationChecker, c1);
		visualizationPanel.add(triangulationColorButton, c2);
		visualizationPanel.add(drawPolygonChecker, c1);
		visualizationPanel.add(polygonColorButton, c2);
		visualizationPanel.add(drawAxesChecker, c1);
		visualizationPanel.add(axesColorButton, c2);
		visualizationPanel.add(drawCurvesOnSurface, c2);
		visualizationPanel.add(useDistanceToCanonicalize, c2);
		visualizationPanel.add(coverToTextureButton, c2);
		visualizationPanel.add(visButtonsPanel, c2);
		shrinkPanel.add(visualizationPanel, c2);

		visButtonsPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
		visButtonsPanel.add(saveTextureButton);
		visButtonsPanel.add(exportHyperbolicButton);
		visButtonsPanel.add(exportHyperbolicSVGButton);
		
		modelPanel.setLayout(new GridBagLayout());
		modelPanel.add(moveToCenterButton, c2);
		modelPanel.setShrinked(true);
		shrinkPanel.add(modelPanel, c2);
		
		customNodesList.getSelectionModel().addListSelectionListener(this);
		customModeCombo.addActionListener(this);
		customQuantizationCombo.addActionListener(this);
		customNodePanel.setShrinked(true);
		useCustomThetaChecker.addActionListener(this);
		customThetaSpinner.addChangeListener(this);
		circularEdgeChecker.addActionListener(this);
		moveToCenterButton.addActionListener(this);
		saveTextureButton.addActionListener(this);
		drawTriangulationChecker.addActionListener(this);
		exportHyperbolicButton.addActionListener(this);
		exportHyperbolicSVGButton.addActionListener(this);
		coverRecursionSpinner.addChangeListener(this);
		drawPolygonChecker.addActionListener(this);
		drawAxesChecker.addActionListener(this);
		drawCurvesOnSurface.addActionListener(this);
		polygonColorButton.addColorChangedListener(this);
		axesColorButton.addColorChangedListener(this);
		triangulationColorButton.addColorChangedListener(this);
		spherizeButton.addActionListener(this);
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
		updateDomainImage();
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
		}
		hif.removeTemporaryGeometry(selectedCustomNodesRoot);
		hif.addTemporaryGeometry(selectedCustomNodesRoot);
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (customThetaSpinner == e.getSource()) {
			if (customNodesList.getSelectedValue() == null) return;
			for (Object s : customNodesList.getSelectedValues()) {
				if (!(s instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)s;
				double thetaDeg = customThetaModel.getNumber().doubleValue();
				v.info.theta = Math.toRadians(thetaDeg);
			}
		}
		if (coverRecursionSpinner == e.getSource()) {
			updateDomainImage();
		}
	}
	
	@Override
	public void jobFinished(Job job) {
		Unwrap unwrapper = (Unwrap)job;
		surface = unwrapper.getSurface();
		surface.revertNormalization();
		genus = unwrapper.genus;
		cutInfo = unwrapper.cutInfo;
		metricErrorAdapter.setLengthMap(unwrapper.lengthMap);
		metricErrorAdapter.setSignature(Pn.EUCLIDEAN);
		createVisualization(surface, genus, cutInfo);
		updateSurface();
		updateDomainImage();
	}
	@Override
	public void jobFailed(Job job, Exception e) {
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
	
	public void createVisualization(CoHDS surface, int genus, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		this.genus = genus;
		this.surface = surface;
		this.cutInfo = cutInfo;
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
				System.out.println("Constructing opposites sides polygon...");
				oppositePolygon = CanonicalFormUtility.canonicalizeOpposite(minimalPolygon);
				System.out.println(oppositePolygon);
				oppositePolygon.checkRelation();	
				System.out.println("Constructing fast canonical polygon...");
				canonicalPolygon = FundamentalPolygonUtility.canonicalize(minimalPolygon, useDistanceToCanonicalize.isSelected());
				System.out.println(canonicalPolygon);
				canonicalPolygon.checkRelation();
				metricErrorAdapter.setSignature(Pn.HYPERBOLIC);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		if (genus == 0) {
			cutInfo = null;
			cuttedPolygon = null;
		}
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object s = e.getSource();
		if (domainCombo == s) {
			updateSurface();
			updateDomainImage();
			return;
		}
		if (coverToTextureButton == s) {
			if (activeDomainImage == null) {
				Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
				JOptionPane.showMessageDialog(w, "No current domain image", "Cannot create texture", WARNING_MESSAGE);
				return;
			}
			AppearanceInspector ai = contentAppearance.getAppearanceInspector();
			TextureInspector ti = ai.getTextureInspector();
			ti.addTexture("Hyperbolic Domain", activeDomainImage);
			ti.setTextureScaleLock(false);
			ti.setTextureUScale(0.5);
			ti.setTextureVScale(-0.5);
			ti.setTextureUTranslation(1.0);
			ti.setTextureVTranslation(1.0);
			ti.setTextureRotation(0.0);
			ti.setTextureShear(0.0);
		}
		if (unwrapBtn == s || spherizeButton == s) {
			CoHDS surface = getLoaderGeometry();
			if (surface == null) return;
			if (isRescaleGeometry()) {
				surface.normalizeCoordinates();
			}
			AdapterSet aSet = hif.getAdapters();
			Unwrap uw = new Unwrap(surface, aSet);
			uw.setSpherize(spherizeButton == s);
			uw.setToleranceExponent(toleranceExpModel.getNumber().intValue());
			uw.setMaxIterations(maxIterationsModel.getNumber().intValue());
			uw.setNumCones(numConesModel.getNumber().intValue());
			uw.setQuantizationMode((QuantizationMode)conesQuantizationModeCombo.getSelectedItem());
			uw.setBoundaryQuantMode((QuantizationMode)boundaryQuantizationCombo.getSelectedItem());
			uw.setBoundaryMode((BoundaryMode)boundaryModeCombo.getSelectedItem());
			uw.setUsePetsc(numericsCombo.getSelectedIndex() == 1);
			uw.setSelectedVertices(hif.getSelection().getVertices(surface));
			uw.addJobListener(this);
			jobQueue.queueJob(uw);
		}
		if (customModeCombo == s) {
			for (Object sel : customNodesList.getSelectedValues()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.boundaryMode = (BoundaryMode)customModeCombo.getSelectedItem();
			}
		}
		if (customQuantizationCombo == s) {
			for (Object sel : customNodesList.getSelectedValues()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.quantizationMode = (QuantizationMode)customQuantizationCombo.getSelectedItem();
			}
		}
		if (useCustomThetaChecker == s) {
			for (Object sel : customNodesList.getSelectedValues()) {
				if (!(sel instanceof CoVertex)) continue;
				CoVertex v = (CoVertex)sel;
				v.info.useCustomTheta = useCustomThetaChecker.isSelected(); 
			}
		}
		if (circularEdgeChecker == s) {
			for (Object sel : customNodesList.getSelectedValues()) {
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
			updateDomainImage();
		}
		if (saveTextureButton == s) {
			Window w = SwingUtilities.getWindowAncestor(this.shrinkPanel);
			int result = pngChooser.showSaveDialog(w);
			if (result != JFileChooser.APPROVE_OPTION) {
				return;
			}
			File file = pngChooser.getSelectedFile();
			if (!file.getName().toLowerCase().endsWith(".png")) {
				file = new File(file.getAbsolutePath() + ".png");
			}
			if (file.exists()) {
				result = JOptionPane.showConfirmDialog(w, "File exists, overwrite?", "File exists", OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
				if (result != JOptionPane.OK_OPTION) {
					return;
				}
			}			
			BufferedImage img = createDomainImage(coverResolution);
			try {
				ImageIO.write(img, "png", pngChooser.getSelectedFile());
			} catch (Exception e2) {
				JOptionPane.showMessageDialog(w, e2.getMessage(), "Error", ERROR_MESSAGE);
			}
		}
		if (exportHyperbolicButton == s) {
			Window w = SwingUtilities.getWindowAncestor(this.shrinkPanel);
			int result = pdfChooser.showSaveDialog(w);
			if (result != JFileChooser.APPROVE_OPTION) {
				return;
			}
			File file = pdfChooser.getSelectedFile();
			if (!file.getName().toLowerCase().endsWith(".pdf")) {
				file = new File(file.getAbsolutePath() + ".pdf");
			}
			if (file.exists()) {
				result = JOptionPane.showConfirmDialog(w, "File exists, overwrite?", "File exists", OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
				if (result != JOptionPane.OK_OPTION) {
					return;
				}
			}
			try {
				exportHyperbolicImageToPDF(file, coverResolution);
			} catch (Exception e2) {
				JOptionPane.showMessageDialog(w, e2.getMessage(), "Error", ERROR_MESSAGE);
			}
		}
		if (exportHyperbolicSVGButton == s) {
			Window w = SwingUtilities.getWindowAncestor(this.shrinkPanel);
			int result = svgChooser.showSaveDialog(w);
			if (result != JFileChooser.APPROVE_OPTION) {
				return;
			}
			File file = svgChooser.getSelectedFile();
			if (!file.getName().toLowerCase().endsWith(".svg")) {
				file = new File(file.getAbsolutePath() + ".svg");
			}
			if (file.exists()) {
				result = JOptionPane.showConfirmDialog(w, "File exists, overwrite?", "File exists", OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
				if (result != JOptionPane.OK_OPTION) {
					return;
				}
			}
			try {
				exportHyperbolicImageToSVG(file, coverResolution);
			} catch (Exception e2) {
				JOptionPane.showMessageDialog(w, e2.getMessage(), "Error", ERROR_MESSAGE);
			}
		}		
		if (drawTriangulationChecker == s ||
			drawAxesChecker == s ||
			drawPolygonChecker == s ||
			drawCurvesOnSurface == s
		) {
			updateSurface();
			updateDomainImage();				
		}
	}
	

	public void updateDomainImage() {
		if (domainVisualisationPlugin == null) return;
		domainVisualisationPlugin.updateVisualization();
		HalfedgeInterface domInterface = domainVisualisationPlugin.getDomainInterface();
		domInterface.removeTemporaryGeometry(domainRoot);
		if (genus > 1) {
			// texture
			domInterface.addTemporaryGeometry(domainRoot);
			boolean showCircle = false;
			showCircle |= vis.getSelectedHyperbolicModel() == HyperbolicModel.Poincaré;
			showCircle |= vis.getSelectedHyperbolicModel() == HyperbolicModel.Klein;
			unitCircle.setVisible(showCircle);
			Image img = createDomainImage(coverResolution);
			ImageData imgData = new ImageData(img);
			Texture2D tex2d = TextureUtility.createTexture(universalCoverAppearance, POLYGON_SHADER, imgData);
			tex2d.setTextureMatrix(polygonTextureMatrix);
			System.out.println("DiscreteConformalPlugin.updateDomainImage()");
		}
		
		// add curves
		boolean drawCurves = drawCurvesOnSurface.isSelected();
		hif.removeTemporaryGeometry(polygonCurvesRoot);
		hif.removeTemporaryGeometry(axesCurvesRoot);
		if (genus > 1 && drawCurves) {
			AdapterSet aSet = hif.getActiveAdapters();
			int depth = coverRecursionModel.getNumber().intValue();
			Domain domain = (Domain)domainCombo.getSelectedItem();
			boolean drawPolygon = drawPolygonChecker.isSelected();
			boolean drawAxes = drawAxesChecker.isSelected();
			Color polygonColor = polygonColorButton.getColor();
			Color axesColor = axesColorButton.getColor();
			ConverterHeds2JR converter = new ConverterHeds2JR();
			if (drawPolygon) {
				CoHDS curves = null;
				switch (domain) {
					case Cut:
						curves = SurfaceCurveUtility.createSurfaceCurves(cuttedPolygon, surface, aSet, depth, true, false);
						break;
					case Minimal:
						curves = SurfaceCurveUtility.createSurfaceCurves(minimalPolygon, surface, aSet, depth, true, false);
						break;
					case Canonical:
						curves = SurfaceCurveUtility.createSurfaceCurves(canonicalPolygon, surface, aSet, depth, true, false);
						break;
					case Opposite:
						curves = SurfaceCurveUtility.createSurfaceCurves(oppositePolygon, surface, aSet, depth, true, false);
						break;				
				}
				IndexedFaceSet curvesGeom = converter.heds2ifs(curves, aSet);
				polygonCurvesAppearance.setAttribute(LINE_SHADER + "." + DIFFUSE_COLOR, polygonColor);
				polygonCurvesRoot.setGeometry(curvesGeom);
				hif.addTemporaryGeometry(polygonCurvesRoot);
			}
			if (drawAxes) {
				CoHDS axesCurves = null;
				switch (domain) {
					case Cut:
						axesCurves = SurfaceCurveUtility.createSurfaceCurves(cuttedPolygon, surface, aSet, depth, false, true);
						break;
					case Minimal:
						axesCurves = SurfaceCurveUtility.createSurfaceCurves(minimalPolygon, surface, aSet, depth, false, true);
						break;
					case Canonical:
						axesCurves = SurfaceCurveUtility.createSurfaceCurves(canonicalPolygon, surface, aSet, depth, false, true);
						break;
					case Opposite:
						axesCurves = SurfaceCurveUtility.createSurfaceCurves(oppositePolygon, surface, aSet, depth, false, true);
						break;				
				}
				IndexedFaceSet axesGeom = converter.heds2ifs(axesCurves, aSet);
				axesCurvesAppearance.setAttribute(LINE_SHADER + "." + DIFFUSE_COLOR, axesColor);
				axesCurvesRoot.setGeometry(axesGeom);
				hif.addTemporaryGeometry(axesCurvesRoot);
			}
		}
	}
	
	
	public CoHDS getLoaderGeometry() {
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
	
	
	public void updateSurface() {
		if (surface == null) {
			return;
		}
		HyperbolicModel model = genus <= 1 ? Klein : Poincaré;
		vis.setHyperbolicModel(model);
		hif.addLayerAdapter(metricErrorAdapter, false);
		vis.setInterpolation(InterpolationMethod.Incircle);
		hif.set(surface);
	}
	

	public BufferedImage createDomainImage(int res) {
		BufferedImage image = new BufferedImage(res, res, TYPE_INT_ARGB);
		Graphics2D gImage = image.createGraphics();
		drawDomainImage(gImage, res);
		activeDomainImage = image;
		return image;
	}
	

	public void drawDomainImage(Graphics2D g2d, int res) {
		int depth = coverRecursionModel.getNumber().intValue();
		HyperbolicModel model = vis.getSelectedHyperbolicModel();
		VisualizationUtility.drawDomainBackground(g2d, res, model);
		if (drawTriangulationChecker.isSelected()) {
			if (surface == null) {
				throw new RuntimeException("No surface available");
			}
			Color triangulationColor = triangulationColorButton.getColor();
			drawTriangulation(surface, model, g2d, res, triangulationColor);
		}
		Color polygonColor = polygonColorButton.getColor();
		Color axesColor = axesColorButton.getColor();
		boolean drawPolygon = drawPolygonChecker.isSelected();
		boolean drawAxes = drawAxesChecker.isSelected();
		Domain domain = (Domain)domainCombo.getSelectedItem();
		switch (domain) {
			case Cut:
				if (cuttedPolygon == null) {
					throw new RuntimeException("No fundamental polygon available");
				}
				drawUniversalCoverImage(cuttedPolygon, drawPolygon, drawAxes, depth, model, g2d, res, polygonColor, axesColor);
				break;
			case Minimal:
				if (minimalPolygon == null) {
					throw new RuntimeException("No fundamental polygon available");
				}
				drawUniversalCoverImage(minimalPolygon, drawPolygon, drawAxes, depth, model, g2d, res, polygonColor, axesColor);
				break;
			case Canonical:
				if (canonicalPolygon == null) {
					throw new RuntimeException("No fundamental polygon available");
				}	
				drawUniversalCoverImage(canonicalPolygon, drawPolygon, drawAxes, depth, model, g2d, res, polygonColor, axesColor);
				break;
			case Opposite:
				if (oppositePolygon == null) {
					throw new RuntimeException("No fundamental polygon available");
				}
				drawUniversalCoverImage(oppositePolygon, drawPolygon, drawAxes, depth, model, g2d, res, polygonColor, axesColor);
				break;
		}
	}
	
	
	public void exportHyperbolicImageToPDF(File file, int res) throws Exception {
		FileOutputStream out = new FileOutputStream(file);
		Rectangle pageSize = new Rectangle(-10, -10, res + 20, res + 20);
		Document doc = new Document(pageSize);
		PdfWriter writer = PdfWriter.getInstance(doc, out);
		doc.open();
		PdfContentByte cb = writer.getDirectContent();
		Graphics2D g2 = cb.createGraphics(res, res);
		try {
			drawDomainImage(g2, res);
		} finally {
			g2.dispose();
			doc.close();
			out.close();
		}
	}
	
	public void exportHyperbolicImageToSVG(File file, int res) throws Exception {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		org.w3c.dom.Document doc = builder.newDocument();
		SVGGraphics2D svg = new SVGGraphics2D(doc);
		svg.setSVGCanvasSize(new Dimension(res, res));
		FileWriter writer = new FileWriter(file);
		try {
			drawDomainImage(svg, res);
			svg.stream(writer);
		} finally {
			svg.dispose();
			writer.close();
		}
	}
	
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		hif.addAdapter(positionAdapter, true);
		hif.addAdapter(texturePositionAdapter, true);
		hif.addSelectionListener(this);
		vis = c.getPlugin(ConformalVisualizationPlugin.class);
		contentAppearance = c.getPlugin(ContentAppearance.class);
		domainVisualisationPlugin = c.getPlugin(DomainVisualisationPlugin.class);
		jobQueue = c.getPlugin(JobQueuePlugin.class);
	}

	
	@Override
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "numCones", numConesModel.getNumber().intValue());
		c.storeProperty(getClass(), "conesQuantizationModeIndex", conesQuantizationModeCombo.getSelectedIndex());
		c.storeProperty(getClass(), "boundaryModeIndex", boundaryModeCombo.getSelectedIndex());
		c.storeProperty(getClass(), "boundaryQuantModeIndex", boundaryQuantizationCombo.getSelectedIndex());
		c.storeProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex());
		c.storeProperty(getClass(), "coverRecursion", coverRecursionModel.getNumber());
		c.storeProperty(getClass(), "useProjectiveTexture", useProjectiveTexture.isSelected());
		c.storeProperty(getClass(), "toleranceExponent", toleranceExpModel.getNumber());
		c.storeProperty(getClass(), "maxIterations", maxIterationsModel.getNumber());
		c.storeProperty(getClass(), "boundaryPanelShrinked", boundaryPanel.isShrinked());	
		c.storeProperty(getClass(), "conesPanelShrinked", coneConfigPanel.isShrinked());	
		c.storeProperty(getClass(), "visualizationPanelShrinked", visualizationPanel.isShrinked());	
		c.storeProperty(getClass(), "modelPanelShrinked", modelPanel.isShrinked());	
		c.storeProperty(getClass(), "coneTexQuantPanelShrinked", texQuantizationPanel.isShrinked());
		c.storeProperty(getClass(), "customVertexPanelShrinked", customNodePanel.isShrinked());
		c.storeProperty(getClass(), "domainModeIndex", domainCombo.getSelectedIndex());
	} 
	
 
	@Override
	public void restoreStates(Controller c) throws Exception { 
		super.restoreStates(c);
		numConesModel.setValue(c.getProperty(getClass(), "numCones", numConesModel.getNumber().intValue()));
		numericsCombo.setSelectedIndex(c.getProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex()));
		coverRecursionModel.setValue(c.getProperty(getClass(), "coverRecursion", coverRecursionModel.getNumber()));
		useProjectiveTexture.setSelected(c.getProperty(getClass(), "useProjectiveTexture", useProjectiveTexture.isSelected()));
		toleranceExpModel.setValue(c.getProperty(getClass(), "toleranceExponent", toleranceExpModel.getNumber()));
		maxIterationsModel.setValue(c.getProperty(getClass(), "maxIterations", maxIterationsModel.getNumber()));
		conesQuantizationModeCombo.setSelectedIndex(c.getProperty(getClass(), "conesQuantizationModeIndex", conesQuantizationModeCombo.getSelectedIndex()));
		boundaryModeCombo.setSelectedIndex(c.getProperty(getClass(), "boundaryModeIndex", boundaryModeCombo.getSelectedIndex()));
		boundaryQuantizationCombo.setSelectedIndex(c.getProperty(getClass(), "boundaryQuantModeIndex", boundaryQuantizationCombo.getSelectedIndex()));
		boundaryPanel.setShrinked(c.getProperty(getClass(), "boundaryPanelShrinked", true));
		coneConfigPanel.setShrinked(c.getProperty(getClass(), "conesPanelShrinked", true));
		visualizationPanel.setShrinked(c.getProperty(getClass(), "visualizationPanelShrinked", true));
		modelPanel.setShrinked(c.getProperty(getClass(), "modelPanelShrinked", true));
		texQuantizationPanel.setShrinked(c.getProperty(getClass(), "coneTexQuantPanelShrinked", true));
		customNodePanel.setShrinked(c.getProperty(getClass(), "customVertexPanelShrinked", customNodePanel.isShrinked()));
		domainCombo.setSelectedIndex(c.getProperty(getClass(), "domainModeIndex", domainCombo.getSelectedIndex()));
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
	
	public boolean isRescaleGeometry() {
		return rescaleChecker.isSelected();
	}
	
}
