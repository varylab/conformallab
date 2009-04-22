package de.varylab.discreteconformal.plugin;

import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.varylab.discreteconformal.util.UniformizationUtility.constructFundamentalPolygon;
import static java.awt.GridBagConstraints.RELATIVE;
import static java.awt.GridBagConstraints.REMAINDER;
import static java.lang.Math.PI;
import geom3d.Point;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.InputStream;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingWorker;

import charlesgunn.jreality.tools.TranslateShapeTool;
import de.jreality.geometry.IndexedFaceSetUtility;
import de.jreality.geometry.IndexedLineSetFactory;
import de.jreality.geometry.Primitives;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.plugin.view.ContentAppearance;
import de.jreality.plugin.view.ContentLoader;
import de.jreality.plugin.view.ManagedContent;
import de.jreality.plugin.view.View;
import de.jreality.plugin.view.ManagedContent.ContentAdapter;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.tool.Tool;
import de.jreality.util.Input;
import de.jtem.halfedge.algorithm.triangulation.Triangulator;
import de.jtem.halfedge.jreality.ConverterHeds2JR;
import de.jtem.halfedge.jreality.adapter.Adapter;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.jtem.halfedge.jreality.adapter.TextCoordsAdapter2Ifs;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
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
import de.varylab.discreteconformal.util.UniformizationUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UniformizationUtility.FundamentalPolygon;
import de.varylab.jrworkspace.plugin.Controller;
import de.varylab.jrworkspace.plugin.PluginInfo;
import de.varylab.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.varylab.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin implements ActionListener, PropertyChangeListener {

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
	private int
		genus = -1;
	
	private Triangulator<CoVertex, CoEdge, CoFace>
		triangulator = new Triangulator<CoVertex, CoEdge, CoFace>();

	private MarkedEdgesAdapter
		cutColorAdapter = new MarkedEdgesAdapter();
	private PointAdapter
		pointAdapter = new PointAdapter();
	
	private Tool
		translateShapeTool = new TranslateShapeTool(),
		hyperbolicCopyTool = new HyperbolicCopyTool(this);
	private SceneGraphComponent
		surfaceRoot = new SceneGraphComponent("Surface Geometry"),
		auxGeometry = new SceneGraphComponent("Aux Geometry"),
		copiedGeometry = new SceneGraphComponent("Click Copy Geometry"),
		fundamentalPolygonRoot = new SceneGraphComponent("Fundamental Polygon"), 
		unitCircle = new SceneGraphComponent("Hyperbolic Boundary");

	// user interface section ------------
	private JButton
		unwrapBtn = new JButton("Unwrap");
	private JPanel
		coneConfigPanel = new JPanel(),
		visualizationPanel = new JPanel();
	private SpinnerNumberModel
		numConesModel = new SpinnerNumberModel(0, 0, 100, 1);
	private JSpinner
		numConesSpinner = new JSpinner(numConesModel);
	private JCheckBox
		quantizeChecker = new JCheckBox("Quantize Cone Angles"),
		showUnwrapped = new JCheckBox("Show Unwrapped Geometry"),
		showUnitCircle = new JCheckBox("Show Unit Cirlce");
	private JRadioButton
		kleinButton = new JRadioButton("Klein"),
		poincareButton = new JRadioButton("Poincaré");
	private JComboBox
		numericsCombo = new JComboBox(new String[] {"Java/MTJ Numerics", "Petsc/Tao Numerics"});
	
	
	public DiscreteConformalPlugin() {
		createLayout();
		unwrapBtn.addActionListener(this);
		showUnwrapped.addActionListener(this);
		showUnitCircle.addActionListener(this);
		kleinButton.addActionListener(this);
		poincareButton.addActionListener(this);
		
		ButtonGroup modelGroup = new ButtonGroup();
		modelGroup.add(kleinButton);
		modelGroup.add(poincareButton);
		
		Appearance circleApp = new Appearance();
		circleApp.setAttribute(EDGE_DRAW, false);
		circleApp.setAttribute(VERTEX_DRAW, false); 
		circleApp.setAttribute(FACE_DRAW, true);
		unitCircle.setAppearance(circleApp);
		unitCircle.setVisible(false);
		MatrixBuilder.euclidean().rotate(PI / 2, 1, 0, 0).assignTo(unitCircle);
		unitCircle.setGeometry(Primitives.torus(1.0, 0.005, 200, 5));
		auxGeometry.addChild(unitCircle);
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
		
		numericsCombo.setLightWeightPopupEnabled(true);
		numericsCombo.setSelectedIndex(0);
		shrinkPanel.add(numericsCombo, c);
		shrinkPanel.add(unwrapBtn, c);
		
		visualizationPanel.setBorder(BorderFactory.createTitledBorder("Visualization"));
		visualizationPanel.setLayout(new GridBagLayout());
		visualizationPanel.add(showUnitCircle, c);
		visualizationPanel.add(showUnwrapped, c);
		c.gridwidth = RELATIVE;
		visualizationPanel.add(kleinButton, c);
		c.gridwidth = REMAINDER;
		visualizationPanel.add(poincareButton, c);
		shrinkPanel.add(visualizationPanel, c);
	}
	
	
	@Override
	public void propertyChange(PropertyChangeEvent evt) {
		if (SwingWorker.StateValue.DONE == evt.getNewValue()) {
			Unwrap unwrapper = (Unwrap)evt.getSource();
			if (unwrapper.isCancelled()) {
				System.out.println("Unwrap jop cancelled: " + unwrapper.getState());
				return;
			}
			surface = unwrapper.getSurface();
			genus = unwrapper.genus;
			if (unwrapper.genus > 1) {
				cutInfo = unwrapper.cutInfo;
				fundamentalPolygon = constructFundamentalPolygon(cutInfo);
				cutColorAdapter.setContext(cutInfo);
				pointAdapter.setContext(cutInfo);
			} else {
				cutInfo = null;
				fundamentalPolygon = null;
			}
			updateViewer(true);
			updateStates();
		}
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object s = e.getSource();
		if (unwrapBtn == s) {
			CoHDS surface = getLoaderGeometry();
			Unwrap unwrapper = new Unwrap(surface);
			unwrapper.setNumCones(numConesModel.getNumber().intValue());
			unwrapper.setQuantizeCones(quantizeChecker.isSelected());
			unwrapper.setUsePetsc(numericsCombo.getSelectedIndex() == 1);
			unwrapper.addPropertyChangeListener(this);
			unwrapper.execute();
		}
		if (showUnwrapped == s || kleinButton == s || poincareButton == s) {
			copiedGeometry.setGeometry(null);
			updateViewer(showUnwrapped == s);
		}
		updateStates();
	}
	
	/*
	 * Set up scene
	 */
	private void updateStates() {
		managedContent.addContentUnique(getClass(), auxGeometry);
		managedContent.addContentUnique(getClass(), surfaceRoot);
		managedContent.addContentUnique(getClass(), copiedGeometry);
		managedContent.addContentUnique(getClass(), fundamentalPolygonRoot);
		managedContent.addToolUnique(getClass(), hyperbolicCopyTool);
		managedContent.addToolUnique(getClass(), translateShapeTool);
		if (genus > 1) {
			unitCircle.setVisible(showUnitCircle.isSelected() && showUnwrapped.isSelected());
			fundamentalPolygonRoot.setVisible(showUnwrapped.isSelected());
		} else {
			unitCircle.setVisible(false);
			fundamentalPolygonRoot.setVisible(false);
		}
	}
	
	
	private CoHDS getLoaderGeometry() {
		CoHDS surface = new CoHDS();
		surface.setTexCoordinatesValid(false);
		hcp.setContentParseRoot(managedContent.getContextRoot(ContentLoader.class));
		hcp.getHalfedgeContent(surface, new PositionAdapter());
		if (surface.numVertices() == 0) {
			return null;
		}
		triangulator.triangulate(surface);
		surface.normalizeCoordinates();
		return surface;
	}
	
	
	private void updateViewer(final boolean align) {
		if (surface == null) {
			return;
		}
		boolean klein = kleinButton.isSelected();
		Adapter texAdapter = new TexCoordAdapter(false, !klein);
		Adapter posAdapter = null;
		if (showUnwrapped.isSelected()) {
			posAdapter = new PositionTexCoordAdapter(!klein);
		} else {
			posAdapter = new PositionAdapter();
		}
		IndexedFaceSet ifs = null;
		if (genus > 1) {
			ifs = hcp.toIndexedFaceSet(surface, true, posAdapter, texAdapter, cutColorAdapter, pointAdapter);
		} else {
			ifs = hcp.toIndexedFaceSet(surface, true, posAdapter, texAdapter);
		}
		surfaceRoot.setGeometry(ifs);
		surfaceRoot.setVisible(true);
		
		if (genus > 1) {
			updateFundamentalPolygon(100);
		}
		
		managedContent.removeAll(ContentLoader.class);
		if (align) {
			managedContent.alignContent();
		}
	}
	
	

	public void updateFundamentalPolygon(int resolution) {
		if (fundamentalPolygon == null || cutInfo == null) {
			return;
		}
		final IndexedLineSetFactory ilsf = new IndexedLineSetFactory();
		final boolean klein = kleinButton.isSelected(); 
		int n = fundamentalPolygon.getLength();
		if (klein) {
			ilsf.setVertexCount(n);
			ilsf.setEdgeCount(n);		
		} else {
			ilsf.setVertexCount(n * resolution);
			ilsf.setEdgeCount(n * resolution);
		}
		double[][] verts = null;
		int[][] edges = null;
		if (klein) {
			verts = new double[n][];
			edges = new int[n][];
		} else {
			verts = new double[n * resolution][];
			edges = new int[n * resolution][];	
		}
		Point pRoot = cutInfo.cutRoot.getTextureCoord();
		
		double[] root = new double[] {pRoot.x(), pRoot.y(), 0.0, pRoot.z()};
		List<double[]> orbit = fundamentalPolygon.getOrbit(root);
		for (int i = 0; i < n; i++) {
			double[] pos = orbit.get(i);
			double[] posNext = orbit.get((i + 1) % n);
			if (klein) {
				verts[i] = new double[] {pos[0], pos[1], 0.0, pos[3]};
				edges[i] = new int[] {i, (i + 1) % (n)};
			} else {
				for (int j = 0; j < resolution; j++) {
					int index = i * resolution + j;
					double t = j / (resolution - 1.0);
					double[] p = Rn.linearCombination(null, t, posNext, 1-t, pos);
					Pn.normalize(p, p, Pn.HYPERBOLIC);
					verts[index] = new double[] {p[0], p[1], 0.0, p[3] + 1};
					edges[index] = new int[] {index, (index + 1) % (n * resolution)};
				}
			}
		}
		ilsf.setVertexCoordinates(verts);
		ilsf.setEdgeIndices(edges);
		ilsf.update();
		fundamentalPolygonRoot.setGeometry(ilsf.getGeometry());
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
		c.getPlugin(ContentAppearance.class);
		hcp = c.getPlugin(HalfedgeConnectorPlugin.class);
		managedContent = c.getPlugin(ManagedContent.class);
		managedContent.addContentListener(new MyContentListener());
		
		// read default scene
		ReaderOBJ reader = new ReaderOBJ();
		InputStream in = getClass().getResourceAsStream("brezelCoarse.obj");
		Input input = Input.getInput("Default OBJ Object", in);
		SceneGraphComponent brezelOBJ = reader.read(input);
		managedContent.setContent(ContentLoader.class, brezelOBJ);
		managedContent.alignContent();
		super.install(c);  
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
		c.storeProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex());
		c.storeProperty(getClass(), "klein", kleinButton.isSelected()); 
		c.storeProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected());
		c.storeProperty(getClass(), "showUnitCircle", showUnitCircle.isSelected());
	} 
	

	@Override
	public void restoreStates(Controller c) throws Exception { 
		super.restoreStates(c);
		numConesModel.setValue(c.getProperty(getClass(), "numCones", numConesModel.getNumber().intValue()));
		quantizeChecker.setSelected(c.getProperty(getClass(), "quantizeCones", quantizeChecker.isSelected()));
		numericsCombo.setSelectedIndex(c.getProperty(getClass(), "numericsMethod", numericsCombo.getSelectedIndex()));
		kleinButton.setSelected(c.getProperty(getClass(), "klein", kleinButton.isSelected()));
		poincareButton.setSelected(!kleinButton.isSelected());
		showUnwrapped.setSelected(c.getProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected()));
		showUnitCircle.setSelected(c.getProperty(getClass(), "showUnitCircle", showUnitCircle.isSelected()));
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
	
	
	private class MyContentListener extends ContentAdapter {
		
		@Override
		public void contentAdded(Class<?> context, SceneGraphComponent c) {
			if (context != ContentLoader.class) {
				return;
			}
			unitCircle.setVisible(false);
			surfaceRoot.setVisible(false);
			copiedGeometry.setVisible(false);
			fundamentalPolygonRoot.setVisible(false);
		}
		
	}

}
