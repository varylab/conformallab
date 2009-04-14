package de.varylab.discreteconformal.plugin;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.varylab.discreteconformal.heds.util.CuttingUtility.cutManifoldToDisk;
import static de.varylab.discreteconformal.heds.util.UniformizationUtility.checkLengthsAndAngles;
import static de.varylab.discreteconformal.heds.util.UniformizationUtility.getLengthMap;
import static de.varylab.discreteconformal.heds.util.UniformizationUtility.toFundamentalPolygon;
import static java.awt.GridBagConstraints.RELATIVE;
import static java.awt.GridBagConstraints.REMAINDER;
import static java.lang.Double.MAX_VALUE;
import static java.lang.Math.PI;
import static javax.swing.JOptionPane.ERROR_MESSAGE;
import geom3d.Point;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.InputStream;
import java.util.Arrays;

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
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import no.uib.cipr.matrix.Vector;
import de.jreality.geometry.Primitives;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.plugin.view.AlignedContent;
import de.jreality.plugin.view.ContentAppearance;
import de.jreality.plugin.view.View;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.Appearance;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedge.algorithm.triangulation.Triangulator;
import de.jtem.halfedge.jreality.adapter.Adapter;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.jtem.halfedge.plugin.HalfedgeDebuggerPlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.adapter.PositionTexCoordAdapter;
import de.varylab.discreteconformal.heds.adapter.TexCoordAdapter;
import de.varylab.discreteconformal.heds.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.heds.util.UniformizationUtility.UAdapter;
import de.varylab.discreteconformal.plugin.adapter.HyperbolicLengthWeightAdapter;
import de.varylab.discreteconformal.plugin.adapter.MarkedEdgesAdapter;
import de.varylab.discreteconformal.plugin.adapter.PointAdapter;
import de.varylab.discreteconformal.unwrapper.CDiskUnwrapper;
import de.varylab.discreteconformal.unwrapper.CDiskUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.CHyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapper;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.CUnwrapper;
import de.varylab.jrworkspace.plugin.Controller;
import de.varylab.jrworkspace.plugin.PluginInfo;
import de.varylab.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.varylab.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin implements ActionListener, ChangeListener {

	private AlignedContent
		content = null;
//	private TranslateShapeTool
//		translateShapeTool = new TranslateShapeTool();
	private HyperbolicCopyTool
		hyperbolicCopyTool = new HyperbolicCopyTool(this);
	public static HalfedgeDebuggerPlugin<CoVertex, CoEdge, CoFace>
		halfedgeDebugger = null;
	private HalfedgeConnectorPlugin
		hcp = null;
	private Triangulator<CoVertex, CoEdge, CoFace>
		triangulator = new Triangulator<CoVertex, CoEdge, CoFace>();
	private CoHDS
		activeGeometry = null,
		unwrappedGeometry = null;
	private CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	private MarkedEdgesAdapter
		cutColorAdapter = new MarkedEdgesAdapter();
	private PointAdapter
		pointAdapter = new PointAdapter();
	private SceneGraphComponent
		auxGeometry = new SceneGraphComponent(),
		unitCircle = new SceneGraphComponent();

	private JButton
		unwrapBtn = new JButton("Unwrap"),
		reduceBtn = new JButton("To Fundamental Polygon");
	private JPanel
		coneConfigPanel = new JPanel(),
		geometryPanel = new JPanel(),
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
		euclideanButton = new JRadioButton("Eucliean"),
		hyperbolicButton = new JRadioButton("Hyperbolic"),
		kleinButton = new JRadioButton("Klein"),
		poincareButton = new JRadioButton("Poincaré");
	private JComboBox
		numericsCombo = new JComboBox(new String[] {"Java/MTJ Numerics", "Petsc/Tao Numerics"});

	private int
		numCones = 0; 
	private boolean
		quantizeCones = true,
		usePetsc = false,
		hyperbolic = true,
		klein = true;
	
	
	public DiscreteConformalPlugin() {
		createLayout();
		unwrapBtn.addActionListener(this);
		numConesSpinner.addChangeListener(this);
		quantizeChecker.addActionListener(this);
		numericsCombo.addActionListener(this);
		euclideanButton.addActionListener(this);
		hyperbolicButton.addActionListener(this);
		showUnwrapped.addActionListener(this);
		showUnitCircle.addActionListener(this);
		kleinButton.addActionListener(this);
		poincareButton.addActionListener(this);
		reduceBtn.addActionListener(this);
		
		ButtonGroup geometryGroup = new ButtonGroup();
		geometryGroup.add(euclideanButton);
		geometryGroup.add(hyperbolicButton);
		
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
		
		geometryPanel.setBorder(BorderFactory.createTitledBorder("Geometry"));
		geometryPanel.setLayout(new GridLayout(1, 2));
		shrinkPanel.add(geometryPanel, c);
		geometryPanel.add(euclideanButton);
		geometryPanel.add(hyperbolicButton);
		
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
		
		shrinkPanel.add(reduceBtn, c);
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object s = e.getSource();
		if (unwrapBtn == s) {
			getGeometry();
			new Thread("Reduce To Fundamental Thread") {
				@Override
				public void run() {
					try {
						unwrapActiveGeometry();
					} catch (Exception e1) {
						Window w = SwingUtilities
								.getWindowAncestor(shrinkPanel);
						JOptionPane.showMessageDialog(w, e1.getMessage(),
								"Unwrap Error", ERROR_MESSAGE);
						e1.printStackTrace();
						return;
					}
					updateViewer();							
				}
			}.start();
		}
		if (quantizeChecker == s) {
			quantizeCones = quantizeChecker.isSelected();
		}
		if (numericsCombo == s) {
			usePetsc = numericsCombo.getSelectedIndex() != 0;
		}
		if (euclideanButton == s) {
			hyperbolic = !euclideanButton.isSelected();
		}
		if (hyperbolicButton == s) {
			hyperbolic = hyperbolicButton.isSelected();
		}
		if (showUnwrapped == s || kleinButton == s || poincareButton == s) {
			klein = kleinButton.isSelected();
//			if (klein) {
////				if (!content.getContent().getTools().contains(translateShapeTool)) {
////					content.getContent().addTool(translateShapeTool);
////				}
//				if (!content.getContent().getTools().contains(hyperbolicCopyTool)) {
//					
//				}
////				content.getContent().setTransformation(new Transformation());
////				content.getContent().getAppearance().setAttribute("metric", -1);
//			} else {
////				if (content.getContent().getTools().contains(translateShapeTool)) {
////					content.getContent().removeTool(translateShapeTool);					
////				}
//				if (content.getContent().getTools().contains(hyperbolicCopyTool)) {
//					content.getContent().removeTool(hyperbolicCopyTool);
//				}
////				content.getContent().setTransformation(new Transformation());
////				content.getContent().getAppearance().setAttribute("metric", 0);
//			}
			updateViewer();
		}
		if (showUnitCircle == s) {
			unitCircle.setVisible(showUnitCircle.isSelected());
		}
		if (reduceBtn == s) {
			getGeometry();
			new Thread("Reduce To Fundamental Thread") {
				@Override
				public void run() {
					try {
						reduceToFundamentalPolygon();
					} catch (Exception e1) {
						e1.printStackTrace(); 
						return;
					}
					updateViewer();
				}
			}.start();
		}
	}
	
	private void getGeometry() {
		activeGeometry = new CoHDS();
		unwrappedGeometry = null;
		activeGeometry.setTexCoordinatesValid(false);
		hcp.getHalfedgeContent(activeGeometry, new PositionAdapter());
		triangulator.triangulate(activeGeometry);
		double[][] bounds = new double[][]{{MAX_VALUE, -MAX_VALUE},{MAX_VALUE, -MAX_VALUE},{MAX_VALUE, -MAX_VALUE}};
		for (CoVertex v : activeGeometry.getVertices()) {
			Point p = v.getPosition();
			bounds[0][0] = p.x() < bounds[0][0] ? p.x() : bounds[0][0];
			bounds[0][1] = p.x() > bounds[0][1] ? p.x() : bounds[0][1];
			bounds[1][0] = p.y() < bounds[1][0] ? p.y() : bounds[1][0];
			bounds[1][1] = p.y() > bounds[1][1] ? p.y() : bounds[1][1];
			bounds[2][0] = p.z() < bounds[2][0] ? p.z() : bounds[2][0];
			bounds[2][1] = p.z() > bounds[2][1] ? p.z() : bounds[2][1];
		}
		double xExtend = Math.abs(bounds[0][1] - bounds[0][0]);
		double yExtend = Math.abs(bounds[1][1] - bounds[1][0]);
		double zExtend = Math.abs(bounds[2][1] - bounds[2][0]);
		double max = Math.max(xExtend, Math.max(yExtend, zExtend));
		double xCenter = (bounds[0][1] + bounds[0][0]) / 2 / max;
		double yCenter = (bounds[1][1] + bounds[1][0]) / 2 / max;
		double zCenter = (bounds[2][1] + bounds[2][0]) / 2 / max;
		for (CoVertex v : activeGeometry.getVertices()) {
			Point p = v.getPosition();
			p.times(1 / max);
			p.move(new geom3d.Vector(-xCenter, -yCenter, -zCenter));
		}
	} 
	
	
	private void updateViewer() {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				if (unwrappedGeometry == null) {
					return;
				}
				Adapter texAdapter = new TexCoordAdapter(false, !klein);
				Adapter posAdapter = null;
				if (showUnwrapped.isSelected()) {
					posAdapter = new PositionTexCoordAdapter(!klein);
				} else {
					posAdapter = new PositionAdapter();
				}
				hcp.updateHalfedgeContent(unwrappedGeometry, true, posAdapter, texAdapter, cutColorAdapter, pointAdapter);
			}
		});
	}
	
	
	private void reduceToFundamentalPolygon() throws Exception {
		unwrappedGeometry = activeGeometry.createCombinatoriallyEquivalentCopy(new CoHDS());
		halfedgeDebugger.setData(unwrappedGeometry);
		for (CoVertex v : activeGeometry.getVertices()) {
			unwrappedGeometry.getVertex(v.getIndex()).setPosition(v.getPosition());
		}
		Vector u = null;
		if (usePetsc) {
			CHyperbolicUnwrapperPETSc unwrapper = new CHyperbolicUnwrapperPETSc();
			u = unwrapper.getConformalFactors(unwrappedGeometry);
		} else {
			CHyperbolicUnwrapper unwrapper = new CHyperbolicUnwrapper();
			u = unwrapper.getConformalFactors(unwrappedGeometry);
		}
		UAdapter uAdapter = new UAdapter(u);
		CoVertex root = unwrappedGeometry.getVertex(0);
		toFundamentalPolygon(unwrappedGeometry, root, uAdapter);
		checkLengthsAndAngles(unwrappedGeometry, getLengthMap(unwrappedGeometry, uAdapter));
		HyperbolicLengthWeightAdapter wa = new HyperbolicLengthWeightAdapter(u);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = cutManifoldToDisk(unwrappedGeometry, root, wa);
		checkLengthsAndAngles(unwrappedGeometry, getLengthMap(unwrappedGeometry, uAdapter));
		CHyperbolicLayout.doLayout(unwrappedGeometry, u);
		cutColorAdapter.setContext(cutInfo);
	}
	
	
	
	private void unwrapActiveGeometry() throws Exception {
		unwrappedGeometry = activeGeometry.createCombinatoriallyEquivalentCopy(new CoHDS());
		halfedgeDebugger.setData(unwrappedGeometry);
		for (CoVertex v : activeGeometry.getVertices()) {
			unwrappedGeometry.getVertex(v.getIndex()).setPosition(v.getPosition());
		}
		if (hyperbolic) {
			unwrappedGeometry.prepareInvariantDataHyperbolic();
			Vector u = null;
			if (usePetsc) {
				CHyperbolicUnwrapperPETSc unwrapper = new CHyperbolicUnwrapperPETSc();
				u = unwrapper.getConformalFactors(unwrappedGeometry);
			} else {
				CHyperbolicUnwrapper unwrapper = new CHyperbolicUnwrapper();
				u = unwrapper.getConformalFactors(unwrappedGeometry);
			}
			HyperbolicLengthWeightAdapter hypWa = new HyperbolicLengthWeightAdapter(u);
//			EuclideanLengthWeightAdapter eucWa = new EuclideanLengthWeightAdapter();
			CoVertex root = unwrappedGeometry.getVertex(getMinUIndex(u));
			cutInfo = cutManifoldToDisk(unwrappedGeometry, root, hypWa);
			
			CHyperbolicLayout.doLayout(unwrappedGeometry, u);
			cutColorAdapter.setContext(cutInfo);
			pointAdapter.setContext(cutInfo);
		} else {
			CUnwrapper unwrapper = null;
			if (usePetsc) {
				unwrapper = new CDiskUnwrapperPETSc(numCones, quantizeCones);
			} else {
				unwrapper = new CDiskUnwrapper(numCones, quantizeCones);
			}
			unwrapper.unwrap(unwrappedGeometry);
		}
	}
	
	
	
	public void copyAtEdge(int edgeIndex) {
		if (unwrappedGeometry == null) {
			return;
		}
		CoEdge edge = null;
		int i = 0;
		for (CoEdge e : unwrappedGeometry.getEdges()) {
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
		
		Point s1 = edge.getTargetVertex().getTextureCoord();
		Point s2 = edge.getStartVertex().getTextureCoord();
		Point t1 = coEdge.getTargetVertex().getTextureCoord();
		Point t2 = coEdge.getStartVertex().getTextureCoord();
		System.out.println("s1: " + s1);
		System.out.println("s2: " + s2);
		System.out.println("t1: " + t1);
		System.out.println("t2: " + t2);
		double ds = Pn.distanceBetween(s1.get(), s2.get(), HYPERBOLIC);
		double dt = Pn.distanceBetween(t1.get(), t2.get(), HYPERBOLIC);
		System.out.println("dist: " + ds + ", " + dt);
		
		geom3d.Vector ns = new geom3d.Vector(s1).cross(s2);
		geom3d.Vector nt = new geom3d.Vector(t1).cross(t2);
		
		double s1s1 = Pn.innerProduct(s1.get(), s1.get(), HYPERBOLIC);
		double t1t1 = Pn.innerProduct(t1.get(), t1.get(), HYPERBOLIC);
		double s1s2 = Pn.innerProduct(s1.get(), s2.get(), HYPERBOLIC);
		double t1t2 = Pn.innerProduct(t1.get(), t2.get(), HYPERBOLIC);
		geom3d.Vector ws = new geom3d.Vector(s2).add(new Point(s1).times(s1s2 / s1s1).times(-1));
		geom3d.Vector wt = new geom3d.Vector(t2).add(new Point(t1).times(t1t2 / t1t1).times(-1));
		System.out.println("Orthogonal check ws: " + Pn.innerProduct(ws.get(), s1.get(), HYPERBOLIC));
		System.out.println("Orthogonal check wt: " + Pn.innerProduct(wt.get(), t1.get(), HYPERBOLIC));
		
		double nss1 = Pn.innerProduct(ns.get(), s1.get(), HYPERBOLIC);
		double ntt1 = Pn.innerProduct(nt.get(), t1.get(), HYPERBOLIC);
		ns.add(new Point(s1).times(nss1 / s1s1).times(-1));
		nt.add(new Point(t1).times(ntt1 / t1t1).times(-1));
		System.out.println("Orthogonal check ns1: " + Pn.innerProduct(ns.get(), s1.get(), HYPERBOLIC));
		System.out.println("Orthogonal check nt1: " + Pn.innerProduct(nt.get(), t1.get(), HYPERBOLIC));
		
		double nsws = Pn.innerProduct(ns.get(), ws.get(), HYPERBOLIC);
		double ntwt = Pn.innerProduct(nt.get(), wt.get(), HYPERBOLIC);	
		double wsws = Pn.innerProduct(ws.get(), ws.get(), HYPERBOLIC);
		double wtwt = Pn.innerProduct(wt.get(), wt.get(), HYPERBOLIC);
		ns.add(new Point(ws).times(nsws / wsws).times(-1));
		nt.add(new Point(wt).times(ntwt / wtwt).times(-1));
		Pn.normalize(ns.get(), ns.get(), HYPERBOLIC);
		Pn.normalize(nt.get(), nt.get(), HYPERBOLIC);
		
		System.out.println("ns: " + ns);
		System.out.println("nt: " + nt);
		System.out.println("Orthogonal check ns1: " + Pn.innerProduct(ns.get(), s1.get(), HYPERBOLIC));
		System.out.println("Orthogonal check ns2: " + Pn.innerProduct(ns.get(), s2.get(), HYPERBOLIC));
		System.out.println("Orthogonal check nt1: " + Pn.innerProduct(nt.get(), t1.get(), HYPERBOLIC));
		System.out.println("Orthogonal check nt2: " + Pn.innerProduct(nt.get(), t2.get(), HYPERBOLIC));
		double[] sa1 = s1.get();
		double[] sa2 = s2.get();
		double[] nsa = ns.get();
		double[] ta1 = t1.get();
		double[] ta2 = t2.get();
		double[] nta = nt.get();
		Matrix E = new Matrix(
			1,		0,		0,		0,
			0,		1, 		0,		0,
			0,		0,		-1,		0,
			0,		0, 		0, 		1
		);
		Matrix S = new Matrix(
			sa1[0], sa2[0], nsa[0], 0,
			sa1[1], sa2[1], nsa[1], 0,
			sa1[2], sa2[2], nsa[2], 0,
			0,		0, 		0, 		1
		);
		Matrix T = new Matrix(
			ta1[0], ta2[0], nta[0], 0,
			ta1[1], ta2[1], nta[1], 0,
			ta1[2], ta2[2], nta[2], 0,
			0,		0, 		0, 		1
		);
		Matrix A = Matrix.times(T, S.getInverse());
		System.out.println("A*s1: " + Arrays.toString(A.multiplyVector(sa1)));
		System.out.println("A*s2: " + Arrays.toString(A.multiplyVector(sa2)));
		System.out.println("A*ns: " + Arrays.toString(A.multiplyVector(nsa)));
		System.out.println(A);
		
		System.out.println("t1 * t1: " + Pn.innerProduct(ta1, ta1, HYPERBOLIC));
		System.out.println("t2 * t2: " + Pn.innerProduct(ta2, ta2, HYPERBOLIC));
		double[] At1 = A.multiplyVector(ta1);
		double[] At2 = A.multiplyVector(ta2);
		System.out.println("At1 * At1: " + Pn.innerProduct(At1, At1, HYPERBOLIC));
		System.out.println("At2 * At2: " + Pn.innerProduct(At2, At2, HYPERBOLIC));
		
		System.out.println(Matrix.times(A.getTranspose(), Matrix.times(E, A)));
	}
	
	
	
	
	private int getMinUIndex(Vector u) {
		int index = 0;
		double iVal = u.get(0);
		for (int i = 1; i < u.size(); i++) {
			double val = u.get(i);
			if (iVal < val) {
				index = i;
				iVal = i;
			}
		}
		return index;
	}
	
	
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (numConesSpinner == e.getSource()) {
			numCones = numConesModel.getNumber().intValue();
		}
	}
	
	
	@SuppressWarnings("unchecked")
	@Override
	public void install(Controller c) throws Exception {
		hcp = c.getPlugin(HalfedgeConnectorPlugin.class);
		halfedgeDebugger = c.getPlugin(HalfedgeDebuggerPlugin.class);
		content = c.getPlugin(AlignedContent.class);
		content.getScalingComponent().addChild(auxGeometry);
		content.getScalingComponent().addTool(hyperbolicCopyTool);
		c.getPlugin(ContentAppearance.class);
		ReaderOBJ reader = new ReaderOBJ();
		InputStream in = getClass().getResourceAsStream("brezelCoarse.obj");
		Input input = Input.getInput("Default OBJ Object", in);
		SceneGraphComponent brezelOBJ = reader.read(input);
		content.setContent(brezelOBJ);
		super.install(c);  
	}
	
	@Override
	public void uninstall(Controller c) throws Exception {
		super.uninstall(c);
		if (content.getScalingComponent().isDirectAncestor(auxGeometry)) {
			content.getScalingComponent().removeChild(auxGeometry);
		} 
		if (content.getScalingComponent().getTools().contains(hyperbolicCopyTool)) {
			content.getScalingComponent().removeTool(hyperbolicCopyTool);
		}
	}
	
	
	@Override
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "numCones", numCones);
		c.storeProperty(getClass(), "quantizeCones", quantizeCones);
		c.storeProperty(getClass(), "hyperbolic", hyperbolic);
		c.storeProperty(getClass(), "usePetsc", usePetsc);
		c.storeProperty(getClass(), "klein", klein); 
		c.storeProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected());
		c.storeProperty(getClass(), "showUnitCircle", showUnitCircle.isSelected());
	} 
	

	@Override
	public void restoreStates(Controller c) throws Exception { 
		super.restoreStates(c);
		numConesModel.setValue(c.getProperty(getClass(), "numCones", numCones));
		quantizeChecker.setSelected(c.getProperty(getClass(), "quantizeCones", quantizeCones));
		hyperbolic = c.getProperty(getClass(), "hyperbolic", hyperbolic);
		usePetsc = c.getProperty(getClass(), "usePetsc", usePetsc);
		klein = c.getProperty(getClass(), "klein", klein);
		showUnwrapped.setSelected(c.getProperty(getClass(), "showUnwrapped", showUnwrapped.isSelected()));
		showUnitCircle.setSelected(c.getProperty(getClass(), "showUnitCircle", showUnitCircle.isSelected()));
		
		numericsCombo.setSelectedIndex(usePetsc ? 1 : 0);
		quantizeCones = quantizeChecker.isSelected();
		numCones = numConesModel.getNumber().intValue();
		euclideanButton.setSelected(!hyperbolic);
		hyperbolicButton.setSelected(hyperbolic);
		kleinButton.setSelected(klein);
		poincareButton.setSelected(!klein);
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

}
