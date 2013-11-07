package de.varylab.discreteconformal.plugin;

import static de.jreality.math.Pn.EUCLIDEAN;
import static de.jreality.math.Pn.HYPERBOLIC;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.AbstractAction;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRootPane;

import cern.colt.Arrays;
import de.jreality.math.Matrix;
import de.jreality.math.P2;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.JRViewerUtility;
import de.jreality.plugin.basic.Scene;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.basic.ViewToolBar;
import de.jreality.plugin.content.ContentAppearance;
import de.jreality.plugin.menu.BackgroundColor;
import de.jreality.plugin.menu.CameraMenu;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.SceneGraphPath;
import de.jreality.scene.Transformation;
import de.jreality.scene.event.AppearanceEvent;
import de.jreality.scene.event.AppearanceListener;
import de.jreality.scene.tool.Tool;
import de.jreality.ui.AppearanceInspector;
import de.jreality.util.SceneGraphUtility;
import de.jreality.util.SystemProperties;
import de.jtem.discretegroup.core.DiscreteGroup;
import de.jtem.discretegroup.core.DiscreteGroupConstraint;
import de.jtem.discretegroup.core.DiscreteGroupElement;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.Adapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeSelection;
import de.jtem.halfedgetools.plugin.SelectionListener;
import de.jtem.halfedgetools.plugin.misc.VertexEditorPlugin;
import de.jtem.halfedgetools.plugin.widget.MarqueeWidget;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoTextureDomainPositionAdapter;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class DomainVisualisationPlugin extends ShrinkPanelPlugin {

	private HalfedgeInterface
		mainHif = null,
		visHif = null;
	private Scene
		domainScene = null;
	private ContentAppearance
		mainAppearance = null,
		visAppearance = null;
	private CoTextureDomainPositionAdapter
		domainAdapter = new CoTextureDomainPositionAdapter();
	private DiscreteConformalPlugin
		conformalPlugin = null;
	private ConformalVisualizationPlugin
		conformalVisualizationPlugin = null;
	private Map<HalfedgeLayer, Matrix>
		copyTransformMap = new HashMap<HalfedgeLayer, Matrix>();
	
	private JPanel
		viewerPanel = new JPanel();
	
	private SceneGraphComponent
		copiesComponent = new SceneGraphComponent("Copies"); 
	
	private SelectionListener
		visSelectionListener = null,
		mainSelectionListener = null;

	public DomainVisualisationPlugin() {
		shrinkPanel.setTitle("Parameterization Domain");
		setInitialPosition(SHRINKER_TOP);
		shrinkPanel.setLayout(new GridLayout());
		shrinkPanel.add(viewerPanel);
		shrinkPanel.setShrinked(true);
	}

	private class UpdateAction extends AbstractAction {
		
		private static final long serialVersionUID = 1L;

		public UpdateAction() {
			putValue(NAME, "Update Domain");
		}
		
		@Override
		public void actionPerformed(ActionEvent arg0) {
			updateVisualization();
		}
		
	}
	
	
	public void updateVisualization() {
		if (shrinkPanel.isShrinked()) {
			return;
		}
		for (HalfedgeLayer l : copyTransformMap.keySet()) {
			visHif.removeLayer(l);
		}
		copyTransformMap.clear();
		updateAdapters();
		visHif.set(mainHif.get());
		addCopyTool(visHif.getActiveLayer());
		JRViewerUtility.encompassEuclidean(domainScene);
	}

	public void updateAdapters() {
		for (Adapter<?> a : mainHif.getAdapters()) {
			visHif.addAdapter(a, true);
		}
	}
	
	private void addCopyTool(HalfedgeLayer layer) {
		List<SceneGraphPath> paths = SceneGraphUtility.getPathsToNamedNodes(layer.getLayerRoot(), "Geometry");
		SceneGraphComponent comp = null;
		for (SceneGraphPath path : paths) {
			comp = path.getLastComponent();
			boolean hasTool = false;
			for (Tool tool : comp.getTools()) {
				if (tool instanceof HyperbolicCopyTool) {
					hasTool = true;
					break;
				}
			}
			if (!hasTool) {
				Tool copyDomainTool = new HyperbolicCopyTool(this, layer);
				comp.addTool(copyDomainTool);
			}
		}
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		conformalPlugin = c.getPlugin(DiscreteConformalPlugin.class);
		conformalVisualizationPlugin = c.getPlugin(ConformalVisualizationPlugin.class);
		viewerPanel.setLayout(new GridLayout());
		// setup viewer, inject viewer system property
		String oldViewerProperty = System.getProperty(SystemProperties.VIEWER); 
		System.setProperty(SystemProperties.VIEWER, SystemProperties.VIEWER_DEFAULT_JOGL);
		JRViewer domainViewer = new JRViewer();
		if (oldViewerProperty != null) {
			System.setProperty(SystemProperties.VIEWER, oldViewerProperty);
		}
		
		domainViewer.addBasicUI();
		domainViewer.addContentSupport(ContentType.Raw);
		domainViewer.addContentUI();
		domainViewer.getController().setPropertyEngineEnabled(false);
		domainViewer.registerPlugin(HalfedgeInterface.class);
		domainViewer.registerPlugin(ContentAppearance.class);
		domainViewer.registerPlugin(MarqueeWidget.class);
		domainViewer.registerPlugin(VertexEditorPlugin.class);
		domainViewer.setShowPanelSlots(false, false, false, false);
		domainViewer.setShowToolBar(true);
		domainViewer.setShowMenuBar(true);
		domainViewer.getController().setRegisterSPIPlugins(false);
		JFrame.setDefaultLookAndFeelDecorated(false);
		JRootPane viewerRoot = domainViewer.startupLocal();
		viewerRoot.setPreferredSize(new Dimension(200, 400));
		viewerPanel.add(viewerRoot);
		mainHif = c.getPlugin(HalfedgeInterface.class);
		mainAppearance = c.getPlugin(ContentAppearance.class);
		mainAppearance.getAppearanceInspector().getAppearance().addAppearanceListener(new AppearanceListener() {
			@Override
			public void appearanceChanged(AppearanceEvent ev) {
				AppearanceInspector vi = visAppearance.getAppearanceInspector();
				AppearanceInspector mi = mainAppearance.getAppearanceInspector();
				transferAppearance(mi, vi);
			}
		});
		visHif = domainViewer.getPlugin(HalfedgeInterface.class);
		visHif.addAdapter(domainAdapter, true);
		visAppearance = domainViewer.getPlugin(ContentAppearance.class);
		visAppearance.getAppearanceInspector().getAppearance().addAppearanceListener(new AppearanceListener() {
			@Override
			public void appearanceChanged(AppearanceEvent ev) {
				AppearanceInspector vi = visAppearance.getAppearanceInspector();
				AppearanceInspector mi = mainAppearance.getAppearanceInspector();
				transferAppearance(vi, mi);
			}
		});
		AppearanceInspector vi = visAppearance.getAppearanceInspector();
		AppearanceInspector mi = mainAppearance.getAppearanceInspector();
		transferAppearance(mi, vi);
		mainSelectionListener = new SelectionListener() {
			@Override
			public void selectionChanged(HalfedgeSelection s, HalfedgeInterface sif) {
				visHif.removeSelectionListener(visSelectionListener);
				try {
					visHif.setSelection(mainHif.getSelection());
				} finally {
					visHif.addSelectionListener(visSelectionListener);
				}
			}
		};
		visSelectionListener = new SelectionListener() {
			@Override
			public void selectionChanged(HalfedgeSelection s, HalfedgeInterface sif) {
				mainHif.removeSelectionListener(mainSelectionListener);
				try {
					mainHif.setSelection(visHif.getSelection());
				} finally {
					mainHif.addSelectionListener(mainSelectionListener);
				}
			}
		};
		mainHif.addSelectionListener(mainSelectionListener);
		visHif.addSelectionListener(visSelectionListener);
		domainScene = domainViewer.getPlugin(Scene.class);
		domainViewer.getPlugin(BackgroundColor.class).setColor("UI Background");
		domainViewer.getPlugin(CameraMenu.class).setZoomEnabled(true);
		ViewToolBar toolBar = domainViewer.getPlugin(ViewToolBar.class);
		toolBar.addSeparator(DomainVisualisationPlugin.class, 9999.0);
		toolBar.addAction(DomainVisualisationPlugin.class, 10000.0, new UpdateAction());
		toolBar.addAction(DomainVisualisationPlugin.class, 10000.0, new MarkBoundariesAction());
	}
	
	boolean transferInProgress = false;
	public void transferAppearance(AppearanceInspector src, AppearanceInspector dst) {
		if (transferInProgress) return;
		try {
			transferInProgress = true;
			System.out.println("transfer ");
			// texture apprearance
			dst.setTextureScaleLock(src.isTextureScaleLock());
			dst.setTextureScaleU(src.getTextureScaleU());
			dst.setTextureScaleV(src.getTextureScaleV());
			dst.setTextureTranslationU(src.getTextureTranslationU());
			dst.setTextureTranslationV(src.getTextureTranslationV());
			dst.setTextureRotationAngle(src.getTextureRotationAngle());
			dst.setTextureShearAngle(src.getTextureShearAngle());
			dst.setTextures(src.getTextures());
			dst.setTexture(src.getTexture());
			
			// vertices, edges, and faces
			dst.setShowPoints(src.isShowPoints());
			dst.setShowLines(src.isShowLines());
			dst.setShowFaces(src.isShowFaces());
			dst.setPointRadius(src.getPointRadius());
			dst.setSpheres(src.isSpheres());
			dst.setTubeRadius(src.getTubeRadius());
			dst.setTubes(src.isTubes());
			dst.setPointColor(src.getPointColor());
			dst.setLineColor(src.getLineColor());
			dst.setFaceColor(src.getFaceColor());
			dst.setFacesFlat(src.isFacesFlat());
		} finally {
			transferInProgress = false;
		}
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}
	
	public HalfedgeInterface getDomainInterface() {
		return visHif;
	}

	
	private class MarkBoundariesAction extends AbstractAction {

		private static final long serialVersionUID = 1L;

		public MarkBoundariesAction() {
			putValue(NAME, "Mark Boundaries");
		}
		
		@Override
		public void actionPerformed(ActionEvent e) {
			for (HalfedgeLayer l : visHif.getAllLayers()) {
				HalfedgeSelection s = new HalfedgeSelection();
				List<CoEdge> bEdges = HalfEdgeUtils.boundaryEdges(l.get(new CoHDS()));
				s.addAll(bEdges);
				l.setSelection(s);
			}
		}
		
	}
	
	
	@Position
	private class TransformDomainAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {
		
		private Matrix
			transform = null;
		
		public TransformDomainAdapter(Matrix transform) {
			super(CoVertex.class, null, null, double[].class, true, false);
			this.transform = transform;
		}
		
		@Override
		public double[] getVertexValue(CoVertex v, AdapterSet a) {
			double[] t = transform.multiplyVector(v.T);
			switch (conformalVisualizationPlugin.getSelectedHyperbolicModel()) {
			case Klein:
				return t;
			case Poincaré: 
			default:
				return new double[] {t[0], t[1], 0.0, t[3] + 1};
			case Halfplane:
				return new double[] {t[1], 1, 0.0, t[3] - t[0]};
			}
		}
		
		@Override
		public double getPriority() {
			return 10000.0;
		}
	}
	
	
	public void copyDomainAtEdge(int pickIndex, HalfedgeLayer layer) {
		CoHDS surface = layer.get(new CoHDS());
		CoFace pickFace = surface.getFace(pickIndex);
		CoEdge edge = null;
		for (CoEdge e : HalfEdgeUtils.boundaryEdges(pickFace)) {
			if (e.getRightFace() == null) {
				edge = e;
				break;
			}
		}
		if (edge == null) {
			System.out.println("no boundary face selected for domain copy");
			return;
		}
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = conformalPlugin.getCurrentCutInfo();
		if (cutInfo == null) {
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

		Matrix T = copyTransformMap.get(layer);
		if (T == null) {
			T = new Matrix();
		}
		
		int genus = conformalPlugin.getCurrentGenus();
		int signature = genus > 1 ? HYPERBOLIC : EUCLIDEAN;
		double[] s1 = Pn.normalize(null, edge.getStartVertex().T, signature); 
		double[] t1 = Pn.normalize(null, edge.getTargetVertex().T, signature); 
		double[] s2 = Pn.normalize(null, coEdge.getStartVertex().T, signature); 
		double[] t2 = Pn.normalize(null, coEdge.getTargetVertex().T, signature); 
		
		double dist1 = Pn.distanceBetween(s1, t1, signature);
		double dist2 = Pn.distanceBetween(s2, t2, signature);
		
		assert Math.abs(dist1 - dist2) < 1E-8 : "corresponding edges have different lengths";
		
		double[] a = P2.makeDirectIsometryFromFrames(null, 
			P2.projectP3ToP2(null, s2), 
			P2.projectP3ToP2(null, t2), 
			P2.projectP3ToP2(null, t1), 
			P2.projectP3ToP2(null, s1), 
			signature
		);
		Matrix A = new Matrix(P2.imbedMatrixP2InP3(null, a));
		A = Matrix.times(T, A);
		
		TransformDomainAdapter adapter = new TransformDomainAdapter(A);
		HalfedgeLayer copy = visHif.createLayer("Isometric Copy");
		copy.addAdapter(adapter, true);
		copy.set(surface);
		addCopyTool(copy);
		copyTransformMap.put(copy, A);
	}
	
	public void createCopies(DiscreteGroup G, final int numCopies) {
		G.setDefaultFundamentalDomain(visHif.getActiveLayer().getGeometry());
		DiscreteGroupConstraint constraint = new DiscreteGroupConstraint() {
			@Override
			public void update() {
			}
			@Override
			public void setMaxNumberElements(int arg0) {
			}
			@Override
			public int getMaxNumberElements() {
				return numCopies;
			}
			@Override
			public boolean acceptElement(DiscreteGroupElement s) {
				return true;
			}
		};
		G.setConstraint(constraint);
		G.generateElements();
		visHif.removeTemporaryGeometry(copiesComponent);
		copiesComponent = new SceneGraphComponent("Copies");
		for (DiscreteGroupElement s : G.getElementList()) {
			SceneGraphComponent element = new SceneGraphComponent(s.getWord());
			element.setGeometry(visHif.getActiveLayer().getGeometry());
			Transformation T = new Transformation(s.getArray());
			element.setTransformation(T);
			copiesComponent.addChild(element);
		}
		visHif.addTemporaryGeometry(copiesComponent);
	}
	
	
	public static void main(String[] args) {
		double[] p0 = {0,0,1};
		double[] p1 = {1,0,1};
		double[] q0 = {0,1,1};
		double[] q1 = {1,1,1};
		double[] a = P2.makeDirectIsometryFromFrames(null, p0, p1, q0, q1, Pn.EUCLIDEAN);
		System.out.println("det: " + Rn.determinant(a));
		double[] checkQ0 = Rn.matrixTimesVector(null, a, p0);
		double[] checkQ1 = Rn.matrixTimesVector(null, a, p1);
		System.out.println(Arrays.toString(checkQ0) + " == " + Arrays.toString(q0));
		System.out.println(Arrays.toString(checkQ1) + " == " + Arrays.toString(q1));
	}

}
