package de.varylab.discreteconformal.plugin;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JPanel;
import javax.swing.JRootPane;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.JRViewerUtility;
import de.jreality.plugin.basic.Scene;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.content.ContentAppearance;
import de.jreality.plugin.menu.BackgroundColor;
import de.jreality.scene.event.AppearanceEvent;
import de.jreality.scene.event.AppearanceListener;
import de.jreality.ui.AppearanceInspector;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.Adapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeListener;
import de.jtem.halfedgetools.plugin.HalfedgeSelection;
import de.jtem.halfedgetools.plugin.SelectionListener;
import de.jtem.halfedgetools.plugin.misc.VertexEditorPlugin;
import de.jtem.halfedgetools.plugin.widget.MarqueeWidget;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;

public class DomainVisualisationPlugin extends ShrinkPanelPlugin implements ActionListener, HalfedgeListener, AppearanceListener, SelectionListener {

	private HalfedgeInterface
		mainHif = null,
		visHif = null;
	private Scene
		domainScene = null;
	private ContentAppearance
		mainAppearance = null,
		visAppearance = null;
	private TextureDomainPositionAdapter
		domainAdapter = new TextureDomainPositionAdapter();
	
	private JRViewer
		domainViewer = new JRViewer();
	private JPanel
		viewerPanel = new JPanel();

	public DomainVisualisationPlugin() {
		setInitialPosition(SHRINKER_TOP);
		shrinkPanel.setLayout(new GridLayout());
		shrinkPanel.add(viewerPanel);
		shrinkPanel.setShrinked(true);
	}

	
	@TexturePosition
	@Position
	public static class TextureDomainPositionAdapter extends AbstractAdapter<double[]> {
		
		public TextureDomainPositionAdapter() {
			super(double[].class, true, true);
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Vertex.class.isAssignableFrom(nodeClass);
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> double[] getV(V v, AdapterSet a) {
			a.setPriorityBound(getPriority());
			double[] tp = a.getD(TexturePosition.class, v);
			a.removePriorityBound();
			return tp;
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> void setV(V v, double[] coords, AdapterSet a) {
			a.setPriorityBound(getPriority());
			a.set(TexturePosition.class, v, coords);
			a.removePriorityBound();
		}
		@Override
		public double getPriority() {
			return 1000.0;
		}
		
	}
	
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		updateVisualization();
	}
	
	public void updateVisualization() {
		if (shrinkPanel.isShrinked()) {
			return;
		}
//		HyperbolicModel model = conformalVisualizationPlugin.getSelectedHyperbolicModel();
//		texturePositionAdapter.setModel(model);
		updateAdapters();
		visHif.set(mainHif.get());
		JRViewerUtility.encompassEuclidean(domainScene);
	}

	public void updateAdapters() {
		for (Adapter<?> a : mainHif.getAdapters()) {
			visHif.addAdapter(a, true);
		}
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		viewerPanel.setLayout(new GridLayout());
		domainViewer.addBasicUI();
		domainViewer.addContentSupport(ContentType.Raw);
		domainViewer.addContentUI();
//		domainViewer.setPropertiesFile("ConformalDomain.jrw");
//		domainViewer.setPropertiesResource(DomainVisualisationPlugin.class, "ConformalDomain.jrw");
		domainViewer.getController().setPropertyEngineEnabled(false);
		domainViewer.registerPlugin(HalfedgeInterface.class);
		domainViewer.registerPlugin(ContentAppearance.class);
		domainViewer.registerPlugin(MarqueeWidget.class);
		domainViewer.registerPlugin(VertexEditorPlugin.class);
		domainViewer.setShowPanelSlots(false, false, false, false);
		domainViewer.setShowToolBar(true);
		domainViewer.getController().setSaveOnExit(true);
		domainViewer.getController().setAskBeforeSaveOnExit(false);
		domainViewer.getController().setLoadFromUserPropertyFile(true);
		domainViewer.getController().setRegisterSPIPlugins(false);
		JRootPane viewerRoot = domainViewer.startupLocal();
		viewerRoot.setJMenuBar(null);
		viewerRoot.setPreferredSize(new Dimension(200, 400));
		viewerPanel.add(viewerRoot);
		mainHif = c.getPlugin(HalfedgeInterface.class);
		mainHif.addHalfedgeListener(this);
		mainAppearance = c.getPlugin(ContentAppearance.class);
		mainAppearance.getAppearanceInspector().getAppearance().addAppearanceListener(this);
		visHif = domainViewer.getPlugin(HalfedgeInterface.class);
		visAppearance = domainViewer.getPlugin(ContentAppearance.class);
		visHif.addAdapter(domainAdapter, true);
		mainHif.addSelectionListener(this);
		visHif.addSelectionListener(this);
		domainScene = domainViewer.getPlugin(Scene.class);
		domainViewer.getPlugin(BackgroundColor.class).setColor("Transparent Black");
	}
	
	@Override
	public void appearanceChanged(AppearanceEvent ev) {
		AppearanceInspector vi = visAppearance.getAppearanceInspector();
		AppearanceInspector mi = mainAppearance.getAppearanceInspector();
		
		// texture apprearance
		vi.setTextureScaleLock(mi.isTextureScaleLock());
		vi.setTextureScaleU(mi.getTextureScaleU());
		vi.setTextureScaleV(mi.getTextureScaleV());
		vi.setTextureTranslationU(mi.getTextureTranslationU());
		vi.setTextureTranslationV(mi.getTextureTranslationV());
		vi.setTextureRotationAngle(mi.getTextureRotationAngle());
		vi.setTextureShearAngle(mi.getTextureShearAngle());
		vi.setTexture(mi.getTexture());
		
		// vertices, edges, and faces
		vi.setShowPoints(mi.isShowPoints());
		vi.setShowLines(mi.isShowLines());
		vi.setShowFaces(mi.isShowFaces());
		vi.setPointRadius(mi.getPointRadius());
		vi.setSpheres(mi.isSpheres());
		vi.setTubeRadius(mi.getTubeRadius());
		vi.setTubes(mi.isTubes());
		vi.setPointColor(mi.getPointColor());
		vi.setLineColor(mi.getLineColor());
		vi.setFaceColor(mi.getFaceColor());
		vi.setFacesFlat(mi.isFacesFlat());
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}
	
	public HalfedgeInterface getDomainInterface() {
		return visHif;
	}

	@Override
	public void dataChanged(HalfedgeLayer layer) {
		updateVisualization();
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
	
	private boolean ignoreSelectionChanged = false;
	@Override
	public void selectionChanged(HalfedgeSelection s, HalfedgeInterface sif) {
		if (shrinkPanel.isShrinked()) {
			return;
		}
		if (ignoreSelectionChanged) return;
		ignoreSelectionChanged = true;
			try {
			if (sif == mainHif) {
				updateAdapters();
				visHif.setSelection(s);
			}
			if (sif == visHif) {
				mainHif.setSelection(s);
			}
		} finally {
			ignoreSelectionChanged = false;
		}
	}

}
