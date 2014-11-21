package de.varylab.discreteconformal.plugin;

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

import de.jreality.math.Matrix;
import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.basic.ViewToolBar;
import de.jreality.plugin.content.ContentAppearance;
import de.jreality.plugin.menu.BackgroundColor;
import de.jreality.plugin.menu.CameraMenu;
import de.jreality.scene.event.AppearanceEvent;
import de.jreality.scene.event.AppearanceListener;
import de.jreality.ui.AppearanceInspector;
import de.jreality.util.SystemProperties;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.Adapter;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeListener;
import de.jtem.halfedgetools.plugin.MarqueeSelectionPlugin;
import de.jtem.halfedgetools.plugin.misc.VertexEditorPlugin;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.halfedgetools.selection.SelectionListener;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.adapter.CoTextureDomainPositionAdapter;
import de.varylab.discreteconformal.plugin.algorithm.MercatorTextureProjection;
import de.varylab.discreteconformal.plugin.algorithm.StereographicTextureProjection;
import de.varylab.discreteconformal.plugin.image.ImageHook;

public class TextureSpaceViewer3D extends ShrinkPanelPlugin {

	private HalfedgeInterface
		mainHif = null,
		visHif = null;
	private ContentAppearance
		mainAppearance = null,
		visAppearance = null;
	private CoTextureDomainPositionAdapter
		domainAdapter = new CoTextureDomainPositionAdapter();
	private Map<HalfedgeLayer, Matrix>
		copyTransformMap = new HashMap<HalfedgeLayer, Matrix>();
	
	private JPanel
		viewerPanel = new JPanel();
	
	private SelectionListener
		visSelectionListener = null,
		mainSelectionListener = null;
	private HalfedgeListener
		mainHalfedgeListener = null;

	private JRViewer 
		domainViewer = null;
	
	private Plugin[] 
			instances = null;
	
	public JRViewer getDomainViewer() {
		return domainViewer;
	}

	public TextureSpaceViewer3D() {
		shrinkPanel.setTitle("Texture Space Viewer 3D");
		setInitialPosition(SHRINKER_TOP);
		shrinkPanel.setLayout(new GridLayout());
		shrinkPanel.add(viewerPanel);
		shrinkPanel.setShrinked(true);
	}

	public TextureSpaceViewer3D(Plugin... P) {
		this();
		this.instances = new Plugin[P.length];
		System.arraycopy(P, 0, instances, 0, instances.length);
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
		visHif.encompassContent();
	}

	public void updateAdapters() {
		for (Adapter<?> a : mainHif.getAdapters()) {
			visHif.addAdapter(a, true);
		}
	}
	
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		mainHif = c.getPlugin(HalfedgeInterface.class);
		viewerPanel.setLayout(new GridLayout());
		// setup viewer, inject viewer system property
		String oldViewerProperty = System.getProperty(SystemProperties.VIEWER); 
		System.setProperty(SystemProperties.VIEWER, SystemProperties.VIEWER_DEFAULT_JOGL);
		domainViewer = new JRViewer();
		if (oldViewerProperty != null) {
			System.setProperty(SystemProperties.VIEWER, oldViewerProperty);
		}
		
		domainViewer.addBasicUI();
		domainViewer.addContentSupport(ContentType.Raw);
		domainViewer.addContentUI();
		domainViewer.getController().setPropertyEngineEnabled(false);
		if(instances != null) {
			for(Plugin p : instances) {
				domainViewer.registerPlugin(p);
			}
		}
		domainViewer.registerPlugin(HalfedgeInterface.class);
		domainViewer.registerPlugin(ContentAppearance.class);
		domainViewer.registerPlugin(MarqueeSelectionPlugin.class);
		domainViewer.registerPlugin(VertexEditorPlugin.class);
		domainViewer.registerPlugin(new StereographicTextureProjection(mainHif));
		domainViewer.registerPlugin(new MercatorTextureProjection(mainHif));
		domainViewer.setShowPanelSlots(false, false, false, false);
		domainViewer.setShowToolBar(true);
		domainViewer.setShowMenuBar(true);
		domainViewer.getController().setRegisterSPIPlugins(false);
		JFrame.setDefaultLookAndFeelDecorated(false);
		JRootPane viewerRoot = domainViewer.startupLocal();
		viewerRoot.setPreferredSize(new Dimension(200, 400));
		viewerPanel.add(viewerRoot);
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
			public void selectionChanged(Selection s, HalfedgeInterface sif) {
				visHif.removeSelectionListener(visSelectionListener);
				try {
					updateAdapters();
					visHif.setSelection(mainHif.getSelection());
				} finally {
					visHif.addSelectionListener(visSelectionListener);
				}
			}
		};
		visSelectionListener = new SelectionListener() {
			@Override
			public void selectionChanged(Selection s, HalfedgeInterface sif) {
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
		mainHalfedgeListener = new HalfedgeListener() {
			@Override
			public void layerRemoved(HalfedgeLayer layer) {
			}
			@Override
			public void layerCreated(HalfedgeLayer layer) {
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
		};
		mainHif.addHalfedgeListener(mainHalfedgeListener);
		domainViewer.getPlugin(BackgroundColor.class).setColor("UI Background");
		domainViewer.getPlugin(CameraMenu.class).setZoomEnabled(true);
		ViewToolBar toolBar = domainViewer.getPlugin(ViewToolBar.class);
		toolBar.addSeparator(TextureSpaceViewer3D.class, 9999.0);
		toolBar.addAction(TextureSpaceViewer3D.class, 10000.0, new UpdateAction());
		toolBar.addAction(TextureSpaceViewer3D.class, 10000.0, new MarkBoundariesAction());
		updateVisualization();
	}
	
	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = super.getPluginInfo();
		info.icon = ImageHook.getIcon("paintbrush.png");
		return info;
	}
	
	
	boolean transferInProgress = false;
	
	public void transferAppearance(AppearanceInspector src, AppearanceInspector dst) {
		if (transferInProgress) return;
		try {
			transferInProgress = true;
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
				Selection s = new Selection();
				List<CoEdge> bEdges = HalfEdgeUtils.boundaryEdges(l.get(new CoHDS()));
				s.addAll(bEdges);
				l.setSelection(s);
			}
		}
		
	}
	
}
