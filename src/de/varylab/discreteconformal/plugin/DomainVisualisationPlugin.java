package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.adapter.HyperbolicModel.Klein;
import static de.varylab.discreteconformal.plugin.InterpolationMethod.Incircle;

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
import de.jreality.scene.event.AppearanceEvent;
import de.jreality.scene.event.AppearanceListener;
import de.jreality.ui.AppearanceInspector;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.HalfedgeListener;
import de.jtem.halfedgetools.plugin.HalfedgeSelection;
import de.jtem.halfedgetools.plugin.SelectionListener;
import de.jtem.halfedgetools.plugin.widget.MarqueeWidget;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;

public class DomainVisualisationPlugin extends ShrinkPanelPlugin implements ActionListener, HalfedgeListener, AppearanceListener, SelectionListener {

	private HalfedgeInterface
		mainHif = null,
		visHif = null;
	private Scene
		domainScene = null;
	private ConformalVisualizationPlugin
		conformalVisualizationPlugin = null;
	private ContentAppearance
		mainAppearance = null,
		visAppearance = null;
	private CoPositionAdapter
		positionAdapter = new CoPositionAdapter(true);
	private CoTexturePositionAdapter
		texturePositionAdapter = new CoTexturePositionAdapter(Klein, Incircle);
	
	private JRViewer
		domainViewer = new JRViewer();
	private JPanel
		viewerPanel = new JPanel();

	public DomainVisualisationPlugin() {
		setInitialPosition(SHRINKER_TOP);
		shrinkPanel.setLayout(new GridLayout());
		shrinkPanel.add(viewerPanel);
	}

	
	@Override
	public void actionPerformed(ActionEvent e) {
		updateVisualization();
	}
	
	
	public void updateVisualization() {
		HyperbolicModel model = conformalVisualizationPlugin.getSelectedHyperbolicModel();
		texturePositionAdapter.setModel(model);
		visHif.set(mainHif.get());
		JRViewerUtility.encompassEuclidean(domainScene);
	}
	
	private void synchronizeApprearances() {
		AppearanceInspector vi = visAppearance.getAppearanceInspector();
		AppearanceInspector mi = mainAppearance.getAppearanceInspector();
		vi.setTextureScaleLock(mi.isTextureScaleLock());
		vi.setTextureScaleU(mi.getTextureScaleU());
		vi.setTextureScaleV(mi.getTextureScaleV());
		vi.setTextureTranslationU(mi.getTextureTranslationU());
		vi.setTextureTranslationV(mi.getTextureTranslationV());
		vi.setTextureRotationAngle(mi.getTextureRotationAngle());
		vi.setTextureShearAngle(mi.getTextureShearAngle());
		vi.setTexture(mi.getTexture());
	}
	
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		viewerPanel.setLayout(new GridLayout());
		domainViewer.addBasicUI();
		domainViewer.addContentSupport(ContentType.Raw);
		domainViewer.addContentUI();
		domainViewer.setPropertiesFile("ConformalDomain.jrw");
		domainViewer.setPropertiesResource(DomainVisualisationPlugin.class, "ConformalDomain.jrw");
		domainViewer.registerPlugin(HalfedgeInterface.class);
		domainViewer.registerPlugin(ContentAppearance.class);
		domainViewer.registerPlugin(MarqueeWidget.class);
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
		visHif.addAdapter(positionAdapter, true);
		visHif.addAdapter(texturePositionAdapter, true);
		conformalVisualizationPlugin = c.getPlugin(ConformalVisualizationPlugin.class);
		mainHif.addSelectionListener(this);
		visHif.addSelectionListener(this);
		domainScene = domainViewer.getPlugin(Scene.class);
	}
	
	@Override
	public void appearanceChanged(AppearanceEvent ev) {
		synchronizeApprearances();
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
		if (ignoreSelectionChanged) return;
		ignoreSelectionChanged = true;
			try {
			if (sif == mainHif) {
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
