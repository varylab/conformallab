package de.varylab.discreteconformal.frontend.controller;

import static de.jreality.geometry.IndexedFaceSetUtility.calculateAndSetNormals;
import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.SMOOTH_SHADING;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static java.awt.Color.WHITE;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.tools.ClickWheelCameraZoomTool;
import de.jreality.ui.viewerapp.ViewerApp;
import de.jreality.util.CameraUtility;
import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.frontend.controller.GeometryController.GeometryChangedListener;
import de.varylab.discreteconformal.heds.CoHDS;

public class UIController implements GeometryChangedListener{

	private SceneGraphComponent
		root = new SceneGraphComponent(),
		meshRoot = new SceneGraphComponent();
	private Appearance
		rootAppearance = new Appearance(),
		meshAppearance = new Appearance();
	private ViewerApp
		viewerApp = new ViewerApp(root);
	private ClickWheelCameraZoomTool
		cameraZoomTool = new ClickWheelCameraZoomTool();
	
	
	public UIController() {
		viewerApp.update();
		viewerApp.setBackgroundColor(WHITE);
		viewerApp.getSceneRoot().addTool(cameraZoomTool);
		viewerApp.setExternalNavigator(true);
		viewerApp.setExternalBeanShell(true);
		cameraZoomTool.setSpeed(1.2);
		ConformalLab.getGeometryController().addChangeListener(this);
		createScene();
	}
	
	
	private void createScene() {
		rootAppearance.setAttribute(SMOOTH_SHADING, true);
		rootAppearance.setAttribute(VERTEX_DRAW, false);
		rootAppearance.setAttribute(EDGE_DRAW, false);
		rootAppearance.setAttribute(FACE_DRAW, true);
		root.setAppearance(rootAppearance);
		
		meshRoot.setAppearance(meshAppearance);
		root.addChild(meshRoot);
	}
	
	
	public void reset() {
		while (root.getChildComponentCount() > 0)
			root.removeChild(root.getChildComponent(0));
	}
	
	private void updateGeometry() {
		IndexedFaceSet ifs = ConformalLab.getGeometryController().getIndexedFaceSet();
		calculateAndSetNormals(ifs);
		meshRoot.setGeometry(ifs);
	}
	
	public void encompass() {
		CameraUtility.encompass(viewerApp.getCurrentViewer());
		viewerApp.getCurrentViewer().render();
	}
	

	public ViewerApp getViewerAppSource() {
		return viewerApp;
	}
	
	public SceneGraphComponent getSceneRoot() {
		return root;
	}

	public SceneGraphComponent getMeshRoot() {
		return meshRoot;
	}
	
	public void geometryChanged(CoHDS heds) {
		updateGeometry();
	}

	
	public Appearance getMeshAppearance(){
		return meshAppearance;
	}
	
}
