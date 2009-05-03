package de.varylab.discreteconformal;

import de.jreality.plugin.ContentViewer;
import de.jreality.plugin.view.ContentGUI;
import de.jreality.plugin.view.ContentLoader;
import de.jreality.plugin.view.InfoOverlayPlugin;
import de.jreality.plugin.view.ViewerKeyListenerPlugin;
import de.jreality.plugin.view.ZoomTool;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.jtem.halfedge.plugin.HalfedgeDebuggerPlugin;
import de.jtem.halfedge.plugin.HalfedgeToolBar;
import de.jtem.halfedge.plugin.buildin.CatmullClarkPlugin;
import de.jtem.halfedge.plugin.buildin.TriangulatePlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;

public class ConformalLab {

	public static void main(String[] args) {
		ContentViewer viewer = new ContentViewer();
		viewer.registerPlugin(new ContentLoader());
		viewer.registerPlugin(new ZoomTool());
		viewer.registerPlugin(new HalfedgeConnectorPlugin());
		viewer.registerPlugin(new HalfedgeDebuggerPlugin<CoVertex, CoEdge, CoFace>());
		viewer.registerPlugin(new HalfedgeToolBar());
		viewer.registerPlugin(new CatmullClarkPlugin());
		viewer.registerPlugin(new TriangulatePlugin());
		viewer.registerPlugin(new DiscreteConformalPlugin());
		viewer.registerPlugin(new ViewerKeyListenerPlugin());
		viewer.registerPlugin(new InfoOverlayPlugin());
		viewer.registerPlugin(new ContentGUI());
		viewer.startup();
	}

}
