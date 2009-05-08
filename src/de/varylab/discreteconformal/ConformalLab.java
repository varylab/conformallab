package de.varylab.discreteconformal;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.view.InfoOverlayPlugin;
import de.jreality.plugin.view.ViewerKeyListenerPlugin;
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
		JRViewer viewer = JRViewer.createViewer();
		viewer.setPropertiesFile("ConformalLab.jrw");
		viewer.registerPlugin(new DiscreteConformalPlugin());
		viewer.registerPlugin(new HalfedgeConnectorPlugin());
		viewer.registerPlugin(new HalfedgeDebuggerPlugin<CoVertex, CoEdge, CoFace>());
		viewer.registerPlugin(new HalfedgeToolBar());
		viewer.registerPlugin(new CatmullClarkPlugin());
		viewer.registerPlugin(new TriangulatePlugin());
		viewer.registerPlugin(new ViewerKeyListenerPlugin());
		viewer.registerPlugin(new InfoOverlayPlugin());
		viewer.startup();
	}

}
