package de.varylab.discreteconformal;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.experimental.ManagedContentGUI;
import de.jreality.plugin.experimental.WebContentLoader;
import de.jtem.halfedgetools.plugin.HalfedgeDebuggerPlugin;
import de.jtem.halfedgetools.plugin.HalfedgeInterfacePlugin;
import de.jtem.halfedgetools.plugin.HalfedgeToolBar;
import de.jtem.halfedgetools.plugin.buildin.CatmullClarkPlugin;
import de.jtem.halfedgetools.plugin.buildin.TriangulatePlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.SurfaceRemeshingPlugin;


public class ConformalLab {

	public static void main(String[] args) {
		JRViewer viewer = new JRViewer();
		viewer.addBasicUI();
		viewer.addContentUI();
		viewer.setShowToolBar(true);
		viewer.setShowPanelSlots(true, true, true, true);
		viewer.addContentSupport(ContentType.CenteredAndScaled);
		viewer.setPropertiesFile("ConformalLab.jrw");
		viewer.registerPlugin(new ManagedContentGUI());
		viewer.registerPlugin(new SurfaceRemeshingPlugin());
		viewer.registerPlugin(new DiscreteConformalPlugin());
		viewer.registerPlugin(new HalfedgeInterfacePlugin<CoVertex, CoEdge, CoFace, CoHDS>());
		viewer.registerPlugin(new HalfedgeDebuggerPlugin<CoVertex, CoEdge, CoFace, CoHDS>());
		viewer.registerPlugin(new HalfedgeToolBar());
		viewer.registerPlugin(new CatmullClarkPlugin<CoVertex, CoEdge, CoFace, CoHDS>());
		viewer.registerPlugin(new TriangulatePlugin<CoVertex, CoEdge, CoFace, CoHDS>());
		viewer.registerPlugin(new WebContentLoader());
		viewer.startup();
	}

}
