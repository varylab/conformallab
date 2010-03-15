package de.varylab.discreteconformal;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.experimental.ManagedContentGUI;
import de.jreality.plugin.experimental.WebContentLoader;
import de.jtem.halfedgetools.plugin.HalfedgeDebuggerPlugin;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeToolBar;
import de.jtem.halfedgetools.plugin.algorithm.subdivision.TriangulatePlugin;
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
		viewer.setPropertiesResource(ConformalLab.class, "ConformalLab.jrw");
		viewer.registerPlugin(new ManagedContentGUI());
		viewer.registerPlugin(new SurfaceRemeshingPlugin());
		viewer.registerPlugin(new DiscreteConformalPlugin());
		viewer.registerPlugin(new HalfedgeInterface());
		viewer.registerPlugin(new HalfedgeDebuggerPlugin());
		viewer.registerPlugin(new HalfedgeToolBar());
		viewer.registerPlugin(new TriangulatePlugin());
		viewer.registerPlugin(new WebContentLoader());
		viewer.startup();
	} 

}
