package de.varylab.discreteconformal;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.basic.ConsolePlugin;
import de.jreality.plugin.experimental.ManagedContentGUI;
import de.jreality.plugin.experimental.WebContentLoader;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.plugin.HalfedgePluginFactory;
import de.jtem.halfedgetools.plugin.HalfedgeDebuggerPlugin;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeToolBar;
import de.jtem.halfedgetools.plugin.algorithm.subdivision.TriangulatePlugin;
import de.jtem.halfedgetools.plugin.misc.VertexEditorPlugin;
import de.jtem.halfedgetools.plugin.visualizers.NodeIndexVisualizer;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.EllipticModulusEngine;
import de.varylab.discreteconformal.plugin.EllipticRootGenerator;
import de.varylab.discreteconformal.plugin.ThetaVisualizer;


public class ConformalLab {

	public static void main(String[] args) {
		NativePathUtility.set("native");
		JRViewer v = new JRViewer();
		v.addBasicUI();
		v.addContentUI();
		v.setShowToolBar(true);
		v.setShowPanelSlots(true, true, true, true);
		v.addContentSupport(ContentType.CenteredAndScaled);
		v.setPropertiesFile("ConformalLab.jrw");
		v.setPropertiesResource(ConformalLab.class, "ConformalLab.jrw");
		v.registerPlugin(new ManagedContentGUI());
		v.registerPlugin(new DiscreteConformalPlugin());
		v.registerPlugin(new HalfedgeInterface());
		v.registerPlugin(new HalfedgeDebuggerPlugin());
		v.registerPlugin(new HalfedgeToolBar());
		v.registerPlugin(new TriangulatePlugin());
		v.registerPlugin(new WebContentLoader());
		v.registerPlugin(new ThetaVisualizer());
		v.registerPlugin(new EllipticRootGenerator());
		v.registerPlugin(new EllipticModulusEngine());
		v.registerPlugin(new ConsolePlugin());
		v.registerPlugin(new VertexEditorPlugin());
		v.registerPlugin(new NodeIndexVisualizer());
		v.registerPlugins(HalfedgePluginFactory.createGeometryPlugins());
		v.registerPlugins(HalfedgePluginFactory.createSubdivisionPlugins());
		v.registerPlugins(HalfedgePluginFactory.createTopologyPlugins());
		v.startup();
	} 

}
