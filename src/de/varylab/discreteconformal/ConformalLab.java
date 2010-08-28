package de.varylab.discreteconformal;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.basic.ConsolePlugin;
import de.jreality.plugin.experimental.WebContentLoader;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.JRHalfedgeViewer;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgePluginFactory;
import de.jtem.halfedgetools.plugin.HalfedgeToolBar;
import de.jtem.halfedgetools.plugin.algorithm.subdivision.TriangulatePlugin;
import de.jtem.halfedgetools.plugin.misc.VertexEditorPlugin;
import de.jtem.halfedgetools.plugin.visualizers.NodeIndexVisualizer;
import de.jtem.jrworkspace.plugin.lnfswitch.LookAndFeelSwitch;
import de.jtem.jrworkspace.plugin.lnfswitch.plugin.CrossPlatformLnF;
import de.jtem.jrworkspace.plugin.lnfswitch.plugin.NimbusLnF;
import de.jtem.jrworkspace.plugin.lnfswitch.plugin.SystemLookAndFeel;
import de.varylab.discreteconformal.plugin.CutAtEdgePlugin;
import de.varylab.discreteconformal.plugin.CutToDiskPlugin;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.EllipticModulusEngine;
import de.varylab.discreteconformal.plugin.EllipticRootGenerator;
import de.varylab.discreteconformal.plugin.FindPathPlugin;
import de.varylab.discreteconformal.plugin.ThetaVisualizer;


public class ConformalLab {

	public static void main(String[] args) {
		JRHalfedgeViewer.initHalfedgeFronted();
		NativePathUtility.set("native");
		JRViewer v = new JRViewer();
		v.addBasicUI();
		v.addContentUI();
		v.setShowToolBar(true);
		v.setShowPanelSlots(true, true, true, true);
		v.addContentSupport(ContentType.Raw);
		v.setPropertiesFile("ConformalLab.jrw");
		v.setPropertiesResource(ConformalLab.class, "ConformalLab.jrw");
		v.getController().setManageLookAndFeel(true);
		v.registerPlugin(new DiscreteConformalPlugin());
		v.registerPlugin(new HalfedgeInterface());
//		v.registerPlugin(new HalfedgeDebuggerPlugin());
		v.registerPlugin(new HalfedgeToolBar());
		v.registerPlugin(new TriangulatePlugin());
		v.registerPlugin(new WebContentLoader());
		v.registerPlugin(new ThetaVisualizer());
		v.registerPlugin(new EllipticRootGenerator());
		v.registerPlugin(new EllipticModulusEngine());
		v.registerPlugin(new ConsolePlugin());
		v.registerPlugin(new VertexEditorPlugin());
		v.registerPlugin(new NodeIndexVisualizer());
		v.registerPlugin(new CutToDiskPlugin());
		v.registerPlugin(new CutAtEdgePlugin());
		v.registerPlugin(new FindPathPlugin());
		v.registerPlugin(new LookAndFeelSwitch());
		v.registerPlugin(new CrossPlatformLnF());
		v.registerPlugin(new NimbusLnF());
		v.registerPlugin(new SystemLookAndFeel());
		v.registerPlugins(HalfedgePluginFactory.createPlugins());
		v.startup();
	} 

}
