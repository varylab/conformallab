package de.varylab.discreteconformal;

import java.util.HashSet;
import java.util.Set;

import javax.swing.Icon;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.basic.ConsolePlugin;
import de.jreality.ui.JRealitySplashScreen;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.JRHalfedgeViewer;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgePluginFactory;
import de.jtem.halfedgetools.plugin.algorithm.vectorfield.CurvatureVectorFields;
import de.jtem.halfedgetools.plugin.misc.VertexEditorPlugin;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.lnfswitch.LookAndFeelSwitch;
import de.jtem.jrworkspace.plugin.lnfswitch.plugin.CrossPlatformLnF;
import de.jtem.jrworkspace.plugin.lnfswitch.plugin.NimbusLnF;
import de.jtem.jrworkspace.plugin.lnfswitch.plugin.SystemLookAndFeel;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.DiscreteRiemannPlugin;
import de.varylab.discreteconformal.plugin.EllipticImageGenerator;
import de.varylab.discreteconformal.plugin.HyperellipticCurvePlugin;
import de.varylab.discreteconformal.plugin.QuasiIsothermicPlugin;
import de.varylab.discreteconformal.plugin.SphereEqualizerPlugin;
import de.varylab.discreteconformal.plugin.algorithm.CutAtEdgePlugin;
import de.varylab.discreteconformal.plugin.algorithm.CutToDiskPlugin;
import de.varylab.discreteconformal.plugin.algorithm.FindPathPlugin;
import de.varylab.discreteconformal.plugin.image.ImageHook;
import de.varylab.discreteconformal.plugin.schottky.SchottkyPlugin;
import de.varylab.discreteconformal.plugin.visualizer.IsothermicityMeasure;
import de.varylab.discreteconformal.plugin.visualizer.ThetaVisualizer;


public class ConformalLab {


	public static Set<Plugin> createConformalPlugins() {
		Set<Plugin> s = new HashSet<Plugin>();
		s.add(new DiscreteConformalPlugin());
		s.add(new DiscreteRiemannPlugin());
		s.add(new SchottkyPlugin());
		s.add(new ThetaVisualizer());
		s.add(new EllipticImageGenerator());
		s.add(new ConsolePlugin());
		s.add(new CutToDiskPlugin());
		s.add(new CutAtEdgePlugin());
		s.add(new FindPathPlugin());
		s.add(new QuasiIsothermicPlugin());
		s.add(new SphereEqualizerPlugin());
		s.add(new IsothermicityMeasure());
		return s;
	}
	
	
	public static void main(String[] args) {
		JRViewer.setApplicationTitle("Conformal Lab");
		JRViewer v = new JRViewer();
		Icon splashImage = ImageHook.getIcon("splash01.png");
		JRealitySplashScreen splash = new JRealitySplashScreen(splashImage);
		splash.setVisible(true);
		v.setSplashScreen(splash);
//		JRHalfedgeViewer.initHalfedgeFronted();
		NativePathUtility.set("native");
		v.addBasicUI();
		v.addContentUI();
		v.setShowToolBar(true);
		v.setShowPanelSlots(true, true, true, true);
		v.addContentSupport(ContentType.Raw);
		v.setPropertiesFile("ConformalLab.jrw");
		v.setPropertiesResource(ConformalLab.class, "ConformalLab.jrw");
		v.getController().setManageLookAndFeel(true);
//		v.registerPlugin(new WebContentLoader());
		v.registerPlugin(new LookAndFeelSwitch());
		v.registerPlugin(new CrossPlatformLnF());
		v.registerPlugin(new NimbusLnF());
		v.registerPlugin(new SystemLookAndFeel());
		v.registerPlugin(new VertexEditorPlugin());
		v.registerPlugin(CurvatureVectorFields.class);
//		v.registerPlugin(new HalfedgeDebuggerPlugin());
		v.registerPlugins(createConformalPlugins());
		v.registerPlugins(HalfedgePluginFactory.createPlugins());
		v.registerPlugin(new HyperellipticCurvePlugin());
		v.startup();
		splash.setVisible(false);
		v.getPlugin(HalfedgeInterface.class).set(new CoHDS());
		v.getPlugin(HalfedgeInterface.class).setTemplateHDS(new CoHDS());
	} 

}
