package de.varylab.discreteconformal;

import java.awt.EventQueue;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import javax.swing.Icon;

import org.pushingpixels.substance.api.SubstanceLookAndFeel;
import org.pushingpixels.substance.api.skin.ModerateSkin;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.basic.ConsolePlugin;
import de.jreality.ui.JRealitySplashScreen;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.JRHalfedgeViewer;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgePluginFactory;
import de.jtem.halfedgetools.plugin.algorithm.geometry.PrimitivesGenerator;
import de.jtem.halfedgetools.plugin.algorithm.vectorfield.CurvatureVectorFields;
import de.jtem.halfedgetools.plugin.misc.VertexEditorPlugin;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.lnfswitch.plugin.SystemLookAndFeel;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.DiscreteRiemannPlugin;
import de.varylab.discreteconformal.plugin.EllipticImageGenerator;
import de.varylab.discreteconformal.plugin.HyperellipticCurvePlugin;
import de.varylab.discreteconformal.plugin.KoebePolyhedronPlugin;
import de.varylab.discreteconformal.plugin.QuasiIsothermicPlugin;
import de.varylab.discreteconformal.plugin.SphereEqualizerPlugin;
import de.varylab.discreteconformal.plugin.algorithm.CutAtEdgePlugin;
import de.varylab.discreteconformal.plugin.algorithm.CutToDiskPlugin;
import de.varylab.discreteconformal.plugin.algorithm.FindPathPlugin;
import de.varylab.discreteconformal.plugin.image.ImageHook;
import de.varylab.discreteconformal.plugin.schottky.SchottkyPlugin;
import de.varylab.discreteconformal.plugin.visualizer.FlippedTriangles;
import de.varylab.discreteconformal.plugin.visualizer.IsothermicityMeasure;
import de.varylab.discreteconformal.plugin.visualizer.ThetaVisualizer;


public class ConformalLab implements Runnable {

	static {
		NativePathUtility.set("native");
		String[] taoCommand = new String[] {
			"-help",
//			"-snes_view",
//			"-snes_type", "ls",
//			"-snes_test_display",
//			"-pc_type", "none",
			"-pc_factor_nonzeros_along_diagonal",
			"-pc_factor_diagonal_fill",
			"-pc_factor_shift_nonzero", "1.0e-10",			
			"-tao_nm_lamda", "0.01", 
			"-tao_nm_mu", "1.0"
		};
		System.out.println("initing tao: " + Arrays.toString(taoCommand));
		Tao.Initialize("Quasiisothermic Parametrization", taoCommand, false);
	}

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
		s.add(new KoebePolyhedronPlugin());
		s.add(new PrimitivesGenerator());
		return s;
	}
	
	public static void installLookAndFeel() {
		try {
			SubstanceLookAndFeel.setSkin(new ModerateSkin());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	@Override
	public void run() {
		JRViewer.setApplicationTitle("Conformal Lab");
		JRViewer v = new JRViewer();
		installLookAndFeel();
		Icon splashImage = ImageHook.getIcon("splash01.png");
		JRealitySplashScreen splash = new JRealitySplashScreen(splashImage);
		splash.setVisible(true);
		v.setSplashScreen(splash);
		JRHalfedgeViewer.initHalfedgeFronted();
		v.addBasicUI();
		v.addContentUI();
		v.setShowToolBar(true);
		v.setShowPanelSlots(true, true, true, true);
		v.addContentSupport(ContentType.Raw);
		v.setPropertiesFile("ConformalLab.jrw");
		v.setPropertiesResource(ConformalLab.class, "ConformalLab.jrw");
		v.getController().setManageLookAndFeel(false);
//		v.registerPlugin(new WebContentLoader());
		
//		v.registerPlugin(new LookAndFeelSwitch());
//		v.registerPlugin(new CrossPlatformLnF());
//		v.registerPlugin(new NimbusLnF());
		v.registerPlugin(new SystemLookAndFeel());
		v.registerPlugin(new VertexEditorPlugin());
		v.registerPlugin(CurvatureVectorFields.class);
		v.registerPlugin(FlippedTriangles.class);
		
//		v.registerPlugin(new HalfedgeDebuggerPlugin());
		v.registerPlugins(createConformalPlugins());
		v.registerPlugins(HalfedgePluginFactory.createPlugins());
		v.registerPlugin(new HyperellipticCurvePlugin());
		v.startup();
		splash.setVisible(false);
		v.getPlugin(HalfedgeInterface.class).set(new CoHDS());
		v.getPlugin(HalfedgeInterface.class).setTemplateHDS(new CoHDS());
		
//		CurvatureLines.addCurvatureLineAdapters(v.getPlugin(HalfedgeInterface.class));
//		CurvatureLines.addBasicAdapters(v.getPlugin(HalfedgeInterface.class));		
	}
	
	
	public static void main(String[] args) {
		EventQueue.invokeLater(new ConformalLab());
	} 

}
