package de.varylab.discreteconformal;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.LogManager;

import javax.swing.JPopupMenu;
import javax.swing.UIManager;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.basic.ConsolePlugin;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedgetools.JRHalfedgeViewer;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgePluginFactory;
import de.jtem.halfedgetools.plugin.algorithm.generator.PrimitivesGenerator;
import de.jtem.halfedgetools.plugin.algorithm.vectorfield.CurvatureVectorFields;
import de.jtem.halfedgetools.plugin.misc.VertexEditorPlugin;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.lnfswitch.plugin.SystemLookAndFeel;
import de.jtem.jrworkspace.plugin.simplecontroller.StartupChain;
import de.jtem.jrworkspace.plugin.simplecontroller.widget.SplashScreen;
import de.varylab.discreteconformal.datasource.CylinderEdgesDataSource;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.plugin.ConformalDataPlugin;
import de.varylab.discreteconformal.plugin.ConformalVisualizationPlugin;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.DiscreteRiemannPlugin;
import de.varylab.discreteconformal.plugin.DomainVisualisationPlugin;
import de.varylab.discreteconformal.plugin.EllipticImageGenerator;
import de.varylab.discreteconformal.plugin.HyperellipticCurvePlugin;
import de.varylab.discreteconformal.plugin.KoebePolyhedronPlugin;
import de.varylab.discreteconformal.plugin.ProjectiveTexturePlugin;
import de.varylab.discreteconformal.plugin.QuasiIsothermicPlugin;
import de.varylab.discreteconformal.plugin.SphereEqualizerPlugin;
import de.varylab.discreteconformal.plugin.algorithm.CutAndGlueConformalDomain;
import de.varylab.discreteconformal.plugin.algorithm.CutAtEdgePlugin;
import de.varylab.discreteconformal.plugin.algorithm.CutToDiskPlugin;
import de.varylab.discreteconformal.plugin.algorithm.FindPathPlugin;
import de.varylab.discreteconformal.plugin.schottky.SchottkyPlugin;
import de.varylab.discreteconformal.plugin.visualizer.FlippedTriangles;
import de.varylab.discreteconformal.plugin.visualizer.IndexMedialGraph;
import de.varylab.discreteconformal.plugin.visualizer.IsothermicityMeasure;
import de.varylab.discreteconformal.plugin.visualizer.ThetaVisualizer;
import de.varylab.discreteconformal.startup.ConformalLabSplashScreen;


public class ConformalLab implements Runnable {

	static {
		NativePathUtility.set("native");
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
		s.add(new ProjectiveTexturePlugin());
		s.add(new ConformalVisualizationPlugin());
		s.add(new DomainVisualisationPlugin());
		s.add(new IndexMedialGraph());
		s.add(new ConformalDataPlugin());
		s.add(new CylinderEdgesDataSource());
		s.add(new CutAndGlueConformalDomain());
		return s;
	}
	
	public static void installLookAndFeel() {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
//			SubstanceLookAndFeel.setSkin(new GraphiteAquaSkin());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	@Override
	public void run() {
		final JRViewer v = new JRViewer();
		final SplashScreen splash = new ConformalLabSplashScreen();
		splash.setVisible(true);
		v.setSplashScreen(splash);

		Runnable jobStaticInit = new Runnable() {
			@Override
			public void run() {
				JRViewer.setApplicationTitle("Discrete Conformal Lab");
				JPopupMenu.setDefaultLightWeightPopupEnabled(false);
				JRHalfedgeViewer.initHalfedgeFronted();
				installLookAndFeel();
			}
		};		
		
		Runnable jobInitViewer = new Runnable() {
			@Override
			public void run() {
				v.setSplashScreen(splash);
				v.addBasicUI();
				v.addContentUI();
				v.addPythonSupport();
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
			}
		};
		
		StartupChain initChain = new StartupChain();
		initChain.appendJob(jobStaticInit);
		initChain.appendJob(jobInitViewer);
		initChain.startQueuedAndWait();
		
		v.startup();
		
		splash.setVisible(false);
		v.getPlugin(HalfedgeInterface.class).set(new CoHDS());
		v.getPlugin(HalfedgeInterface.class).setTemplateHDS(new CoHDS());
		
//		CurvatureLines.addCurvatureLineAdapters(v.getPlugin(HalfedgeInterface.class));
//		CurvatureLines.addBasicAdapters(v.getPlugin(HalfedgeInterface.class));		
	}
	
	
	public static void main(final String[] args) throws Exception {
		initLogging();
		if (args.length == 0) { // gui mode
			new ConformalLab().run();
		} else { // batch mode
			ConformalLabBatch cl = new ConformalLabBatch();
			cl.process(args);
		}
	} 
	
	public static void initLogging() {
		LogManager lm = LogManager.getLogManager();
		File localConf = new File("logging.properties");
		File localCustomConf = new File("logcustom.properties");
		@SuppressWarnings("unused")
		String fileName = null;
		InputStream confIn = null;
		if (localConf.exists() | localCustomConf.exists()) {
			File confFile = null;
			if (localCustomConf.exists()) {
				confFile = localCustomConf;
			} else {
				confFile = localConf;
			}
			fileName = confFile.getAbsolutePath();
			try {
				confIn = new FileInputStream(confFile);
			} catch (FileNotFoundException e) {}
		} else {
			fileName = ConformalLab.class.getResource("logging.properties").getFile();
			confIn = ConformalLab.class.getResourceAsStream("logging.properties");
		}
		assert confIn != null;
		try {
			lm.readConfiguration(confIn);
		} catch (Exception e) {
			System.out.println(e);
		}
	}

}
