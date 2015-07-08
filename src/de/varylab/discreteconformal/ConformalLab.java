package de.varylab.discreteconformal;

import java.awt.Image;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JPopupMenu;
import javax.swing.UIDefaults;
import javax.swing.UIManager;

import org.pushingpixels.substance.api.SubstanceLookAndFeel;
import org.pushingpixels.substance.api.fonts.FontPolicy;
import org.pushingpixels.substance.api.fonts.FontSet;
import org.pushingpixels.substance.api.skin.GraphiteSkin;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.basic.ConsolePlugin;
import de.jreality.plugin.basic.View;
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
import de.varylab.discreteconformal.datasource.ConicalEdgesDataSource;
import de.varylab.discreteconformal.datasource.CylinderEdgesDataSource;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.plugin.ConformalDataPlugin;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.DiscreteRiemannPlugin;
import de.varylab.discreteconformal.plugin.HyperIdealPlugin;
import de.varylab.discreteconformal.plugin.HyperIdealVisualizationPlugin;
import de.varylab.discreteconformal.plugin.HyperellipticCurveGenerator;
import de.varylab.discreteconformal.plugin.HyperellipticCurvePlugin;
import de.varylab.discreteconformal.plugin.KoebePolyhedronPlugin;
import de.varylab.discreteconformal.plugin.ProjectiveTexturePlugin;
import de.varylab.discreteconformal.plugin.QuasiIsothermicPlugin;
import de.varylab.discreteconformal.plugin.SphericalEqualizerPlugin;
import de.varylab.discreteconformal.plugin.SphericalNormalizationPlugin;
import de.varylab.discreteconformal.plugin.TextureSpaceViewer3D;
import de.varylab.discreteconformal.plugin.UniformizationDomainPlugin;
import de.varylab.discreteconformal.plugin.algorithm.AddConeOfRevolutionCommand;
import de.varylab.discreteconformal.plugin.algorithm.ContractShortEdges;
import de.varylab.discreteconformal.plugin.algorithm.CutAndGlueConformalDomain;
import de.varylab.discreteconformal.plugin.algorithm.CutAtEdgePlugin;
import de.varylab.discreteconformal.plugin.algorithm.CutToDiskPlugin;
import de.varylab.discreteconformal.plugin.algorithm.FindPathPlugin;
import de.varylab.discreteconformal.plugin.algorithm.MapToConeCommand;
import de.varylab.discreteconformal.plugin.algorithm.SelectHomotopyGeneratorsPlugin;
import de.varylab.discreteconformal.plugin.image.ImageHook;
import de.varylab.discreteconformal.plugin.schottky.SchottkyPlugin;
import de.varylab.discreteconformal.plugin.visualizer.DiscreteConformalEquivalencemMeasure;
import de.varylab.discreteconformal.plugin.visualizer.FlippedTriangles;
import de.varylab.discreteconformal.plugin.visualizer.IndexMedialGraph;
import de.varylab.discreteconformal.plugin.visualizer.IsothermicityMeasure;
import de.varylab.discreteconformal.plugin.visualizer.ThetaVisualizer;
import de.varylab.discreteconformal.startup.ConformalLabSplashScreen;
import de.varylab.discreteconformal.startup.TrebuchetFontSet;

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
		s.add(new HyperellipticCurveGenerator());
		s.add(new ConsolePlugin());
		s.add(new CutToDiskPlugin());
		s.add(new CutAtEdgePlugin());
		s.add(new FindPathPlugin());
		s.add(new QuasiIsothermicPlugin());
		s.add(new SphericalEqualizerPlugin());
		s.add(new IsothermicityMeasure());
		s.add(new KoebePolyhedronPlugin());
		s.add(new PrimitivesGenerator());
		s.add(new ProjectiveTexturePlugin());
		s.add(new TextureSpaceViewer3D());
		s.add(new IndexMedialGraph());
		s.add(new ConformalDataPlugin());
		s.add(new CylinderEdgesDataSource());
		s.add(new ConicalEdgesDataSource());
		s.add(new CutAndGlueConformalDomain());
		s.add(new MapToConeCommand());
		s.add(new AddConeOfRevolutionCommand());
		s.add(new ContractShortEdges());
		s.add(new UniformizationDomainPlugin());
		s.add(new DiscreteConformalEquivalencemMeasure());
		s.add(new SphericalNormalizationPlugin());
		s.add(new SelectHomotopyGeneratorsPlugin());
		s.add(new HyperIdealVisualizationPlugin());
		return s;
	}
	
	public static void installLookAndFeel() {
		try {
//			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			SubstanceLookAndFeel.setSkin(new GraphiteSkin());
			SubstanceLookAndFeel.setToUseConstantThemesOnDialogs(true);
			JFrame.setDefaultLookAndFeelDecorated(true);
			UIManager.put(SubstanceLookAndFeel.SHOW_EXTRA_WIDGETS, Boolean.TRUE);
			FontPolicy newFontPolicy = new FontPolicy() {
				@Override
				public FontSet getFontSet(String lafName, UIDefaults table) {
					return new TrebuchetFontSet(12);
				}
			};
			SubstanceLookAndFeel.setFontPolicy(newFontPolicy);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static List<Image> getMainIconList() {
		List<Image> iconList = new LinkedList<Image>();
		iconList.add(ImageHook.getImage("logo16.png"));
		iconList.add(ImageHook.getImage("logo24.png"));
		iconList.add(ImageHook.getImage("logo32.png"));
		iconList.add(ImageHook.getImage("logo48.png"));
		iconList.add(ImageHook.getImage("logo64.png"));
		iconList.add(ImageHook.getImage("logo128.png"));
		iconList.add(ImageHook.getImage("logo256.png"));
		iconList.add(ImageHook.getImage("logo512.png"));
		iconList.add(ImageHook.getImage("logo1024.png"));
		return iconList;
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
				v.registerPlugin(HyperIdealPlugin.class);
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
		View.setTitle("Discrete Conformal Lab");
		View.setIcon(ImageHook.getIcon("logo512.png"));
		View.setIconList(getMainIconList());
		JRViewer.setApplicationTitle("Discrete Conformal Lab");
		JRViewer.setApplicationIcon(ImageHook.getImage("logo512.png"));
		LoggingUtility.initLogging();
		if (args.length == 0) { // gui mode
			new ConformalLab().run();
		} else { // batch mode
			ConformalLabBatch cl = new ConformalLabBatch();
			cl.process(args);
		}
	}

}
