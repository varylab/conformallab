package de.varylab.discreteconformal.plugin;

import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.AbstractAction;

import de.jreality.plugin.AlignedContent;
import de.jreality.plugin.Background;
import de.jreality.plugin.CameraStand;
import de.jreality.plugin.ContentAppearance;
import de.jreality.plugin.ContentLoader;
import de.jreality.plugin.ContentTools;
import de.jreality.plugin.Export;
import de.jreality.plugin.HeadUpDisplay;
import de.jreality.plugin.Lights;
import de.jreality.plugin.View;
import de.jreality.plugin.ViewMenuBar;
import de.jreality.plugin.ViewPreferences;
import de.jreality.plugin.ZoomTool;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.varylab.jrworkspace.plugin.lookandfeelswitch.LookAndFeelSwitch;
import de.varylab.jrworkspace.plugin.lookandfeelswitch.plugin.CrossPlatformLnF;
import de.varylab.jrworkspace.plugin.lookandfeelswitch.plugin.FHLookAndFeel;
import de.varylab.jrworkspace.plugin.lookandfeelswitch.plugin.JGoodiesPlasticLnF;
import de.varylab.jrworkspace.plugin.lookandfeelswitch.plugin.NimbusLnF;
import de.varylab.jrworkspace.plugin.lookandfeelswitch.plugin.SubstanceLnF;
import de.varylab.jrworkspace.plugin.lookandfeelswitch.plugin.SyntheticaStandardLnf;
import de.varylab.jrworkspace.plugin.lookandfeelswitch.plugin.SystemLookAndFeel;
import de.varylab.jrworkspace.plugin.simplecontroller.SimpleController;

public class ConformalPluginTest {
	
	private static class ExitAction extends AbstractAction {

		private static final long 
			serialVersionUID = 1L;

		public ExitAction() {
			putValue(AbstractAction.NAME, "Exit");
		}
		
		@Override
		public void actionPerformed(ActionEvent e) {
			System.exit(0);
		}
		
	}

	
	public static void main(String[] args) {
		File props = new File("workspace.xml");
		SimpleController c = new SimpleController(props);

		ViewMenuBar viewerMenu = new ViewMenuBar();
		viewerMenu.addMenuSeparator(ConformalPluginTest.class, 19.0, "File");
		viewerMenu.addMenuItem(ConformalPluginTest.class, 20.0, new ExitAction(), "File");
		
		c.registerPlugin(new DiscreteConformalPlugin());
		c.registerPlugin(new HalfedgeConnectorPlugin());
		c.registerPlugin(new View());
		c.registerPlugin(new HeadUpDisplay());
		c.registerPlugin(new CameraStand());
		c.registerPlugin(new Lights());
		c.registerPlugin(new Background());
		c.registerPlugin(viewerMenu);
		c.registerPlugin(new AlignedContent());
		c.registerPlugin(new ContentLoader());
		c.registerPlugin(new ViewPreferences());
		c.registerPlugin(new ContentAppearance());
		c.registerPlugin(new ContentTools());
		c.registerPlugin(new Export());
		c.registerPlugin(new ZoomTool());
		
		// look and feel switch
		c.registerPlugin(new LookAndFeelSwitch());
		c.registerPlugin(new CrossPlatformLnF());
		c.registerPlugin(new JGoodiesPlasticLnF());
		c.registerPlugin(new NimbusLnF());
		c.registerPlugin(new SystemLookAndFeel());
		c.registerPlugin(new FHLookAndFeel());
		c.registerPlugin(new SyntheticaStandardLnf());
		c.registerPlugin(new SubstanceLnF());
		
//		Avatar avatarPlugin = new Avatar();
//		avatarPlugin.setShowPanel(false);
//		c.registerPlugin(avatarPlugin);
//		c.registerPlugin(new Terrain());
//		c.registerPlugin(new Audio());
//		c.registerPlugin(new Physics());
//		c.registerPlugin(new ContentPhysics());
//		c.registerPlugin(new ProbeBodies());
//		c.registerPlugin(new ApplyImpulse());
		
		c.startup();
	}
}
