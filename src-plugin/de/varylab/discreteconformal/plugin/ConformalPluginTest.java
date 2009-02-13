package de.varylab.discreteconformal.plugin;

import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.AbstractAction;

import de.jreality.ui.plugin.view.AlignedContent;
import de.jreality.ui.plugin.view.Background;
import de.jreality.ui.plugin.view.CameraStand;
import de.jreality.ui.plugin.view.ContentAppearance;
import de.jreality.ui.plugin.view.ContentLoader;
import de.jreality.ui.plugin.view.ContentTools;
import de.jreality.ui.plugin.view.Export;
import de.jreality.ui.plugin.view.Inspector;
import de.jreality.ui.plugin.view.Lights;
import de.jreality.ui.plugin.view.View;
import de.jreality.ui.plugin.view.ViewMenuBar;
import de.jreality.ui.plugin.view.ViewPreferences;
import de.jreality.ui.plugin.view.ZoomTool;
import de.jreality.vr.plugin.HeadUpDisplay;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
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
		c.registerPlugin(new Inspector());
		
		c.startup();
	}
}
