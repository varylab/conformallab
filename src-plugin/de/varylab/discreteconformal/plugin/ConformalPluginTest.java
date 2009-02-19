package de.varylab.discreteconformal.plugin;

import java.io.File;

import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

import com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel;

import de.jreality.plugin.view.AlignedContent;
import de.jreality.plugin.view.Background;
import de.jreality.plugin.view.CameraStand;
import de.jreality.plugin.view.ContentAppearance;
import de.jreality.plugin.view.ContentLoader;
import de.jreality.plugin.view.ContentTools;
import de.jreality.plugin.view.Export;
import de.jreality.plugin.view.Lights;
import de.jreality.plugin.view.View;
import de.jreality.plugin.view.ViewMenuBar;
import de.jreality.plugin.view.ViewPreferences;
import de.jreality.plugin.view.ZoomTool;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.jtem.halfedge.plugin.HalfedgeDebuggerPlugin;
import de.jtem.halfedge.plugin.buildin.CatmullClarkPlugin;
import de.jtem.halfedge.plugin.buildin.TriangulatePlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.jrworkspace.plugin.simplecontroller.SimpleController;

public class ConformalPluginTest {
	
	static {
		try {
			UIManager.setLookAndFeel(new NimbusLookAndFeel());
		} catch (UnsupportedLookAndFeelException e) {
			e.printStackTrace();
		}
	}
	
	
	public static void main(String[] args) {
		File props = new File("workspace.xml");
		SimpleController c = new SimpleController(props);
		c.setManageLookAndFeel(false);
		
		// conformal maps
		c.registerPlugin(new DiscreteConformalPlugin());
		c.registerPlugin(new HalfedgeConnectorPlugin());
		c.registerPlugin(new HalfedgeDebuggerPlugin<CoVertex, CoEdge, CoFace>());
		c.registerPlugin(new CatmullClarkPlugin());
		c.registerPlugin(new TriangulatePlugin());
		
		// standard viewer
		c.registerPlugin(new View());
		c.registerPlugin(new CameraStand());
		c.registerPlugin(new Lights());
		c.registerPlugin(new Background());
		c.registerPlugin(new ViewMenuBar());
		c.registerPlugin(new AlignedContent());
		c.registerPlugin(new ContentLoader());
		c.registerPlugin(new ViewPreferences());
		c.registerPlugin(new ContentAppearance());
		c.registerPlugin(new ContentTools());
		c.registerPlugin(new Export());
		c.registerPlugin(new ZoomTool());
		
		c.startup();
	}
}
