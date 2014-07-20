package de.varylab.discreteconformal.plugin;

import de.jtem.halfedgetools.plugin.texturespace.TextureSpacePlugin;
import de.jtem.java2d.SceneComponent;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.sidecontainer.widget.ShrinkPanel;

public class UniformizationTextureSpacePlugin extends Plugin implements TextureSpacePlugin {

	private ShrinkPanel
		options = new ShrinkPanel("Uniformization");
	private SceneComponent
		scene = new SceneComponent();
	
	
	
	@Override
	public SceneComponent getSceneComponent() {
		return scene;
	}

	@Override
	public ShrinkPanel getOptionPanel() {
		return options;
	}

	@Override
	public boolean getRenderOnTop() {
		return true;
	}

}
