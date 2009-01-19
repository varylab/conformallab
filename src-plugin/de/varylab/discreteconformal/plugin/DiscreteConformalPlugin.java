package de.varylab.discreteconformal.plugin;

import java.awt.GridLayout;

import de.jreality.plugin.View;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.varylab.jrworkspace.plugin.Controller;
import de.varylab.jrworkspace.plugin.PluginInfo;
import de.varylab.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.varylab.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin {

	private UnwrapShrinker
		unwrapShrinker = null;
	
	public DiscreteConformalPlugin() {

	}
	
	@Override
	public void install(Controller c) throws Exception {
		unwrapShrinker = new UnwrapShrinker(c.getPlugin(HalfedgeConnectorPlugin.class));
		shrinkPanel.setLayout(new GridLayout());
		shrinkPanel.add(unwrapShrinker.getContentPanel());
		super.install(c);
	}
	
	
	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = new PluginInfo();
		info.name = "Discrete Conformal Parametrization";
		info.vendorName = "Stefan Sechelmann";
		info.email = "sechel@math.tu-berlin.de";
		return info;
	}

	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

}
