package de.varylab.discreteconformal.plugin;

import de.jreality.scene.pick.PickResult;
import de.jreality.scene.tool.AbstractTool;
import de.jreality.scene.tool.InputSlot;
import de.jreality.scene.tool.ToolContext;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;

public class HyperbolicCopyTool extends AbstractTool {

	private TextureSpaceViewer3D
		domainPlugin = null;
	private HalfedgeLayer
		layer = null;
	
	public HyperbolicCopyTool(TextureSpaceViewer3D dcp, HalfedgeLayer layer) {
		super(InputSlot.getDevice("PrimaryAction"));
		this.domainPlugin = dcp;
		this.layer = layer;
	}
	
	@Override
	public void activate(ToolContext tc) {
		PickResult pr = tc.getCurrentPick();
		if (pr == null || PickResult.PICK_TYPE_FACE != pr.getPickType()) {
			return;
		}
		domainPlugin.copyDomainAtEdge(pr.getIndex(), layer);
	}
	
}