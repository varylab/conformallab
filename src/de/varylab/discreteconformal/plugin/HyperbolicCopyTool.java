package de.varylab.discreteconformal.plugin;

import static de.jreality.scene.pick.PickResult.PICK_TYPE_LINE;
import de.jreality.scene.pick.PickResult;
import de.jreality.scene.tool.AbstractTool;
import de.jreality.scene.tool.InputSlot;
import de.jreality.scene.tool.ToolContext;

public class HyperbolicCopyTool extends AbstractTool {

	private DiscreteConformalPlugin
		discreteConformalPlugin = null;
	
	public HyperbolicCopyTool(DiscreteConformalPlugin dcp) {
		super(InputSlot.getDevice("PrimaryAction"));
		discreteConformalPlugin = dcp;
	}
	
	
	@Override
	public void activate(ToolContext tc) {
		PickResult pr = tc.getCurrentPick();
		if (pr == null || PICK_TYPE_LINE != pr.getPickType()) {
			return;
		}
		discreteConformalPlugin.copyAtEdge(pr.getIndex());
	}
	
	
	
}
