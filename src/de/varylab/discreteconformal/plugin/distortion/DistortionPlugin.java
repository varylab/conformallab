package de.varylab.discreteconformal.plugin.distortion;

import java.awt.GridLayout;

import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmDialogPlugin;

public class DistortionPlugin extends AlgorithmDialogPlugin {
	
	private JPanel
		dialogPanel = new JPanel();
	private JRadioButton
		xButton = new JRadioButton("X"),
		yButton = new JRadioButton("Y"),
		zButton = new JRadioButton("Z", true);
	private ButtonGroup
		group = new ButtonGroup();
	
	public DistortionPlugin() {
		group.add(xButton);
		group.add(yButton);
		group.add(zButton);
		dialogPanel.setLayout(new GridLayout(3, 1));
		dialogPanel.add(xButton);
		dialogPanel.add(yButton);
		dialogPanel.add(zButton);
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void executeAfterDialog(HDS hds, AdapterSet a, HalfedgeInterface hi) throws Exception {
		for (V v : hds.getVertices()) {
			// fix normal
			a.set(Normal.class, v, a.getD(Normal.class, v));
			// project to cylinder
			double[] pos = a.getD(Position3d.class, v);
			if (xButton.isSelected()) {
				double[] pos2d = new double[] {pos[1], pos[2]};
				Rn.normalize(pos2d, pos2d);
				pos[1] = pos2d[0]; pos[2] = pos2d[1];
			}
			if (yButton.isSelected()) {
				double[] pos2d = new double[] {pos[0], pos[2]};
				Rn.normalize(pos2d, pos2d);
				pos[0] = pos2d[0]; pos[2] = pos2d[1];
			}
			if (zButton.isSelected()) {
				double[] pos2d = new double[] {pos[0], pos[1]};
				Rn.normalize(pos2d, pos2d);
				pos[0] = pos2d[0]; pos[1] = pos2d[1];
			}
			a.set(Position.class, v, pos);
		}
		hi.update();
	}
	
	@Override
	protected JPanel getDialogPanel() {
		return dialogPanel;
	}
	
	@Override
	public String getAlgorithmName() {
		return "Project to Cylinder";
	}
	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.DDG;
	}

}
