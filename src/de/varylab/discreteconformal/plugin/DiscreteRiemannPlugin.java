package de.varylab.discreteconformal.plugin;

import static de.jtem.halfedge.util.HalfEdgeUtils.getGenus;
import static javax.swing.JOptionPane.showMessageDialog;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;
import java.util.Set;

import javax.swing.JButton;

import de.jreality.plugin.basic.View;
import de.jreality.ui.LayoutFactory;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeSelection;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.adapter.EuclideanLengthWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.HomologyUtility;

public class DiscreteRiemannPlugin extends ShrinkPanelPlugin implements ActionListener {

	private HalfedgeInterface
		hif = null;
	private JButton
		calculateButton = new JButton("Calculate Differentials"); 
	
	public DiscreteRiemannPlugin() {
		shrinkPanel.setLayout(new GridBagLayout());
		GridBagConstraints cr = LayoutFactory.createRightConstraint();
		shrinkPanel.add(calculateButton, cr);
		
		calculateButton.addActionListener(this);
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		CoHDS S = hif.get(new CoHDS());
		int g = getGenus(S);
		if (g == 0) {
			showMessageDialog(
				shrinkPanel, 
				"No holomorphic differentials on genus 0 surfaces"
			);
			return;
		}
		CoVertex root = S.getVertex(0);
		EuclideanLengthWeightAdapter w = new EuclideanLengthWeightAdapter(null);
		List<Set<CoEdge>> paths = HomologyUtility.getGeneratorPaths(root, w);
		System.out.println("Got " + paths.size() + " paths");
		HalfedgeSelection s = new HalfedgeSelection();
		for (Set<CoEdge> path : paths) {
			s.addAll(path);
		}
		for (Set<CoEdge> path : paths) {
			EdgeVectorAdapter eva = new EdgeVectorAdapter(path, "Homology Path " + path.size());
			hif.addLayerAdapter(eva, false);
		}
		
		hif.setSelection(s);
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

}
