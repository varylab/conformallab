package de.varylab.discreteconformal.plugin;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import de.jreality.plugin.basic.View;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.plugin.image.ImageHook;

public class ConformalVisualizationPlugin extends ShrinkPanelPlugin implements ActionListener {

	private HalfedgeInterface
		hif = null;
	private DiscreteConformalPlugin
		dcp = null;
	
	private JComboBox
		interpolationCombo = new JComboBox(InterpolationMethod.values());
	private JRadioButton
		kleinButton = new JRadioButton("Klein"),
		poincareButton = new JRadioButton("Poincaré", true),
		halfplaneButton = new JRadioButton("Half-Plane"); 	
	private JPanel
		modelPanel = new JPanel();
	
	public ConformalVisualizationPlugin() {
		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(2, 2, 2, 2);
		c.weightx = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;
		shrinkPanel.setLayout(new GridBagLayout());
		c.gridwidth = GridBagConstraints.RELATIVE;
		shrinkPanel.add(new JLabel("Interpolation"), c);
		c.gridwidth = GridBagConstraints.REMAINDER;
		shrinkPanel.add(interpolationCombo, c);
		shrinkPanel.add(modelPanel, c);
		
		modelPanel.setBorder(BorderFactory.createTitledBorder("Hyperbolic Model"));
		modelPanel.setLayout(new GridBagLayout()); c.gridwidth = 1;
		modelPanel.add(kleinButton, c);
		modelPanel.add(poincareButton, c);
		modelPanel.add(halfplaneButton, c);
		
		interpolationCombo.addActionListener(this);
		kleinButton.addActionListener(this);
		poincareButton.addActionListener(this);
		halfplaneButton.addActionListener(this);

		ButtonGroup modelGroup = new ButtonGroup();
		modelGroup.add(kleinButton);
		modelGroup.add(poincareButton);
		modelGroup.add(halfplaneButton);
		
		shrinkPanel.setIcon(ImageHook.getIcon("palette.png"));
	}
	
	
	public InterpolationMethod getSelectedInterpolation() {
		return (InterpolationMethod)interpolationCombo.getSelectedItem();
	}

	public HyperbolicModel getSelectedHyperbolicModel() {
		if (kleinButton.isSelected()) {
			return HyperbolicModel.Klein;
		}
		if (poincareButton.isSelected()) {
			return HyperbolicModel.Poincaré;
		}
		if (halfplaneButton.isSelected()) {
			return HyperbolicModel.Halfplane;
		}
		return HyperbolicModel.Klein;
	}
	
	public void setInterpolation(InterpolationMethod method) {
		interpolationCombo.setSelectedItem(method);
		updateAdapters();
	}

	public void setHyperbolicModel(HyperbolicModel model) {
		kleinButton.setSelected(model == HyperbolicModel.Klein);
		poincareButton.setSelected(model == HyperbolicModel.Poincaré);
		halfplaneButton.setSelected(model == HyperbolicModel.Halfplane);
		updateAdapters();
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		updateAdapters();
		hif.updateNoUndo();
	}


	public void updateAdapters() {
		CoTexturePositionAdapter texturePositionAdapter = dcp.texturePositionAdapter;
		texturePositionAdapter.setInterpolationMethod(getSelectedInterpolation());
		texturePositionAdapter.setModel(getSelectedHyperbolicModel());
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		dcp = c.getPlugin(DiscreteConformalPlugin.class);
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

}
