package de.varylab.discreteconformal.plugin;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.ProgressMonitor;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import de.jtem.halfedge.jreality.adapter.Adapter;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.adapter.TexCoordAdapter;
import de.varylab.discreteconformal.unwrapper.CDiskUnwrapper;
import de.varylab.discreteconformal.unwrapper.CDiskUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapper;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.CUnwrapper;
import de.varylab.discreteconformal.unwrapper.UnwrapException;
import de.varylab.jrworkspace.plugin.sidecontainer.widget.ShrinkPanel;

public class UnwrapShrinker extends ShrinkPanel implements ActionListener, ChangeListener {

	private static final long 
		serialVersionUID = 1L;
	private JButton
		unwrapBtn = new JButton("Unwrap");
	private JPanel
		coneConfigPanel = new JPanel(),
		geometryPanel = new JPanel();
	private SpinnerNumberModel
		numConesModel = new SpinnerNumberModel(0, 0, 100, 1);
	private JSpinner
		numConesSpinner = new JSpinner(numConesModel);
	private JCheckBox
		quantizeChecker = new JCheckBox("Quantize Cone Angles");
	private JRadioButton
		euclideanButton = new JRadioButton("Eucliean"),
		hyperbolicButton = new JRadioButton("Hyperbolic");
	private JComboBox
		numericsCombo = new JComboBox(new String[] {"Java/MTJ Numerics", "Petsc/Tao Numerics"});
	
	private int
		numCones = 0; 
	private boolean
		quantizeCones = true,
		usePetsc = false,
		hyperbolic = true;

	private HalfedgeConnectorPlugin
		halfedgeConnectorPlugin = null;
	
	
	public UnwrapShrinker(HalfedgeConnectorPlugin hcp) {
		super("Unwrap");
		this.halfedgeConnectorPlugin = hcp;
		createLayout();
		unwrapBtn.addActionListener(this);
		numConesSpinner.addChangeListener(this);
		quantizeChecker.addActionListener(this);
		numericsCombo.addActionListener(this);
		euclideanButton.addActionListener(this);
		hyperbolicButton.addActionListener(this);
		
		ButtonGroup geometryGroup = new ButtonGroup();
		geometryGroup.add(euclideanButton);
		geometryGroup.add(hyperbolicButton);
	}

	
	
	
	private void createLayout() {
		setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(2,2,2,2);
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 1.0;
		c.gridwidth = GridBagConstraints.REMAINDER;
		
		coneConfigPanel.setBorder(BorderFactory.createTitledBorder("Cone Singularities"));
		coneConfigPanel.setLayout(new GridBagLayout());
		add(coneConfigPanel, c);
		
		JLabel numConesLabel = new JLabel("Cones");
		c.weightx = 0.0;
		c.gridwidth = GridBagConstraints.RELATIVE;
		coneConfigPanel.add(numConesLabel, c);
		c.weightx = 1.0;
		c.gridwidth = GridBagConstraints.REMAINDER;
		coneConfigPanel.add(numConesSpinner, c);
		quantizeChecker.setSelected(true);
		coneConfigPanel.add(quantizeChecker, c);
		
		geometryPanel.setBorder(BorderFactory.createTitledBorder("Geometry"));
		geometryPanel.setLayout(new GridLayout(1, 2));
		add(geometryPanel, c);
		euclideanButton.setSelected(!hyperbolic);
		geometryPanel.add(euclideanButton);
		hyperbolicButton.setSelected(hyperbolic);
		geometryPanel.add(hyperbolicButton);
		
		numericsCombo.setLightWeightPopupEnabled(true);
		numericsCombo.setSelectedIndex(0);
		add(numericsCombo, c);
		add(unwrapBtn, c);
	}


	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object s = e.getSource();
		if (unwrapBtn == s) {
			new ComputationThread().execute();
		}
		if (quantizeChecker == s) {
			quantizeCones = quantizeChecker.isSelected();
		}
		if (numericsCombo == s) {
			usePetsc = numericsCombo.getSelectedIndex() != 0;
		}
		if (euclideanButton == s) {
			hyperbolic = !euclideanButton.isSelected();
		}
		if (hyperbolicButton == s) {
			hyperbolic = hyperbolicButton.isSelected();
		}
	}
	
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (numConesSpinner == e.getSource()) {
			numCones = numConesModel.getNumber().intValue();
		}
	}
	
	
	private class ComputationThread extends SwingWorker<Object, Object>{

		private ProgressMonitor
			mon = null;
		
		
		public ComputationThread() {
			Window w = SwingUtilities.getWindowAncestor(UnwrapShrinker.this);
			mon = new ProgressMonitor(w, "Unwrapping", "Initializing", 0, 5);
		}
		
		
		@Override
		public Object doInBackground() {
			System.out.println("Unwrapping...");
			CoHDS hds = new CoHDS();
			hds.setTexCoordinatesValid(false);
			Adapter ca = new PositionAdapter();
			halfedgeConnectorPlugin.getHalfedgeContent(hds, ca);
			
			CUnwrapper unwrapper = null;
			if (hyperbolic) {
				if(usePetsc) {
					unwrapper = new CHyperbolicUnwrapperPETSc();
				} else {
					unwrapper = new CHyperbolicUnwrapper();
				}
			} else {
				if(usePetsc) {
					unwrapper = new CDiskUnwrapperPETSc(numCones, quantizeCones);
				} else {
					unwrapper = new CDiskUnwrapper(numCones, quantizeCones);
				}
			}
			
			// unwrap
			try {
				unwrapper.unwrap(hds, null);
			} catch (UnwrapException ue) {
				ue.printStackTrace();
				mon.close();
				return null;
			}
			
			hds.setTexCoordinatesValid(true);			
			TexCoordAdapter tca = new TexCoordAdapter();
			halfedgeConnectorPlugin.updateHalfedgeContent(hds, true, ca, tca);
			return null;
		}
		
	}
	
}