package de.varylab.discreteconformal.plugin;

import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static java.awt.GridBagConstraints.RELATIVE;
import static java.awt.GridBagConstraints.REMAINDER;
import static java.lang.Math.PI;
import static javax.swing.JOptionPane.ERROR_MESSAGE;

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
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import no.uib.cipr.matrix.Vector;
import de.jreality.geometry.Primitives;
import de.jreality.math.MatrixBuilder;
import de.jreality.scene.Appearance;
import de.jreality.scene.proxy.scene.SceneGraphComponent;
import de.jreality.ui.plugin.AlignedContent;
import de.jreality.ui.plugin.View;
import de.jtem.halfedge.algorithm.triangulation.Triangulator;
import de.jtem.halfedge.jreality.adapter.Adapter;
import de.jtem.halfedge.plugin.HalfedgeConnectorPlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.adapter.PositionTexCoordAdapter;
import de.varylab.discreteconformal.heds.adapter.TexCoordAdapter;
import de.varylab.discreteconformal.heds.util.UniformizationUtility;
import de.varylab.discreteconformal.heds.util.UniformizationUtility.UAdapter;
import de.varylab.discreteconformal.unwrapper.CDiskUnwrapper;
import de.varylab.discreteconformal.unwrapper.CDiskUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.CHyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapper;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.CUnwrapper;
import de.varylab.discreteconformal.unwrapper.CHyperbolicLayout.HyperbolicLayoutContext;
import de.varylab.jrworkspace.plugin.Controller;
import de.varylab.jrworkspace.plugin.PluginInfo;
import de.varylab.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.varylab.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;

public class DiscreteConformalPlugin extends ShrinkPanelPlugin implements ActionListener, ChangeListener {

	private AlignedContent
		content = null;
	private HalfedgeConnectorPlugin
		hcp = null;
	private Triangulator<CoVertex, CoEdge, CoFace>
		triangulator = new Triangulator<CoVertex, CoEdge, CoFace>();
	private CoHDS
		activeGeometry = null,
		unwrappedGeometry = null;
	private MarkedEdgesAdapter
		pointColorAdapter = new MarkedEdgesAdapter();
	private SceneGraphComponent
		auxGeometry = new SceneGraphComponent(),
		unitCircle = new SceneGraphComponent();

	private JLabel
		geometryLabel = new JLabel();
	private JButton
		unwrapBtn = new JButton("Unwrap"),
		getGeometryBtn = new JButton("Retrieve Geometry"),
		reduceBtn = new JButton("To Fundamental Polygon");
	private JPanel
		coneConfigPanel = new JPanel(),
		geometryPanel = new JPanel(),
		visualizationPanel = new JPanel();
	private SpinnerNumberModel
		numConesModel = new SpinnerNumberModel(0, 0, 100, 1);
	private JSpinner
		numConesSpinner = new JSpinner(numConesModel);
	private JCheckBox
		quantizeChecker = new JCheckBox("Quantize Cone Angles"),
		showUnwrapped = new JCheckBox("Show Unwrapped Geometry"),
		showUnitCircle = new JCheckBox("Show Unit Cirlce");
	private JRadioButton
		euclideanButton = new JRadioButton("Eucliean"),
		hyperbolicButton = new JRadioButton("Hyperbolic"),
		kleinButton = new JRadioButton("Klein"),
		poincareButton = new JRadioButton("Poincaré");
	private JComboBox
		numericsCombo = new JComboBox(new String[] {"Java/MTJ Numerics", "Petsc/Tao Numerics"});

	private int
		numCones = 0; 
	private boolean
		quantizeCones = true,
		usePetsc = false,
		hyperbolic = true,
		klein = true;
	
	
	public DiscreteConformalPlugin() {
		createLayout();
		unwrapBtn.addActionListener(this);
		numConesSpinner.addChangeListener(this);
		quantizeChecker.addActionListener(this);
		numericsCombo.addActionListener(this);
		euclideanButton.addActionListener(this);
		hyperbolicButton.addActionListener(this);
		getGeometryBtn.addActionListener(this);
		showUnwrapped.addActionListener(this);
		showUnitCircle.addActionListener(this);
		kleinButton.addActionListener(this);
		poincareButton.addActionListener(this);
		reduceBtn.addActionListener(this);
		
		ButtonGroup geometryGroup = new ButtonGroup();
		geometryGroup.add(euclideanButton);
		geometryGroup.add(hyperbolicButton);
		
		ButtonGroup modelGroup = new ButtonGroup();
		modelGroup.add(kleinButton);
		modelGroup.add(poincareButton);
		
		Appearance circleApp = new Appearance();
		circleApp.setAttribute(EDGE_DRAW, false);
		circleApp.setAttribute(VERTEX_DRAW, false); 
		circleApp.setAttribute(FACE_DRAW, true);
		unitCircle.setAppearance(circleApp);
		unitCircle.setVisible(false);
		MatrixBuilder.euclidean().rotate(PI / 2, 1, 0, 0).assignTo(unitCircle);
		unitCircle.setGeometry(Primitives.torus(1.0, 0.005, 200, 5));
		auxGeometry.addChild(unitCircle);
	}

	
	
	
	private void createLayout() {
		shrinkPanel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(2,2,2,2);
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 1.0;
		c.gridwidth = GridBagConstraints.REMAINDER;
		
		shrinkPanel.add(getGeometryBtn, c);
		shrinkPanel.add(geometryLabel, c);
		
		coneConfigPanel.setBorder(BorderFactory.createTitledBorder("Cone Singularities"));
		coneConfigPanel.setLayout(new GridBagLayout());
		shrinkPanel.add(coneConfigPanel, c);
		
		JLabel numConesLabel = new JLabel("Cones");
		c.weightx = 0.0;
		c.gridwidth = RELATIVE;
		coneConfigPanel.add(numConesLabel, c);
		c.weightx = 1.0;
		c.gridwidth = REMAINDER;
		coneConfigPanel.add(numConesSpinner, c);
		quantizeChecker.setSelected(true);
		coneConfigPanel.add(quantizeChecker, c);
		
		geometryPanel.setBorder(BorderFactory.createTitledBorder("Geometry"));
		geometryPanel.setLayout(new GridLayout(1, 2));
		shrinkPanel.add(geometryPanel, c);
		geometryPanel.add(euclideanButton);
		geometryPanel.add(hyperbolicButton);
		
		numericsCombo.setLightWeightPopupEnabled(true);
		numericsCombo.setSelectedIndex(0);
		shrinkPanel.add(numericsCombo, c);
		shrinkPanel.add(unwrapBtn, c);
		
		visualizationPanel.setBorder(BorderFactory.createTitledBorder("Visualization"));
		visualizationPanel.setLayout(new GridBagLayout());
		visualizationPanel.add(showUnitCircle, c);
		visualizationPanel.add(showUnwrapped, c);
		c.gridwidth = RELATIVE;
		visualizationPanel.add(kleinButton, c);
		c.gridwidth = REMAINDER;
		visualizationPanel.add(poincareButton, c);
		shrinkPanel.add(visualizationPanel, c);
		
		shrinkPanel.add(reduceBtn, c);
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object s = e.getSource();
		if (unwrapBtn == s) {
			if (activeGeometry == null) {
				getGeometry();
			}
			try {
				unwrapActiveGeometry();
			} catch (Exception e1) {
				Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
				JOptionPane.showMessageDialog(w, e1.getMessage(), "Unwrap Error", ERROR_MESSAGE);
				e1.printStackTrace();
				return;
			}
			updateViewer();
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
		if (getGeometryBtn == s) {
			getGeometry();
		}
		if (showUnwrapped == s || kleinButton == s || poincareButton == s) {
			klein = kleinButton.isSelected();
			updateViewer();
		}
		if (showUnitCircle == s) {
			unitCircle.setVisible(showUnitCircle.isSelected());
		}
		if (reduceBtn == s) {
			try {
				reduceToFundamentalPolygon();
			} catch (Exception e1) {
//				Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
//				JOptionPane.showMessageDialog(w, e1.getMessage(), "Reduce Error", ERROR_MESSAGE);
				e1.printStackTrace(); 
				return;
			}
			updateViewer();
		}
	}
	
	private void getGeometry() {
		activeGeometry = new CoHDS();
		unwrappedGeometry = null;
		activeGeometry.setTexCoordinatesValid(false);
		hcp.getHalfedgeContent(activeGeometry, new PositionAdapter());
		triangulator.triangulate(activeGeometry);
		geometryLabel.setText("Geometry: " + hcp.getHalfedgeContentName());
	} 
	
	
	private void updateViewer() {
		if (unwrappedGeometry == null) {
			return;
		}
		Adapter texAdapter = new TexCoordAdapter(false, !klein);
		Adapter posAdapter = null;
		if (showUnwrapped.isSelected()) {
			posAdapter = new PositionTexCoordAdapter(!klein);
		} else {
			posAdapter = new PositionAdapter();
		}
		hcp.updateHalfedgeContent(unwrappedGeometry, true, posAdapter, texAdapter, pointColorAdapter); 
	}
	
	
	private void reduceToFundamentalPolygon() throws Exception {
		unwrappedGeometry = activeGeometry.createCombinatoriallyEquivalentCopy(new CoHDS());
		for (CoVertex v : activeGeometry.getVertices()) {
			unwrappedGeometry.getVertex(v.getIndex()).setPosition(v.getPosition());
		}
		Vector u = null;
		if (usePetsc) {
			CHyperbolicUnwrapperPETSc unwrapper = new CHyperbolicUnwrapperPETSc();
			u = unwrapper.getConformalFactors(unwrappedGeometry);
		} else {
			CHyperbolicUnwrapper unwrapper = new CHyperbolicUnwrapper();
			u = unwrapper.getConformalFactors(unwrappedGeometry);
		}
		UAdapter uAdapter = new UAdapter(u);
		CoVertex root = unwrappedGeometry.getVertex(0);
		UniformizationUtility.reduceToFundamentalPolygon(unwrappedGeometry, root, uAdapter);
	}
	
	
	
	private void unwrapActiveGeometry() throws Exception {
		unwrappedGeometry = activeGeometry.createCombinatoriallyEquivalentCopy(new CoHDS());
		for (CoVertex v : activeGeometry.getVertices()) {
			unwrappedGeometry.getVertex(v.getIndex()).setPosition(v.getPosition());
		}
		if (hyperbolic) {
			unwrappedGeometry.prepareInvariantDataHyperbolic();
			Vector u = null;
			if (usePetsc) {
				CHyperbolicUnwrapperPETSc unwrapper = new CHyperbolicUnwrapperPETSc();
				u = unwrapper.getConformalFactors(unwrappedGeometry);
			} else {
				CHyperbolicUnwrapper unwrapper = new CHyperbolicUnwrapper();
				u = unwrapper.getConformalFactors(unwrappedGeometry);
			}
			HyperbolicLayoutContext context = CHyperbolicLayout.doLayout(unwrappedGeometry, u);
			pointColorAdapter.setContext(context);
		} else {
			CUnwrapper unwrapper = null;
			if (usePetsc) {
				unwrapper = new CDiskUnwrapperPETSc(numCones, quantizeCones);
			} else {
				unwrapper = new CDiskUnwrapper(numCones, quantizeCones);
			}
			unwrapper.unwrap(unwrappedGeometry);
		}
	}
	
	
	
	@Override
	public void stateChanged(ChangeEvent e) {
		if (numConesSpinner == e.getSource()) {
			numCones = numConesModel.getNumber().intValue();
		}
	}
	
	
	@Override
	public void install(Controller c) throws Exception {
		hcp = c.getPlugin(HalfedgeConnectorPlugin.class);
		content = c.getPlugin(AlignedContent.class);
		content.getScalingComponent().addChild(auxGeometry);
		super.install(c); 
	}
	
	@Override
	public void uninstall(Controller c) throws Exception {
		super.uninstall(c);
		if (content.getScalingComponent().isDirectAncestor(auxGeometry)) {
			content.getScalingComponent().removeChild(auxGeometry);
		} 
	}
	
	
	@Override
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "numCones", numCones);
		c.storeProperty(getClass(), "quantizeCones", quantizeCones);
		c.storeProperty(getClass(), "hyperbolic", hyperbolic);
		c.storeProperty(getClass(), "usePetsc", usePetsc);
		c.storeProperty(getClass(), "klein", klein);
	}
	

	@Override
	public void restoreStates(Controller c) throws Exception { 
		super.restoreStates(c);
		numConesModel.setValue(c.getProperty(getClass(), "numCones", numCones));
		quantizeChecker.setSelected(c.getProperty(getClass(), "quantizeCones", quantizeCones));
		hyperbolic = c.getProperty(getClass(), "hyperbolic", hyperbolic);
		usePetsc = c.getProperty(getClass(), "usePetsc", usePetsc);
		klein = c.getProperty(getClass(), "klein", klein);
		
		numericsCombo.setSelectedIndex(usePetsc ? 1 : 0);
		quantizeCones = quantizeChecker.isSelected();
		numCones = numConesModel.getNumber().intValue();
		euclideanButton.setSelected(!hyperbolic);
		hyperbolicButton.setSelected(hyperbolic);
		kleinButton.setSelected(klein);
		poincareButton.setSelected(!klein);
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
