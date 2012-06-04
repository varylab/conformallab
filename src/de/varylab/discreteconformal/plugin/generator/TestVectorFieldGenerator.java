package de.varylab.discreteconformal.plugin.generator;

import static java.lang.Math.PI;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.VectorField;
import de.jtem.halfedgetools.adapter.type.generic.BaryCenter3d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmDialogPlugin;

public class TestVectorFieldGenerator extends AlgorithmDialogPlugin {

	private JCheckBox
		k1Radio = new JCheckBox("K1", true),
		k2Radio = new JCheckBox("K2"),
		nRadio = new JCheckBox("N"),
		onBoundaryChecker = new JCheckBox("On Boundary");
	private JComboBox
		nodeTypeCombo = new JComboBox(new String[] {"Vertices", "Edges", "Faces"});
	private SpinnerNumberModel
		radiusModel = new SpinnerNumberModel(6.0, 0.1, 10.0, 0.1);
	private JSpinner
		radiusSpinner = new JSpinner(radiusModel);
	private JPanel
		panel = new JPanel(),
		vecPanel = new JPanel();
	
	public TestVectorFieldGenerator() {
//		panel.setLayout(new GridBagLayout());
//		GridBagConstraints c = new GridBagConstraints();
//		c.fill = GridBagConstraints.BOTH;
//		c.weightx = 1.0;
//		c.gridwidth = GridBagConstraints.RELATIVE;
//		panel.add(new JLabel("Radius"), c);
//		c.gridwidth = GridBagConstraints.REMAINDER;
//		panel.add(radiusSpinner, c);
//		panel.add(vecPanel, c);
//		c.gridwidth = GridBagConstraints.RELATIVE;
//		panel.add(new JLabel("On"), c);
//		c.gridwidth = GridBagConstraints.REMAINDER;
//		panel.add(nodeTypeCombo, c);
//		panel.add(onBoundaryChecker, c);
//		nodeTypeCombo.setSelectedIndex(0);
//		
//		vecPanel.setLayout(new GridLayout(1, 3));
//		vecPanel.add(k1Radio);
//		vecPanel.add(k2Radio);
//		vecPanel.add(nRadio);
	}
	
	
	@VectorField
	@CurvatureFieldMin
	private class TestVectorField extends AbstractAdapter<double[]> {
		
		public TestVectorField() {
			super(double[].class, true, false);
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> double[] getE(E e, AdapterSet a) {
			double[] c = a.getD(BaryCenter3d.class, e);
			double r = Rn.maxNorm(new double[] {c[0], c[2]});
			double l = 0.5 + Math.cos(r * PI)/2.0;
			double phi = l * PI/2.0;
			double[] vec = {Math.cos(phi), 0, Math.sin(phi)};
			return vec;
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}
		
	}
	
	
	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void executeAfterDialog (
		HDS hds,
		AdapterSet a,
		HalfedgeInterface hi
	){
		hi.addLayerAdapter(new TestVectorField(), false);
	}
	
	@Override
	protected JPanel getDialogPanel() {
		return panel;
	}

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.VectorField;
	}

	@Override
	public String getAlgorithmName() {
		return "Test Vector Field";
	}

}
