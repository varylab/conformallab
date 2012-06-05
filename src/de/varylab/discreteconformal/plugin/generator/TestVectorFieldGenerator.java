package de.varylab.discreteconformal.plugin.generator;

import static java.lang.Math.PI;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JLabel;
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

	private SpinnerNumberModel
		angleModel = new SpinnerNumberModel(0.25, -0.5, 0.5, 0.1);
	private JSpinner
		radiusSpinner = new JSpinner(angleModel);
	private JPanel
		panel = new JPanel();
	
	public TestVectorFieldGenerator() {
		panel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 1.0;
		c.gridwidth = 1;
		panel.add(new JLabel("Center Rotation"), c);
		c.gridwidth = GridBagConstraints.RELATIVE;
		panel.add(radiusSpinner, c);
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 0.0;
		panel.add(new JLabel("Pi"));
	}
	
	
	@VectorField
	@CurvatureFieldMin
	private class TestVectorField extends AbstractAdapter<double[]> {
		
		private double 
			angle = PI/2;
		
		public TestVectorField(double angle) {
			super(double[].class, true, false);
			this.angle = angle;
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
			double phi = l * angle;
			double[] vec = {Math.cos(phi), 0, Math.sin(phi)};
			return vec;
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}
		
		@Override
		public String toString() {
			return "Test Vector Field " + angle/PI + "Pi";
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
		double angle = angleModel.getNumber().doubleValue() * PI;
		hi.addLayerAdapter(new TestVectorField(angle), false);
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
