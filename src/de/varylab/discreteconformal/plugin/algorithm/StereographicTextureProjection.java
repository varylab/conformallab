package de.varylab.discreteconformal.plugin.algorithm;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition3d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.math.ComplexUtility;

public class StereographicTextureProjection extends AlgorithmPlugin {

	private HalfedgeInterface
		masterInterface = null;
	
	public StereographicTextureProjection(HalfedgeInterface master) {
		this.masterInterface = master;
	}

	@Override
	public String getAlgorithmName() {
		return "Stereographic Projection";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hi) {
		for (V v : hds.getVertices()) {
			double[] pos = a.getD(TexturePosition3d.class, v);
			Complex c = ComplexUtility.stereographic(pos);
			a.set(TexturePosition.class, v, new double[]{c.re, c.im});
		}
		hi.update();
		masterInterface.update();
	}

}
