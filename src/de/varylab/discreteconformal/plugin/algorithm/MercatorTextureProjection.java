package de.varylab.discreteconformal.plugin.algorithm;

import java.util.Iterator;
import java.util.Set;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition2d;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition3d;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.math.CP1;
import de.varylab.discreteconformal.math.ComplexUtility;

public class MercatorTextureProjection extends AlgorithmPlugin {

	private HalfedgeInterface
		masterInterface = null;
	
	public MercatorTextureProjection(HalfedgeInterface master) {
		this.masterInterface = master;
	}

	@Override
	public String getAlgorithmName() {
		return "Mercator Projection";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hi) {
		// north and south pole normalization
		Set<V> poles = hi.getSelection().getVertices(hds);
		if (poles.size() == 3) {
			Iterator<V> it = poles.iterator();
			normalize(hds, it.next(), it.next(), it.next(), a);
		}
		// stereographic projection
		for (V v : hds.getVertices()) {
			double[] pos = a.getD(TexturePosition3d.class, v);
			Complex c = ComplexUtility.stereographic(pos);
			a.set(TexturePosition.class, v, new double[]{c.re, c.im});
		}
		// complex logarithm
		for (V v : hds.getVertices()) {
			double[] pos = a.getD(TexturePosition2d.class, v);
			Complex c = new Complex(pos[0], pos[1]);
			c = c.log();
			a.set(TexturePosition.class, v, new double[]{c.re, c.im});
		}
		hi.update();
		masterInterface.update();
	}
	
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void normalize(HDS hds, V n, V s, V mid, AdapterSet a) {
		double[] pn = a.getD(TexturePosition3d.class, n);
		double[] ps = a.getD(TexturePosition3d.class, s);
		double[] pmid = a.getD(TexturePosition3d.class, mid);
		Complex cn = ComplexUtility.stereographic(pn);
		Complex cs = ComplexUtility.stereographic(ps);
		Complex cmid = ComplexUtility.stereographic(pmid);
		Complex[] Tcp = CP1.standardProjectivity(null, cn, cmid, cs);
		double[] Tp = CP1.convertPSL2CToSO31(null, Tcp);
		for (V v : hds.getVertices()) {
			double[] tp = a.getD(TexturePosition4d.class, v);
			Rn.matrixTimesVector(tp, Tp, tp);
			a.set(TexturePosition.class, v, tp);
		}
	}

}
