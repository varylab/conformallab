package de.varylab.discreteconformal.plugin.algorithm;

import java.util.Iterator;
import java.util.Set;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.math.CP1;
import de.varylab.discreteconformal.math.ComplexUtility;

public class Moebius3d extends AlgorithmPlugin {

	public Moebius3d() {
	}

	@Override
	public String getAlgorithmName() {
		return "Moebius Transformation 3D";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds2, AdapterSet a, HalfedgeInterface hi) {
		// north and south pole normalization
		CoHDS hds = hi.get(new CoHDS());
		Set<CoVertex> vertexSel = hi.getSelection().getVertices(hds);
		if (vertexSel.size() == 3) {
			Iterator<CoVertex> vertexIt = vertexSel.iterator();
			normalize(hds, vertexIt.next(), vertexIt.next(), vertexIt.next(), a);
		} else {
			throw new RuntimeException(
				"Please select three vertices to define the north pole, "
				+ "the equator, and the south pole."
			);
		}
		hi.update();
	}
	
	
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void normalize(CoHDS hds, V s, V mid, V n, AdapterSet a) {
		double[] ps = a.getD(Position3d.class, s);
		double[] pmid = a.getD(Position3d.class, mid);
		double[] pn = a.getD(Position3d.class, n);
		Rn.normalize(ps, ps);
		Rn.normalize(pmid, pmid);
		Rn.normalize(pn, pn);
		Complex cs = ComplexUtility.stereographic(ps);
		Complex cmid = ComplexUtility.stereographic(pmid);
		Complex cn = ComplexUtility.stereographic(pn);
		Complex[] Tcp = CP1.standardProjectivity(null, cs, cmid, cn);
		double[] Tp = CP1.convertPSL2CToSO31(null, Tcp);
		for (CoVertex v : hds.getVertices()) {
			double[] tp = a.getD(Position4d.class, v);
			Rn.matrixTimesVector(tp, Tp, tp);
			a.set(Position.class, v, tp);
		}
	}

}
