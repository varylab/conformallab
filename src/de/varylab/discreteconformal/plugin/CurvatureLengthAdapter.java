package de.varylab.discreteconformal.plugin;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMax;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.PrincipalCurvatureMax;
import de.jtem.halfedgetools.adapter.type.PrincipalCurvatureMin;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;

@Length
public class CurvatureLengthAdapter extends AbstractAdapter<Double> {

	public CurvatureLengthAdapter() {
		super(Double.class, true, false);
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Double getE(E e, AdapterSet a) {
		// get curvature information on e
		double kMax = a.get(PrincipalCurvatureMax.class, e, Double.class);
		double kMin = a.get(PrincipalCurvatureMin.class, e, Double.class);
		
		double[] eVec = a.getD(EdgeVector.class, e);
		double[] vMax = a.getD(CurvatureFieldMax.class, e);
		double[] vMin = a.getD(CurvatureFieldMin.class, e);
		Rn.normalize(eVec, eVec);
		double tauMax = Math.abs(Rn.innerProduct(eVec, vMax));
		double tauMin = Math.abs(Rn.innerProduct(eVec, vMin));

		double kRatio = Math.abs(kMin*kMax);
		
		// get real length
		double[] s = a.getD(Position3d.class, e.getStartVertex());
		double[] t = a.getD(Position3d.class, e.getTargetVertex());
		double l = Rn.euclideanDistance(s, t) * kRatio;
		return l;
	}
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return Edge.class.isAssignableFrom(nodeClass);
	}

	@Override
	public double getPriority() {
		return 10;
	}
	
}
