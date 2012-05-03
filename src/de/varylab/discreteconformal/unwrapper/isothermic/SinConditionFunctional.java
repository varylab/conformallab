package de.varylab.discreteconformal.unwrapper.isothermic;

import static java.lang.Math.sin;

import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddObjective;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.unwrapper.IsothermicUtility;

public class SinConditionFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	HDS extends HalfEdgeDataStructure<V, E, F>
> extends TaoApplication implements TaoAppAddObjective {

	private HDS
		hds = null;
	private Map<Integer, Integer>
		egdeIndices = null;
	
	public SinConditionFunctional(HDS hds, Map<Integer, Integer> undirectedIndexMap) {
		this.hds = hds;
		this.egdeIndices = undirectedIndexMap;
	}
	
	public void calculateAndSetInitionSolution(AdapterSet a) {
		Vec alpha = new Vec(hds.numEdges() / 2);
		for (E e : hds.getEdges()) {
			int index = egdeIndices.get(e.getIndex());
			double[] N = a.getD(Normal.class, e);
			double[] Kmin = a.getD(CurvatureFieldMin.class, e);
			double[] E = a.getD(EdgeVector.class, e);
			double ae = IsothermicUtility.getSignedAngle(N, Kmin, E);
			alpha.setValue(index, ae, InsertMode.INSERT_VALUES);
		}
		alpha.assemble();
		setInitialSolutionVec(alpha);
	}
	
	@Override
	public double evaluateObjective(Vec x) {
		double E = 0.0;
		for (V v : hds.getVertices()) {
			E += getVertexEnergy(v, x);
		}
		return E;
	}
	
	protected double getVertexEnergy(V v, Vec aVec) {
		if (HalfEdgeUtils.isBoundaryVertex(v)) {
			return 0;
		}
		double sl = 1.0;
		double sr = 1.0;
		for (E eIn : HalfEdgeUtils.incomingEdges(v)) {
			int ki = egdeIndices.get(eIn.getIndex());
			int ij = egdeIndices.get(eIn.getNextEdge().getIndex());
			int jk = egdeIndices.get(eIn.getPreviousEdge().getIndex());
			double alpha_ki = aVec.getValue(ki);
			double alpha_ij = aVec.getValue(ij);
			double alpha_jk = aVec.getValue(jk);
			double alphaLeft = IsothermicUtility.calculateTriangleAngle(alpha_ij, alpha_jk, alpha_ki);
			double alphaRight = IsothermicUtility.calculateTriangleAngle(alpha_jk, alpha_ki, alpha_ij);
			sl *= sin(alphaLeft);
			sr *= sin(alphaRight);
		}
		return (sl - sr) * (sl - sr);
	}
	
	
}
