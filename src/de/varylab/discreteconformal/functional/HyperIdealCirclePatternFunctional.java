package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.HyperIdealUtility.calculateTetrahedronVolume;
import static java.lang.Math.PI;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;

public class HyperIdealCirclePatternFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> implements Functional<V, E, F> {
	
	private Variable<V, E> 
		var = null;
	private Theta<V, E>
		theta = null;

	public HyperIdealCirclePatternFunctional(Variable<V, E> var, Theta<V, E> theta) {
		this.var = var;
		this.theta = theta;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<V, E, F>> void evaluate(
		HDS hds, 
		DomainValue x, 
		Energy E, 
		Gradient G, 
		Hessian H
	) {
		if (E != null) {
			E.setZero();
			for (F f : hds.getFaces()) {
				E.add(U(f, x));
			}
			for (E e : hds.getPositiveEdges()) {
				int ie = var.getVarIndex(e);
				double a = var.isVariable(e) ? 0.0 : x.get(ie);
				double θ = theta.getTheta(e);
				E.add(-θ*a);
			}
			for (V v : hds.getVertices()) {
				int iv = var.getVarIndex(v);
				double b = var.isVariable(v) ? 0.0 : x.get(iv);
				E.add(-2 * PI * b);
			}
		}
		if (G != null) {
			G.setZero();
			
		}
	}
	
	private double U(F f, DomainValue x) {
		E e1 = f.getBoundaryEdge();
		E e2 = e1.getNextEdge();
		E e3 = e2.getNextEdge();
		V v1 = e1.getStartVertex();
		V v2 = e2.getStartVertex();
		V v3 = e3.getStartVertex();
		boolean e1b = var.isVariable(e1);
		boolean e2b = var.isVariable(e2);
		boolean e3b = var.isVariable(e3);
		boolean v1b = var.isVariable(v1);
		boolean v2b = var.isVariable(v2);
		boolean v3b = var.isVariable(v3);
		double a1 = e1b ? 0.0 : x.get(var.getVarIndex(e1));
		double a2 = e2b ? 0.0 : x.get(var.getVarIndex(e2));
		double a3 = e3b ? 0.0 : x.get(var.getVarIndex(e3));
		double b1 = v1b ? 0.0 : x.get(var.getVarIndex(v1));
		double b2 = v2b ? 0.0 : x.get(var.getVarIndex(v2));
		double b3 = v3b ? 0.0 : x.get(var.getVarIndex(v3));
		double αβ[] = αβ(a1, a2, a3, b1, b2, b3, v1b, v2b, v3b);
		double aa = a1*αβ[0] + a2*αβ[1] + a3*αβ[2];
		double bb = b1*αβ[3] + b2*αβ[4] + b3*αβ[5];
		double V = calculateTetrahedronVolume(αβ[3], αβ[4], αβ[5], αβ[0], αβ[1], αβ[2]);
		return aa + bb + V;
	}
	
	private double[] αβ(
		double a1, double a2, double a3,
		double b1, double b2, double b3,
		boolean v1, boolean v2, boolean v3 
	) {
		int index = 1;
		index += v1 ? 0 : 1;
		index += v2 ? 0 : 1;
		index += v3 ? 0 : 1;
		switch (index) {
		case 1:
			
			break;
		case 2:
			
			break;
		case 3:
			
			break;
		case 4:
			
			break;
		}
		return new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	}
	
	@Override
	public boolean hasHessian() {
		return false;
	}

	@Override
	public boolean hasGradient() {
		return true;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<V, E, F>> int getDimension(HDS hds) {
		return hds.numEdges() / 2;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<V, E, F>> int[][] getNonZeroPattern(HDS hds) {
		return null;
	}

}
