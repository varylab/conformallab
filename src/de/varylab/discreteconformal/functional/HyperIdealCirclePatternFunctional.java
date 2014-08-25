package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.HyperIdealUtility.calculateTetrahedronVolume;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ_13;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ_14;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ_15;
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
		E e12 = f.getBoundaryEdge();
		E e23 = e12.getNextEdge();
		E e31 = e23.getNextEdge();
		V v1 = e12.getStartVertex();
		V v2 = e23.getStartVertex();
		V v3 = e31.getStartVertex();
		boolean e12b = var.isVariable(e12);
		boolean e23b = var.isVariable(e23);
		boolean e31b = var.isVariable(e31);
		boolean v1b = var.isVariable(v1);
		boolean v2b = var.isVariable(v2);
		boolean v3b = var.isVariable(v3);
		double a12 = e12b ? 0.0 : x.get(var.getVarIndex(e12));
		double a23 = e23b ? 0.0 : x.get(var.getVarIndex(e23));
		double a31 = e31b ? 0.0 : x.get(var.getVarIndex(e31));
		double b1 = v1b ? 0.0 : x.get(var.getVarIndex(v1));
		double b2 = v2b ? 0.0 : x.get(var.getVarIndex(v2));
		double b3 = v3b ? 0.0 : x.get(var.getVarIndex(v3));
		double αβ[] = αβ(a12, a23, a31, b1, b2, b3, v1b, v2b, v3b);
		double aa = a12*αβ[0] + a23*αβ[1] + a31*αβ[2];
		double bb = b1*αβ[3] + b2*αβ[4] + b3*αβ[5];
		double V = calculateTetrahedronVolume(αβ[3], αβ[4], αβ[5], αβ[0], αβ[1], αβ[2]);
		return aa + bb + V;
	}
	
	private double[] αβ(
		double a12, double a23, double a31,
		double b1, double b2, double b3,
		boolean v1b, boolean v2b, boolean v3b 
	) {
		double l12 = lij(b1, b2, a12, v1b, v2b);
		double l23 = lij(b2, b3, a23, v2b, v3b);
		double l31 = lij(b3, b1, a31, v3b, v1b);
		double β1 = ζ(l23, l12, l31);
		double β2 = ζ(l31, l12, l23);
		double β3 = ζ(l12, l23, l31);
		return new double[] {0.0, 0.0, 0.0, β1, β2, β3};
	}

	private double lij(double bi, double bj, double aij, boolean vib, boolean vjb) {
		if (vib && vjb) {
			return ζ_13(bi, bj, aij);
		} else if (vib) {
			return ζ_14(aij, bi);
		} else if (vjb) {
			return ζ_14(aij, bj);
		} else {
			return ζ_15(aij);
		}
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
