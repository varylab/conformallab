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
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Beta;
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
	private Alpha<E>
		alpha = null;
	private Beta<E>
		beta = null;

	public HyperIdealCirclePatternFunctional(Variable<V, E> var, Theta<V, E> theta, Alpha<E> alpha, Beta<E> beta) {
		this.var = var;
		this.theta = theta;
		this.alpha = alpha;
		this.beta = beta;
	}

	@Override
	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void evaluate(
		HDS hds, 
		DomainValue x, 
		Energy E, 
		Gradient G, 
		Hessian H
	) {
		// calculate angles
		for (F f : hds.getFaces()) {
			calculateAndSetAngles(f, x);
		}
		// evaluate functional
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
				double Θ = theta.getTheta(v);
				E.add(-Θ * b);
			}
		}
		if (G != null) {
			G.setZero();
			for (V v : hds.getVertices()) {
				if (!var.isVariable(v)) continue;
				int index = var.getVarIndex(v);
				for (E e : HalfEdgeUtils.incomingEdges(v)) {
					G.add(index, beta.getBeta(e.getNextEdge()));
				}
				G.add(index, -theta.getTheta(v));
			}
			for (E e : hds.getPositiveEdges()) {
				if (!var.isVariable(e)) continue;
				int index = var.getVarIndex(e);
				G.add(index, alpha.getAlpha(e));
				G.add(index, alpha.getAlpha(e.getOppositeEdge()));
				G.add(index, -theta.getTheta(e));
			}
		}
	}
	
	private void calculateAndSetAngles(F f, DomainValue x) {
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
		double l12 = lij(b1, b2, a12, v1b, v2b);
		double l23 = lij(b2, b3, a23, v2b, v3b);
		double l31 = lij(b3, b1, a31, v3b, v1b);
		double β1 = ζ(l23, l12, l31);
		double β2 = ζ(l31, l12, l23);
		double β3 = ζ(l12, l23, l31);
		double α12 = αij(a12, a23, a31, b1, b2, b3, β1, β2, β3, v1b, v2b, v3b);
		double α23 = αij(a23, a31, a12, b2, b3, b1, β2, β3, β1, v2b, v3b, v1b);
		double α31 = αij(a31, a12, a23, b3, b1, b2, β3, β1, β2, v3b, v1b, v2b);
		alpha.setAlpha(e12, α12);
		alpha.setAlpha(e23, α23);
		alpha.setAlpha(e31, α31);
		beta.setBeta(e23, β1);
		beta.setBeta(e31, β2);
		beta.setBeta(e12, β3);
	}
	
	private double U(F f, DomainValue x) {
		E e12 = f.getBoundaryEdge();
		E e23 = e12.getNextEdge();
		E e31 = e23.getNextEdge();
		V v1 = e12.getStartVertex();
		V v2 = e23.getStartVertex();
		V v3 = e31.getStartVertex();
		double a12 = var.isVariable(e12) ? 0.0 : x.get(var.getVarIndex(e12));
		double a23 = var.isVariable(e23) ? 0.0 : x.get(var.getVarIndex(e23));
		double a31 = var.isVariable(e31) ? 0.0 : x.get(var.getVarIndex(e31));
		double b1 = var.isVariable(v1) ? 0.0 : x.get(var.getVarIndex(v1));
		double b2 = var.isVariable(v2) ? 0.0 : x.get(var.getVarIndex(v2));
		double b3 = var.isVariable(v3) ? 0.0 : x.get(var.getVarIndex(v3));
		double β1 = beta.getBeta(e23);
		double β2 = beta.getBeta(e31);
		double β3 = beta.getBeta(e12);
		double α12 = alpha.getAlpha(e12);
		double α23 = alpha.getAlpha(e23);
		double α31 = alpha.getAlpha(e31);
		double aa = a12*α12 + a23*α23 + a31*α31;
		double bb = b1*β1 + b2*β2 + b3*β3;
		double V = calculateTetrahedronVolume(β1, β2, β3, α12, α23, α31);
		return aa + bb + V;
	}
	
	public double αij(
		double aij, double ajk, double aki, 
		double bi, double bj, double bk, 
		double βi, double βj, double βk, 
		boolean vib, boolean vjb, boolean vkb
	) {
		if (vib) { // case 1
			double σi = σi(aij, ajk, aki, vib, vjb);
			double σij = σij(bi, bj, aij, vjb);
			double σik = σij(bi, bk, aki, vkb);
			return ζ(σi, σij, σik);
		} else if (vjb) {
			double σj = σi(ajk, aki, aij, vjb, vkb);
			double σjk = σij(bj, bk, ajk, vkb);
			double σji = σij(bj, bi, aij, vib);
			return ζ(σj, σjk, σji);
		} else if (vkb) {
			double αjk = αij(ajk, aki, aij, bj, bk, bi, βj, βk, βi, vjb, vkb, vib);
			return PI - αjk - βj;
		} else {
			return 0.5 * (PI + βk - βi - βj);
		}
	}
	
	private double σi(
		double aij, double ajk, double aki, 
		boolean vjb, boolean vkb
	) {
		if (vjb && vkb) {
			return ζ_13(aij, aki, ajk);
		} else if (vjb) {
			return ζ_14(ajk - aki, aij);
		} else if (vkb) {
			return ζ_14(ajk - aij, aki);
		} else {
			return ζ_15(ajk - aij - aki);
		}
	}
	
	private double σij(double bi, double bj, double aij, boolean vjb) {
		if (vjb) {
			return ζ_13(aij, bi, bj);
		} else {
			return ζ_14(-aij, bi);
		}
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
		int dim = 0;
		for (V v : hds.getVertices()) {
			if (var.isVariable(v)) {
				dim++;
			}
		}
		for (E e : hds.getPositiveEdges()) {
			if (var.isVariable(e)) {
				dim++;
			}
		}
		return dim;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<V, E, F>> int[][] getNonZeroPattern(HDS hds) {
		return null;
	}

}
