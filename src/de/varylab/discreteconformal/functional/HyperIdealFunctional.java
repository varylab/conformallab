package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.HyperIdealUtility.calculateTetrahedronVolume;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.calculateTetrahedronVolumeWithIdealVertexAtGamma;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ_13;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ_14;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ_15;
import static java.lang.Math.PI;

import java.util.logging.Logger;

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

public class HyperIdealFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> implements Functional<V, E, F> {
	
	private Logger
		log = Logger.getLogger(HyperIdealFunctional.class.getName());
	private Variable<V, E> 
		var = null;
	private Theta<V, E>
		theta = null;
	private Alpha<E>
		alpha = null;
	private Beta<E>
		beta = null;

	public HyperIdealFunctional(Variable<V, E> var, Theta<V, E> theta, Alpha<E> alpha, Beta<E> beta) {
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
				double a = var.isVariable(e) ? x.get(ie) : 0.0;
				double θ = theta.getTheta(e);
				E.add(-θ * a);
			}
			for (V v : hds.getVertices()) {
				int iv = var.getVarIndex(v);
				double b = var.isVariable(v) ? x.get(iv) : 0.0;
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
					G.add(index, beta.getBeta(e.getPreviousEdge()));
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
		double a12 = e12b ? x.get(var.getVarIndex(e12)) : 0.0;
		double a23 = e23b ? x.get(var.getVarIndex(e23)) : 0.0;
		double a31 = e31b ? x.get(var.getVarIndex(e31)) : 0.0;
		double b1 = v1b ? x.get(var.getVarIndex(v1)) : 0.0;
		double b2 = v2b ? x.get(var.getVarIndex(v2)) : 0.0;
		double b3 = v3b ? x.get(var.getVarIndex(v3)) : 0.0;
		if (e12b && a12 < 0 && v1b && v2b) {
			log.warning("a12 must not be negative if v1 and v2 are variables");
			a12 = 0.0;
		}
		if (e23b && a23 < 0 && v2b && v3b) {
			log.warning("a23 must not be negative if v2 and v3 are variables");
			a23 = 0.0;
		}
		if (e31b && a31 < 0 && v3b && v1b) {
			log.warning("a31 must not be negative if v3 and v1 are variables");
			a31 = 0.0;
		}
		if (v1b && b1 < 0) {
			log.warning("argument b1 must not be negative");
			b1 = 0.01;
		}
		if (v2b && b2 < 0) {
			log.warning("argument b2 must not be negative");
			b2 = 0.01;	
		}
		if (v3b && b3 < 0) {
			log.warning("argument b3 must not be negative");
			b3 = 0.01;	
		}
		double l12 = lij(b1, b2, a12, v1b, v2b);
		double l23 = lij(b2, b3, a23, v2b, v3b);
		double l31 = lij(b3, b1, a31, v3b, v1b);
		double β1, β2, β3;
		double α12, α23, α31;
		// check triangle inequalities and extend
		if (l12 > l23 + l31) {
			β1 = 0;
			β2 = 0;
			β3 = PI;
			α12 = PI;
			α23 = 0;
			α31 = 0;
//			log.info("triangle inequality violated at face " + f);
		} else if (l23 > l12 + l31) {
			β1 = PI;
			β2 = 0;
			β3 = 0;
			α12 = 0;
			α23 = PI;
			α31 = 0;
//			log.info("triangle inequality violated at face " + f);
		} else if (l31 > l12 + l23) {
			β1 = 0;
			β2 = PI;
			β3 = 0;
			α12 = 0;
			α23 = 0;
			α31 = PI;
//			log.info("triangle inequality violated at face " + f);
		} else {
			β1 = ζ(l12, l31, l23);
			β2 = ζ(l23, l12, l31);
			β3 = ζ(l31, l23, l12);
			α12 = αij(a12, a23, a31, b1, b2, b3, β1, β2, β3, v1b, v2b, v3b);
			α23 = αij(a23, a31, a12, b2, b3, b1, β2, β3, β1, v2b, v3b, v1b);
			α31 = αij(a31, a12, a23, b3, b1, b2, β3, β1, β2, v3b, v1b, v2b);
		}
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
		double a12 = var.isVariable(e12) ? x.get(var.getVarIndex(e12)) : 0.0;
		double a23 = var.isVariable(e23) ? x.get(var.getVarIndex(e23)) : 0.0;
		double a31 = var.isVariable(e31) ? x.get(var.getVarIndex(e31)) : 0.0;
		double b1 = var.isVariable(v1) ? x.get(var.getVarIndex(v1)) : 0.0;
		double b2 = var.isVariable(v2) ? x.get(var.getVarIndex(v2)) : 0.0;
		double b3 = var.isVariable(v3) ? x.get(var.getVarIndex(v3)) : 0.0;
		double β1 = beta.getBeta(e23);
		double β2 = beta.getBeta(e31);
		double β3 = beta.getBeta(e12);
		double α12 = alpha.getAlpha(e12);
		double α23 = alpha.getAlpha(e23);
		double α31 = alpha.getAlpha(e31);
		double aa = a12*α12 + a23*α23 + a31*α31;
		double bb = b1*β1 + b2*β2 + b3*β3;
		double V = 0.0;
		if (var.isVariable(v1) && var.isVariable(v2) && var.isVariable(v3)) {
			V = calculateTetrahedronVolume(β1, β2, β3, α23, α31, α12);
		} else {
			if (!var.isVariable(v1)) {
				V = calculateTetrahedronVolumeWithIdealVertexAtGamma(β1, α31, α12, α23, β2, β3);
			} else
			if (!var.isVariable(v2)) {
				V = calculateTetrahedronVolumeWithIdealVertexAtGamma(β2, α12, α23, α31, β3, β1);
			} else
			if (!var.isVariable(v3)) {
				V = calculateTetrahedronVolumeWithIdealVertexAtGamma(β3, α23, α31, α12, β1, β2);
			}
		}
		return aa + bb + 2*V;
	}
	
	private double αij(
		double aij, double ajk, double aki, 
		double bi, double bj, double bk, 
		double βi, double βj, double βk, 
		boolean vib, boolean vjb, boolean vkb
	) {
		if (vib) { // case 1
			double σi = σi(aij, aki, ajk, vjb, vkb);
			double σij = σij(aij, bi, bj, vjb);
			double σik = σij(aki, bi, bk, vkb);
			return ζ(σi, σij, σik);
		} else if (vjb) {
			double σj = σi(ajk, aij, aki, vkb, vib);
			double σjk = σij(ajk, bj, bk, vkb);
			double σji = σij(aij, bj, bi, vib);
			return ζ(σj, σji, σjk);
		} else if (vkb) {
			double αjk = αij(ajk, aki, aij, bj, bk, bi, βj, βk, βi, vjb, vkb, vib);
			return PI - αjk - βj;
		} else {
			return 0.5 * (PI + βk - βi - βj);
		}
	}
	
	private double σi(
		double aij, double aki, double ajk, 
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
	
	private double σij(double aij, double bi, double bj, boolean vjb) {
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

	public double getEdgeLength(E e, DomainValue x) {
		V vs = e.getStartVertex();
		V vt = e.getTargetVertex();
		boolean eb = var.isVariable(e);
		boolean vsb = var.isVariable(vs);
		boolean vtb = var.isVariable(vt);
		double ast = eb ? x.get(var.getVarIndex(e)) : 0.0;
		double bs = vsb ? x.get(var.getVarIndex(vs)) : 0.0;
		double bt = vtb ? x.get(var.getVarIndex(vt)) : 0.0;
		return lij(bs, bt, ast, vsb, vtb);
	}

}
