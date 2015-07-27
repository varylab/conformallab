package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.Clausen.Л;
import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Lambda;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Phi;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;

public class EuclideanCyclicFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> implements ConformalFunctional<V, E, F> {
	
	private Variable<V, E> 
		var = null;
	private Theta<V, E> 
		theta = null;
	private Phi<E>
		phi = null;
	private Lambda<E> 
		lambda = null;
	private Alpha<E> 
		alpha = null;
	private InitialEnergy<F> 
		initE = null;
	
	
	public EuclideanCyclicFunctional(
		Variable<V, E> var,
		Theta<V, E> theta,
		Phi<E> phi,
		Lambda<E> lambda,
		Alpha<E> alpha,
		InitialEnergy<F> energy
	) {
		this.var = var;
		this.theta = theta;
		this.phi = phi;
		this.lambda = lambda;
		this.alpha = alpha;
		this.initE = energy;
	}
	
	
	@Override
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> void evaluate(
		HDS hds, 
		DomainValue u,
		Energy E, 
		Gradient G, 
		Hessian H
	) {
		if (E != null || G != null) {
			conformalEnergyAndGradient(hds, u, E, G);
		}
		if (H != null) {
			conformalHessian(hds, u, H);
		}
	};
	
	
	@Override
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> int getDimension(HDS hds) {
		int dim = 0;
		for (V v : hds.getVertices()) {
			if (var.isVariable(v)) dim++;
		}
		for (E e : hds.getPositiveEdges()) {
			if (var.isVariable(e)) dim++;
		}
		return dim;
	}

	
	public void conformalEnergyAndGradient(
		// combinatorics
			final HalfEdgeDataStructure<V, E, F> hds,
		// input
			final DomainValue u,
		// output
			final Energy E, 
			final Gradient G
	) {
		// Vertex Energy
		if (E != null) E.setZero();
		if (G != null) G.setZero();
		for (final V v : hds.getVertices()) {
			if (!var.isVariable(v)) continue;
			int i = var.getVarIndex(v);
			if (E != null) E.add(theta.getTheta(v) * u.get(i));
			if (G != null) G.add(var.getVarIndex(v), theta.getTheta(v));
		}
		// Face Energy
		for (final F t : hds.getFaces()) {
			final E 
				eij = t.getBoundaryEdge(),
				ejk = eij.getNextEdge(),
				eki = eij.getPreviousEdge();
			final V 
				vi = eki.getTargetVertex(),
				vj = eij.getTargetVertex(),
				vk = ejk.getTargetVertex();
			final int
				i = var.getVarIndex(vi),
				j = var.getVarIndex(vj),
				k = var.getVarIndex(vk);
			final int
				ij = var.getVarIndex(eij),
				jk = var.getVarIndex(ejk),
				ki = var.getVarIndex(eki);
			triangleEnergyAndAlphas(u, t, E, initE);
			final double
				αi = alpha.getAlpha(ejk),
				αj = alpha.getAlpha(eki),
				αk = alpha.getAlpha(eij);
			if (G != null) {
				if (var.isVariable(vi)) G.add(i, -αi);
				if (var.isVariable(vj)) G.add(j, -αj);
				if (var.isVariable(vk)) G.add(k, -αk);
				// circular edges
				if (var.isVariable(eij)) {
					G.add(ij, αk);
					if (eij.isPositive()) G.add(ij, -phi.getPhi(eij));
				}
				if (var.isVariable(ejk)) {
					G.add(jk, αi);
					if (ejk.isPositive()) G.add(jk, -phi.getPhi(ejk));
				}
				if (var.isVariable(eki)) {
					G.add(ki, αj);
					if (eki.isPositive()) G.add(ki, -phi.getPhi(eki));
				}
				
			}
		}
	}
	
	
	public void conformalHessian(
		// combinatorics
			final HalfEdgeDataStructure<V, E, F> hds,
		// input
			final DomainValue u,
			final Hessian H
	) {
		H.setZero();
		for (final E eij : hds.getEdges()) {
			if (eij.getLeftFace() == null) continue;
			final E 
				ejk = eij.getNextEdge(),
				eki = eij.getPreviousEdge();
			final V 
				vi = eij.getTargetVertex(),
				vj = eij.getStartVertex();
			final int
				i = var.getVarIndex(vi),
				j = var.getVarIndex(vj),
				jk = var.getVarIndex(ejk),
				ki = var.getVarIndex(eki);
			if (!triangleEnergyAndAlphas(u, eij.getLeftFace(), null, initE)) continue;
			final double
				αk = alpha.getAlpha(eij);
			final double 
				cot = 1 / tan(αk);			
			final double 
				w = 0.5 * cot;
			// edge hessian
			if (var.isVariable(vi)) {
				H.add(i, i, w);
			}
			if (var.isVariable(vj)) {
				H.add(j, j, w);
			}
			if (var.isVariable(vi) && var.isVariable(vj)) {
				H.add(i, j, -w);
				H.add(j, i, -w);
			}
			// quadratic lambda terms
			if (var.isVariable(ejk)) {
				H.add(jk, jk, w);
			}
			if (var.isVariable(eki)) {
				H.add(ki, ki, w);
			}
			// mixed lambda terms
			if (var.isVariable(ejk) && var.isVariable(eki)) {
				H.add(jk, ki, -w);
				H.add(ki, jk, -w);
			}
			// mixed lambda and u terms
			if (var.isVariable(vi) && var.isVariable(eki)) {
				H.add(i, ki, -w);
				H.add(ki, i, -w);
			}
			if (var.isVariable(vi) && var.isVariable(ejk)) {
				H.add(i, jk, w);
				H.add(jk, i, w);
			}
			if (var.isVariable(vj) && var.isVariable(eki)) {
				H.add(j, ki, w);
				H.add(ki, j, w);
			}
			if (var.isVariable(vj) && var.isVariable(ejk)) {
				H.add(j, jk, -w);
				H.add(jk, j, -w);
			}
		}
	}
	
	@Override
	public boolean triangleEnergyAndAlphas(
		// input	
			final DomainValue u, 
			final F f,
		// output
			final Energy E,
			final InitialEnergy<F> initialEnergy
	) {
		final E 
			e1 = f.getBoundaryEdge(),
			e2 = e1.getNextEdge(),
			e3 = e1.getPreviousEdge();
		final V 
			v1 = e1.getTargetVertex(),
			v2 = e2.getTargetVertex(),
			v3 = e3.getTargetVertex();
		double 
			α1 = 0.0,
			α2 = 0.0,
			α3 = 0.0;
		final double 
			u1 = var.isVariable(v1) ? u.get(var.getVarIndex(v1)) : 0.0,
			u2 = var.isVariable(v2) ? u.get(var.getVarIndex(v2)) : 0.0,
			u3 = var.isVariable(v3) ? u.get(var.getVarIndex(v3)) : 0.0,
			umean = (u1 + u2 + u3) / 3;
		final double
			λ1 = var.isVariable(e1) ? u.get(var.getVarIndex(e1)) : lambda.getLambda(e1),
			λ2 = var.isVariable(e2) ? u.get(var.getVarIndex(e2)) : lambda.getLambda(e2),
			λ3 = var.isVariable(e3) ? u.get(var.getVarIndex(e3)) : lambda.getLambda(e3);
		final double
			Φ1 = phi.getPhi(e1),
			Φ2 = phi.getPhi(e2),
			Φ3 = phi.getPhi(e3);
		final double 
			λ̃1 = λ2 + u1 + u2,
			λ̃2 = λ3 + u2 + u3,
			λ̃3 = λ1 + u3 + u1;
		final double 
			x12 = λ2 + u1 + u2 - 2*umean,
			x23 = λ3 + u2 + u3 - 2*umean,
			x31 = λ1 + u3 + u1 - 2*umean;
		final double 
			l12 = exp(x12/2),
			l23 = exp(x23/2),
			l31 = exp(x31/2);
		final double 
			t31 = +l12+l23-l31,
			t23 = +l12-l23+l31,
			t12 = -l12+l23+l31;
		boolean valid = t31 > 0 && t23 > 0 && t12 > 0;
		if (valid) {
			final double 
				l123 = l12 + l23 + l31,
				denom = sqrt(t12 * t23 * t31 * l123);
			α1 = 2 * atan2(t12 * t31, denom);
			α2 = 2 * atan2(t23 * t12, denom);
			α3 = 2 * atan2(t31 * t23, denom);
		} else if (t31 <= 0) {
			α2 = PI;
		} else if (t23 <= 0) {
			α1 = PI;
		} else if (t12 <= 0) {
			α3 = PI;
		}
		if (E != null) {
			E.add(α1*λ̃2 + α2*λ̃3 + α3*λ̃1);
			E.add(2*Л(α1) + 2*Л(α2) + 2*Л(α3));
			E.add(-(Φ1*λ1 + Φ2*λ2 + Φ3*λ3) / 2 - PI * (u1 + u2 + u3));
			E.add(-initialEnergy.getInitialEnergy(f));
		}
		alpha.setAlpha(e1, α2);
		alpha.setAlpha(e2, α3);
		alpha.setAlpha(e3, α1);
		return valid;
	}
	
	
	@Override
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> int[][] getNonZeroPattern(HDS hds) {
		int n = getDimension(hds);
		Map<Integer, TreeSet<Integer>> nonZeros = new HashMap<Integer, TreeSet<Integer>>();
		for (int i = 0; i < n; i++) {
			nonZeros.put(i, new TreeSet<Integer>());
		}
		for (V v : hds.getVertices()) {
			if (!var.isVariable(v)) continue;
			List<E> star = HalfEdgeUtils.incomingEdges(v);
			Set<Integer> nonZeroIndices = nonZeros.get(var.getVarIndex(v));
			nonZeroIndices.add(var.getVarIndex(v));
			for (E e : star) {
				V connectedVertex = e.getOppositeEdge().getTargetVertex();
				if (var.isVariable(connectedVertex)) {
					nonZeroIndices.add(var.getVarIndex(connectedVertex));
				}
				if (var.isVariable(e)) {
					nonZeroIndices.add(var.getVarIndex(e));
				}
				if (var.isVariable(e.getPreviousEdge())) {
					nonZeroIndices.add(var.getVarIndex(e.getPreviousEdge()));
				}
			}
		}
		for (E e : hds.getEdges()) {
			if (!var.isVariable(e)) continue;
			Set<Integer> nonZeroIndices = nonZeros.get(var.getVarIndex(e));
			
			// quadratic derivative
			nonZeroIndices.add(var.getVarIndex(e));
			
			// mixed edge derivatives
			if (var.isVariable(e.getNextEdge())) {
				nonZeroIndices.add(var.getVarIndex(e.getNextEdge()));
			}
			if (var.isVariable(e.getPreviousEdge())) {
				nonZeroIndices.add(var.getVarIndex(e.getPreviousEdge()));
			}
			
			// mixed vertex derivatives
			if (var.isVariable(e.getTargetVertex())) {
				nonZeroIndices.add(var.getVarIndex(e.getTargetVertex()));
			}
			if (var.isVariable(e.getNextEdge().getTargetVertex())) {
				nonZeroIndices.add(var.getVarIndex(e.getNextEdge().getTargetVertex()));
			}
		}
		int[][] nz = new int[n][];
		for (int j = 0; j < n; j++) {
			Set<Integer> nonZeroIndices = nonZeros.get(j);
			nz[j] = new int[nonZeroIndices.size()];
			int i = 0;
			for (Integer index : nonZeroIndices) {
				nz[j][i++] = index;
			}
		}
		return nz;

	}
	
	@Override
	public boolean hasGradient() {
		return true;
	}
	
	@Override
	public boolean hasHessian() {
		return true;
	}
	
	@Override
	public double getLambda(double length) {
		return 2 * log(length);
	}
	@Override
	public double getLength(double lambda) {
		return exp(lambda / 2);
	}
	
	@Override
	public double getNewLength(E e, DomainValue u) {
		V v1 = e.getStartVertex();
		V v2 = e.getTargetVertex();
		int i1 = var.getVarIndex(v1);
		int i2 = var.getVarIndex(v2);
		int ei = var.getVarIndex(e);
		Double u1 = var.isVariable(v1) ? u.get(i1) : 0.0; 
		Double u2 = var.isVariable(v2) ? u.get(i2) : 0.0;
		double l2 = (var.isVariable(e) ? u.get(ei) : lambda.getLambda(e)) + u1 + u2;
		return getLength(l2);
	}
	@Override
	public double getVertexU(V v, DomainValue u) {
		int i = var.getVarIndex(v);
		return var.isVariable(v) ? u.get(i) : 0.0; 
	};
	
	
}