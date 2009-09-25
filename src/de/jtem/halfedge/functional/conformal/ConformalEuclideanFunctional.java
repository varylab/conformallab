package de.jtem.halfedge.functional.conformal;

import static de.jtem.halfedge.functional.conformal.Clausen.lob;
import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static java.lang.Math.exp;
import static java.lang.Math.sqrt;

import java.util.LinkedList;
import java.util.List;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.functional.DomainValue;
import de.jtem.halfedge.functional.Energy;
import de.jtem.halfedge.functional.Functional;
import de.jtem.halfedge.functional.Gradient;
import de.jtem.halfedge.functional.Hessian;
import de.jtem.halfedge.functional.conformal.ConformalAdapters.Alpha;
import de.jtem.halfedge.functional.conformal.ConformalAdapters.InitialEnergy;
import de.jtem.halfedge.functional.conformal.ConformalAdapters.Lambda;
import de.jtem.halfedge.functional.conformal.ConformalAdapters.Theta;
import de.jtem.halfedge.functional.conformal.ConformalAdapters.Variable;

public class ConformalEuclideanFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	X extends DomainValue
> implements Functional<V, E, F, X> {
	
	private Variable<V> 
		var = null;
	private Theta<V> 
		theta = null;
	private Lambda<E> 
		lambda = null;
	private Alpha<E> 
		alpha = null;
	private InitialEnergy<F> 
		energy = null;
	
	
	public ConformalEuclideanFunctional(
		Variable<V> var,
		Theta<V> theta,
		Lambda<E> lambda,
		Alpha<E> alpha,
		InitialEnergy<F> energy
	) {
		this.var = var;
		this.theta = theta;
		this.lambda = lambda;
		this.alpha = alpha;
		this.energy = energy;
	}
	
	
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> void evaluate(
		HDS hds, 
		X u,
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
	
	
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> int getDimension(HDS hds) {
		int dim = 0;
		for (V v : hds.getVertices()) {
			if (var.isVariable(v)) {
				dim++;
			}
		}
		return dim;
	}

	
	public void conformalEnergyAndGradient(
		// combinatorics
			final HalfEdgeDataStructure<V, E, F> hds,
		// input
			final X u,
		// output
			final Energy E, 
			final Gradient G
	) {
		// Vertex Energy
		if (E != null) {
			E.setZero();
		}
		if (G != null) {
			G.setZero();
		}
		for (final V v : hds.getVertices()) {
			if (!var.isVariable(v)) {
				continue;
			}
			int i = var.getVarIndex(v);
			if (E != null) {
				E.add(theta.getTheta(v) * u.get(i));
			}
			if (G != null) {
				G.add(var.getVarIndex(v), theta.getTheta(v));
			}
		}
		// Face Energy
		for (final F t : hds.getFaces()) {
			final E 
				e1 = t.getBoundaryEdge(),
				e2 = e1.getNextEdge(),
				e3 = e1.getPreviousEdge();
			final V 
				v1 = e1.getTargetVertex(),
				v2 = e2.getTargetVertex(),
				v3 = e3.getTargetVertex();
			final int
				v1i = var.getVarIndex(v1),
				v2i = var.getVarIndex(v2),
				v3i = var.getVarIndex(v3);
			triangleEnergyAndAlphas(hds, u, t, E);
			if (G != null) {
				if (var.isVariable(v1)) {
					G.add(v1i, -alpha.getAlpha(e3));
				}
				if (var.isVariable(v2)) {
					G.add(v2i, -alpha.getAlpha(e1));
				}
				if (var.isVariable(v3)) {
					G.add(v3i, -alpha.getAlpha(e2));
				}
			}
		}
	}
	
	
	public void conformalHessian(
		// combinatorics
			final HalfEdgeDataStructure<V, E, F> hds,
		// input
			final X u,
			final Hessian H
	) {
		H.setZero();
		// Face Energy
		for (final F t : hds.getFaces()) {
			final E 
				e1 = t.getBoundaryEdge(),
				e2 = e1.getNextEdge(),
				e3 = e1.getPreviousEdge();
			final V 
				v1 = e1.getTargetVertex(),
				v2 = e2.getTargetVertex(),
				v3 = e3.getTargetVertex();
			final int
				v1i = var.getVarIndex(v1),
				v2i = var.getVarIndex(v2),
				v3i = var.getVarIndex(v3);
			final double[] 
			     cotE = {0, 0, 0},
			     cotV = {0, 0, 0};
			triangleHessian(hds, u, t, cotE, cotV);
			// edge hessian
			if (var.isVariable(v1) && var.isVariable(v3)) {
				H.add(v1i, v3i, -cotE[0]);
				H.add(v3i, v1i, -cotE[0]);
			}
			if (var.isVariable(v2) && var.isVariable(v1)) {
				H.add(v2i, v1i, -cotE[1]);
				H.add(v1i, v2i, -cotE[1]);
			}
			if (var.isVariable(v3) && var.isVariable(v2)) {
				H.add(v2i, v3i, -cotE[2]);
				H.add(v3i, v2i, -cotE[2]);
			}
			// vertex hessian
			if (var.isVariable(v1)) {
				H.add(v1i, v1i, cotV[0]);
			}
			if (var.isVariable(v2)) {
				H.add(v2i, v2i, cotV[1]);
			}
			if (var.isVariable(v3)) {
				H.add(v3i, v3i, cotV[2]);
			}
		}
	}
	
	
	
	public void triangleEnergyAndAlphas(
		// combinatorics	
			final HalfEdgeDataStructure<V, E, F> hds,
		// input	
			final X u, 
			final F f,
		// output
			final Energy E
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
			a1 = 0.0,
			a2 = 0.0,
			a3 = 0.0;
		final double 
			u1 = var.isVariable(v1) ? u.get(var.getVarIndex(v1)) : 0.0,
			u2 = var.isVariable(v2) ? u.get(var.getVarIndex(v2)) : 0.0,
			u3 = var.isVariable(v3) ? u.get(var.getVarIndex(v3)) : 0.0;
		final double 
			umean = (u1+u2+u3)/3;
		final double 
			x12 = lambda.getLambda(e2) + (u1+u2 - 2*umean),
			x23 = lambda.getLambda(e3) + (u2+u3 - 2*umean),
			x31 = lambda.getLambda(e1) + (u3+u1 - 2*umean);
		final double 
			l12 = exp(x12),
			l23 = exp(x23),
			l31 = exp(x31);
		final double 
			t31 = +l12+l23-l31,
			t23 = +l12-l23+l31,
			t12 = -l12+l23+l31;
		if (t31 > 0 && t23 > 0 && t12 > 0) {
			final double 
				l123 = l12 + l23 + l31,
				denom = sqrt(t12 * t23 * t31 * l123);
			a1 = 2 * atan2(t12 * t31, denom);
			a2 = 2 * atan2(t23 * t12, denom);
			a3 = 2 * atan2(t31 * t23, denom);
		} else if (t31 <= 0) {
			a2 = PI;
		} else if (t23 <= 0) {
			a1 = PI;
		} else if (t12 <= 0) {
			a3 = PI;
		}
		if (E != null) {
			E.add(a1*x23 + a2*x31 + a3*x12);
			E.add(lob(a1) + lob(a2) + lob(a3));
			E.add(- PI * umean - energy.getInitialEnergy(f));
		}
		alpha.setAlpha(e1, a2);
		alpha.setAlpha(e2, a3);
		alpha.setAlpha(e3, a1);
	}
	
	
	
	
	public void triangleHessian(
		// combinatorics	
			final HalfEdgeDataStructure<V, E, F> hds,
		// input	
			final X u,
			final F f,
		// output	
			final double[] cotE, 
			final double[] cotV
	) {
		final E
			e1 = f.getBoundaryEdge(),
			e2 = e1.getNextEdge(),
			e3 = e1.getPreviousEdge();
		final V 
			v1 = e1.getTargetVertex(),
			v2 = e2.getTargetVertex(),
			v3 = e3.getTargetVertex();
		final double 
			u1 = var.isVariable(v1) ? u.get(var.getVarIndex(v1)) : 0.0,
			u2 = var.isVariable(v2) ? u.get(var.getVarIndex(v2)) : 0.0,
			u3 = var.isVariable(v3) ? u.get(var.getVarIndex(v3)) : 0.0;
		final double 
			x12 = lambda.getLambda(e2) + u1 + u2,
			x23 = lambda.getLambda(e3) + u2 + u3,
			x31 = lambda.getLambda(e1) + u3 + u1;
		final double 
			xmean = (x12 + x23 + x31) / 3;
		final double 
			l12 = exp(x12 - xmean),
			l23 = exp(x23 - xmean),
			l31 = exp(x31 - xmean);
		final double
			t31 = +l12+l23-l31,
			t23 = +l12-l23+l31,
			t12 = -l12+l23+l31;
		double 
			cot1 = 0.0,
			cot2 = 0.0,
			cot3 = 0.0;
		if (t31 > 0 && t23 > 0 && t12 > 0) {
			final double
				l123 = l12 + l23 + l31,
				denom = sqrt(t12 * t23 * t31 * l123) * 2;
			cot1 = (t23*l123 - t31*t12) / denom;
			cot2 = (t31*l123 - t12*t23) / denom;
			cot3 = (t12*l123 - t23*t31) / denom;
		}
		cotE[0] = cot2;
		cotE[1] = cot3;
		cotE[2] = cot1;
		cotV[0] = cot2 + cot3;
		cotV[1] = cot3 + cot1;
		cotV[2] = cot1 + cot2;
	}
	
	
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> int[][] getNonZeroPattern(HDS hds) {
		int n = getDimension(hds);
		int[][] nz = new int[n][];
		for (V v : hds.getVertices()) {
			if (!var.isVariable(v)) {
				continue;
			}
			int i = var.getVarIndex(v);
			List<E> star = incomingEdges(v);
			List<Integer> nzList = new LinkedList<Integer>();
			nzList.add(var.getVarIndex(v));
			for (E e : star) {
				V sv = e.getOppositeEdge().getTargetVertex();
				if (var.isVariable(sv)) {
					nzList.add(var.getVarIndex(sv));
				}
			}
			nz[i] = new int[nzList.size()];
			int j = 0;
			for (Integer index : nzList) {
				nz[i][j++] = index;
			}
		}
		return nz;
	}
	
}
