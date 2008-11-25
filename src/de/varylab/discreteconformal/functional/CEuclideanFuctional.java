package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.math.Lob.lob;
import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static java.lang.Math.exp;
import static java.lang.Math.sqrt;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.varylab.discreteconformal.functional.CAdapters.Alpha;
import de.varylab.discreteconformal.functional.CAdapters.Energy;
import de.varylab.discreteconformal.functional.CAdapters.Gradient;
import de.varylab.discreteconformal.functional.CAdapters.Hessian;
import de.varylab.discreteconformal.functional.CAdapters.Lambda;
import de.varylab.discreteconformal.functional.CAdapters.Theta;
import de.varylab.discreteconformal.functional.CAdapters.U;
import de.varylab.discreteconformal.functional.CAdapters.Variable;

public class CEuclideanFuctional {

	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void conformalEnergyAndGradient(
		// combinatorics
			final HalfEdgeDataStructure<V, E, F> hds,
		// input
			final U<V> u,
		// output
			final double[] E, 
			final Gradient G,
		// adapters
			final Variable<V> var,
			final Theta<V> theta,
			final Lambda<E> lambda,
			final Alpha<E> alpha,
			final Energy<F> energy
	) {
		// Vertex Energy
		if (E != null) {
			E[0] = 0.0;
		}
		if (G != null) {
			G.setZero();
		}
		for (final V v : hds.getVertices()) {
			if (!var.isVariable(v)) {
				continue;
			}
			if (E != null) {
				E[0] += theta.getTheta(v) * u.getU(v);
			}
			if (G != null) {
				G.addGradient(var.getVarIndex(v), theta.getTheta(v));
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
			final double e = triangleEnergyAndAlphas(hds, u, t, var, lambda, alpha, energy);
			if (E != null) {
				E[0] += e;
			}
			if (G != null) {
				if (var.isVariable(v1)) {
					G.addGradient(v1i, -alpha.getAlpha(e3));
				}
				if (var.isVariable(v2)) {
					G.addGradient(v2i, -alpha.getAlpha(e1));
				}
				if (var.isVariable(v3)) {
					G.addGradient(v3i, -alpha.getAlpha(e2));
				}
			}
		}
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void conformalHessian(
		// combinatorics
			final HalfEdgeDataStructure<V, E, F> hds,
		// input
			final U<V> u,
			final Hessian H,
		// adapters
			final Variable<V> var,
			final Lambda<E> lambda
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
			triangleHessian(hds, u, t, cotE, cotV, var, lambda);
			// edge hessian
			if (var.isVariable(v1) && var.isVariable(v3)) {
				H.addHessian(v1i, v3i, -cotE[0]);
				H.addHessian(v3i, v1i, -cotE[0]);
			}
			if (var.isVariable(v2) && var.isVariable(v1)) {
				H.addHessian(v2i, v1i, -cotE[1]);
				H.addHessian(v1i, v2i, -cotE[1]);
			}
			if (var.isVariable(v3) && var.isVariable(v2)) {
				H.addHessian(v2i, v3i, -cotE[2]);
				H.addHessian(v3i, v2i, -cotE[2]);
			}
			// vertex hessian
			if (var.isVariable(v1)) {
				H.addHessian(v1i, v1i, cotV[0]);
			}
			if (var.isVariable(v2)) {
				H.addHessian(v2i, v2i, cotV[1]);
			}
			if (var.isVariable(v3)) {
				H.addHessian(v3i, v3i, cotV[2]);
			}
		}
	}
	
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double triangleEnergyAndAlphas(
		// combinatorics	
			final HalfEdgeDataStructure<V, E, F> hds,
		// input	
			final U<V> u, 
			final F f, 
		// adapters	
			final Variable<V> var,
			final Lambda<E> lambda,
			final Alpha<E> alpha,
			final Energy<F> energy
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
			u1 = var.isVariable(v1) ? u.getU(v1) : 0.0,
			u2 = var.isVariable(v2) ? u.getU(v2) : 0.0,
			u3 = var.isVariable(v3) ? u.getU(v3) : 0.0;
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
		final double 
			E1 = a1*x23 + a2*x31 + a3*x12,
			E2 = lob(a1) + lob(a2) + lob(a3),
			E3 = - PI * umean - energy.getEnergy(f);
		alpha.setAlpha(e1, a2);
		alpha.setAlpha(e2, a3);
		alpha.setAlpha(e3, a1);
		return E1 + E2 + E3; 
	}
	
	
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void triangleHessian(
		// combinatorics	
			final HalfEdgeDataStructure<V, E, F> hds,
		// input	
			final U<V> u,
			final F f,
		// output	
			final double[] cotE, 
			final double[] cotV,
		// adapters
			final Variable<V> var,
			final Lambda<E> lambda
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
			u1 = var.isVariable(v1) ? u.getU(v1) : 0.0,
			u2 = var.isVariable(v2) ? u.getU(v2) : 0.0,
			u3 = var.isVariable(v3) ? u.getU(v3) : 0.0;
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
	
}
