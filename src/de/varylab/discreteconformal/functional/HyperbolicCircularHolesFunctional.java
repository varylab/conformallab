package de.varylab.discreteconformal.functional;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.varylab.discreteconformal.functional.Clausen.lob;
import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.tanh;

import java.util.LinkedList;
import java.util.List;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Lambda;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;

public class HyperbolicCircularHolesFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> implements Functional<V, E, F> {

	private Variable<V, E> 
		var = null;
	private Theta<V> 
		theta = null;
	private Lambda<E> 
		lambda = null;
	private Alpha<E> 
		alpha = null;
	private InitialEnergy<F> 
		energy = null;
	
	
	public HyperbolicCircularHolesFunctional(
		Variable<V, E> var,
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
				eij = t.getBoundaryEdge(),
				ejk = eij.getNextEdge(),
				eki = eij.getPreviousEdge();
			final V 
				vi = eij.getStartVertex(),
				vj = ejk.getStartVertex(),
				vk = eki.getStartVertex();
			final int
				ivi = var.getVarIndex(vi),
				ivj = var.getVarIndex(vj),
				ivk = var.getVarIndex(vk);
			triangleEnergyAndAlphas(u, t, E);
			final double
				αi = alpha.getAlpha(ejk),
				αj = alpha.getAlpha(eki),
				αk = alpha.getAlpha(eij),
				αij = 0.5 * (PI - αi - αj + αk),
				αjk = 0.5 * (PI - αj - αk + αi),
				αki = 0.5 * (PI - αk - αi + αj);
			if (G != null) {
				if (var.isVariable(vi)) {
					if (var.isVariable(eij) || var.isVariable(eki)) {
						if (!var.isVariable(eij)) {
							G.add(ivi, αij - PI/2);
						}
						if (!var.isVariable(eki)) {
							G.add(ivi, αki - PI/2);
						}
					} else {
						G.add(ivi, -alpha.getAlpha(ejk));
					}
				}
				if (var.isVariable(vj)) {
					if (var.isVariable(eij) || var.isVariable(ejk)) {
						if (!var.isVariable(eij)) {
							G.add(ivj, αij - PI/2);
						}
						if (!var.isVariable(ejk)) {
							G.add(ivj, αjk - PI/2);
						}
					} else {
						G.add(ivj, -alpha.getAlpha(eki));
					}
				}
				if (var.isVariable(vk)) {
					if (var.isVariable(ejk) || var.isVariable(eki)) {
						if (!var.isVariable(ejk)) {
							G.add(ivk, αjk - PI/2);
						}
						if (!var.isVariable(eki)) {
							G.add(ivk, αki - PI/2);
						}
					} else {
						G.add(ivk, -alpha.getAlpha(eij));
					}	
				}
			}
		}
		// Circular Edges Gradient
		if (G != null) {
			for (final E eij : hds.getPositiveEdges()) {
				if (!var.isVariable(eij)) continue;
				E eji = eij.getOppositeEdge();
				E ejl = eij.getNextEdge();
				E eli = ejl.getNextEdge();
				E eik = eji.getNextEdge();
				E ekj = eik.getNextEdge();
				int i = var.getVarIndex(eij);
				double αij = alpha.getAlpha(eij);
				double αji = alpha.getAlpha(eji);
				double αjl = alpha.getAlpha(ejl);
				double αli = alpha.getAlpha(eli);
				double αik = alpha.getAlpha(eik);
				double αkj = alpha.getAlpha(ekj);
				G.add(i, 0.5 * (αij + αji - αjl - αli - αik - αkj));
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
		// Face Energy
		for(final E e : hds.getPositiveEdges()) {
			double cot = 0.0;
			if (e.getLeftFace() != null) {
				E ek = e;
				boolean valid = triangleEnergyAndAlphas(u, ek.getLeftFace(), null); 
				if (valid) {
					final double
						αi = alpha.getAlpha(ek.getPreviousEdge()),
						αj = alpha.getAlpha(ek.getNextEdge()),
						αk = alpha.getAlpha(ek);
					final double
						βk = 0.5 * (PI - αi - αj + αk);
					cot += 1 / tan(βk);
				}
			}
			if (e.getRightFace() != null) {
				E el = e.getOppositeEdge();
				boolean valid = triangleEnergyAndAlphas(u, el.getLeftFace(), null);
				if (valid) {
					final double
						αi = alpha.getAlpha(el.getPreviousEdge()),
						αj = alpha.getAlpha(el.getNextEdge()),
						αl = alpha.getAlpha(el);
					final double
						βl = 0.5 * (PI - αi - αj + αl);
					cot += 1 / tan(βl);
				}
			}
			final V
				vi = e.getStartVertex(),
				vj = e.getTargetVertex();
			final int
				i = var.getVarIndex(vi),
				j = var.getVarIndex(vj);
			final double 
				ui = var.isVariable(vi) ? u.get(var.getVarIndex(vi)) : 0.0,
				uj = var.isVariable(vj) ? u.get(var.getVarIndex(vj)) : 0.0;
			final double 
				λij = lambda.getLambda(e) + ui + uj;
			final double
				lij = 2*arsinh(exp(λij / 2));
			double 
				tan2 = tanh(lij / 2);
			tan2 *= tan2;
			final double 
				Hii = 0.5 * cot * (1 + tan2),
				Hij = 0.5 * cot * (tan2 - 1);
			if (var.isVariable(vi)) {
				H.add(i, i, Hii);
			}
			if (var.isVariable(vj)) {
				H.add(j, j, Hii);
			}
			if (var.isVariable(vi) && var.isVariable(vj)) {
				H.add(i, j, Hij);
				H.add(j, i, Hij);
			}
		}
	}
	
	
	
	public boolean triangleEnergyAndAlphas(
		// input	
			final DomainValue u, 
			final F f, 
		// output
			final Energy E
	) {
		final E 
			eij = f.getBoundaryEdge(),
			ejk = eij.getNextEdge(),
			eki = eij.getPreviousEdge();
		final V 
			vi = eij.getStartVertex(),
			vj = ejk.getStartVertex(),
			vk = eki.getStartVertex();
		final double 
			ui = var.isVariable(vi) ? u.get(var.getVarIndex(vi)) : 0.0,
			uj = var.isVariable(vj) ? u.get(var.getVarIndex(vj)) : 0.0,
			uk = var.isVariable(vk) ? u.get(var.getVarIndex(vk)) : 0.0;
		final double
			λk = var.isVariable(eij) ? u.get(var.getVarIndex(eij)) : lambda.getLambda(eij),
			λi = var.isVariable(ejk) ? u.get(var.getVarIndex(ejk)) : lambda.getLambda(ejk),
			λj = var.isVariable(eki) ? u.get(var.getVarIndex(eki)) : lambda.getLambda(eki);
		final double 
			λij = λk + (var.isVariable(eij) ? 0 : ui + uj),
			λjk = λi + (var.isVariable(ejk) ? 0 : uj + uk),
			λki = λj + (var.isVariable(eki) ? 0 : uk + ui);
		final double
			lij = 2 * arsinh(exp(λij / 2)),
			ljk = 2 * arsinh(exp(λjk / 2)),
			lki = 2 * arsinh(exp(λki / 2));
		final double
			Δij = - lij + ljk + lki,
			Δjk = + lij - ljk + lki,
			Δki = + lij + ljk - lki,
			Δijk = + lij + ljk + lki;
		final double
			si = sinh(Δjk / 2),
			sj = sinh(Δki / 2),
			sk = sinh(Δij / 2),
			s = sinh(Δijk / 2);
		double
			αi = 0.0,
			αj = 0.0,
			αk = 0.0;
		final boolean 
			valid = Δij > 0 && Δjk > 0 && Δki > 0;
		if (valid) {
			αi = 2 * atan2(sqrt(sj * sk), sqrt(si * s));
			αj = 2 * atan2(sqrt(sk * si), sqrt(sj * s));
			αk = 2 * atan2(sqrt(si * sj), sqrt(sk * s));
		} else if (Δij <= 0) {
			αk = PI;
		} else if (Δjk <= 0) {
			αi = PI;
		} else if (Δki <= 0) {
			αj = PI;
		}
		final double
			βi = 0.5 * (PI + αi - αj - αk),
			βj = 0.5 * (PI - αi + αj - αk),
			βk = 0.5 * (PI - αi - αj + αk);
		if (E != null) {
			E.add(βi*λjk + βj*λki + βk*λij);
			E.add(+ lob(αi) + lob(αj) + lob(αk) + lob(βi) + lob(βj) + lob(βk));
			E.add(+ lob(0.5 * (PI - αi - αj - αk)));
			E.add(-0.5 * PI * (λjk + λki + λij));
			E.add(-energy.getInitialEnergy(f));
		}
		if (alpha != null) {
			alpha.setAlpha(eij, αk);
			alpha.setAlpha(ejk, αi);
			alpha.setAlpha(eki, αj);
		}
		return valid;
	}
	
	
	
	private double arsinh(double x) {
		double r = x + sqrt(x*x + 1);
		return log(r);
	}
	
	
	@Override
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
	
	@Override
	public boolean hasHessian() {
		return true;
	}
	
	
}

