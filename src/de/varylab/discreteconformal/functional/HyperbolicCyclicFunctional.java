package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.Clausen.Л;
import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.tanh;

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
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;
import de.varylab.discreteconformal.math.MathUtility;

public class HyperbolicCyclicFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> implements ConformalFunctional<V, E, F> {

	private Variable<V, E> 
		var = null;
	private Theta<V, E> 
		theta = null;
	private Lambda<E> 
		lambda = null;
	private Alpha<E> 
		alpha = null;
	private InitialEnergy<F> 
		initE = null;
	
	
	public HyperbolicCyclicFunctional(
		Variable<V, E> var,
		Theta<V, E> theta,
		Lambda<E> lambda,
		Alpha<E> alpha,
		InitialEnergy<F> energy
	) {
		this.var = var;
		this.theta = theta;
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
			triangleEnergyAndAlphas(u, t, E, initE);
			final double
				αi = alpha.getAlpha(ejk),
				αj = alpha.getAlpha(eki),
				αk = alpha.getAlpha(eij);
			if (G != null) {
				if (var.isVariable(vi)) {
					G.add(ivi, -αi);
				}
				if (var.isVariable(vj)) {
					G.add(ivj, -αj);
				}
				if (var.isVariable(vk)) {
					G.add(ivk, -αk);
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
			double cot_k = 0.0;
			double cot_m = 0.0;
			if (e.getLeftFace() != null) {
				E ek = e;
				boolean valid = triangleEnergyAndAlphas(u, ek.getLeftFace(), null, initE); 
				if (valid) {
					final double
						αi = alpha.getAlpha(ek.getPreviousEdge()),
						αj = alpha.getAlpha(ek.getNextEdge()),
						αk = alpha.getAlpha(ek);
					final double
						βk = 0.5 * (PI - αi - αj + αk);
					cot += (cot_k = 1 / tan(βk));
				}
			}
			if (e.getRightFace() != null) {
				E em = e.getOppositeEdge();
				boolean valid = triangleEnergyAndAlphas(u, em.getLeftFace(), null, initE);
				if (valid) {
					final double
						αi = alpha.getAlpha(em.getPreviousEdge()),
						αj = alpha.getAlpha(em.getNextEdge()),
						αm = alpha.getAlpha(em);
					final double
						βm = 0.5 * (PI - αi - αj + αm);
					cot += (cot_m = 1 / tan(βm));
				}
			}
			final V
				vi = e.getStartVertex(),
				vj = e.getTargetVertex();
			final int
				i = var.getVarIndex(vi),
				j = var.getVarIndex(vj);
			final E
				ejk = e.getNextEdge(),
				eki = e.getPreviousEdge(),
				eim = e.getOppositeEdge().getNextEdge(),
				emj = e.getOppositeEdge().getPreviousEdge();			
			final int
				ij = var.getVarIndex(e),
				jk = var.getVarIndex(ejk),
				ki = var.getVarIndex(eki),
				im = var.getVarIndex(eim),
				mj = var.getVarIndex(emj);
			final double 
				ui = var.isVariable(vi) ? u.get(i) : 0.0,
				uj = var.isVariable(vj) ? u.get(j) : 0.0;
			final double 
				λij = (var.isVariable(e) ? u.get(ij) : lambda.getLambda(e)) + ui + uj;
			final double
				lij = 2*MathUtility.arsinh(exp(λij / 2));
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
			
			// quadratic lambda terms
			if (var.isVariable(e)) {
				H.add(ij, ij, 0.5*cot*tan2);
			}
			if (var.isVariable(ejk)) {
				H.add(jk, jk, 0.5*cot_k);
			}
			if (var.isVariable(eki)) {
				H.add(ki, ki, 0.5*cot_k);
			}
			if (var.isVariable(eim)) {
				H.add(im, im, 0.5*cot_m);
			}
			if (var.isVariable(emj)) {
				H.add(mj, mj, 0.5*cot_m);
			}
			// mixed lambda terms
			if (var.isVariable(ejk) && var.isVariable(eki)) {
				H.add(jk, ki, -0.5*cot_k);
				H.add(ki, jk, -0.5*cot_k);
			}
			if (var.isVariable(eim) && var.isVariable(emj)) {
				H.add(im, mj, -0.5*cot_m);
				H.add(mj, im, -0.5*cot_m);
			}			
			// mixed lambda and u terms for both sides
			if (var.isVariable(vi) && var.isVariable(e)) {
				H.add(i, ij, 0.5*cot*tan2);
				H.add(ij, i, 0.5*cot*tan2);
			}	
			if (var.isVariable(vj) && var.isVariable(e)) {
				H.add(j, ij, 0.5*cot*tan2);
				H.add(ij, j, 0.5*cot*tan2);
			}
			// mixed lambda and u terms for the k side of edge
			if (var.isVariable(vi) && var.isVariable(eki)) {
				H.add(i, ki, 0.5*cot_k);
				H.add(ki, i, 0.5*cot_k);
			}
			if (var.isVariable(vi) && var.isVariable(ejk)) {
				H.add(i, jk, -0.5*cot_k);
				H.add(jk, i, -0.5*cot_k);
			}
			if (var.isVariable(vj) && var.isVariable(eki)) {
				H.add(j, ki, -0.5*cot_k);
				H.add(ki, j, -0.5*cot_k);
			}
			if (var.isVariable(vj) && var.isVariable(ejk)) {
				H.add(j, jk, 0.5*cot_k);
				H.add(jk, j, 0.5*cot_k);
			}
			// mixed lambda and u terms for the m side of edge
			if (var.isVariable(vi) && var.isVariable(eim)) {
				H.add(i, im, 0.5*cot_m);
				H.add(im, i, 0.5*cot_m);
			}
			if (var.isVariable(vi) && var.isVariable(emj)) {
				H.add(i, mj, -0.5*cot_m);
				H.add(mj, i, -0.5*cot_m);
			}
			if (var.isVariable(vj) && var.isVariable(eim)) {
				H.add(j, im, -0.5*cot_m);
				H.add(im, j, -0.5*cot_m);
			}
			if (var.isVariable(vj) && var.isVariable(emj)) {
				H.add(j, mj, 0.5*cot_m);
				H.add(mj, j, 0.5*cot_m);
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
			λij = λk + ui + uj,
			λjk = λi + uj + uk,
			λki = λj + uk + ui;
		final double
			lij = 2 * MathUtility.arsinh(exp(λij / 2)),
			ljk = 2 * MathUtility.arsinh(exp(λjk / 2)),
			lki = 2 * MathUtility.arsinh(exp(λki / 2));
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
			E.add(+ Л(αi) + Л(αj) + Л(αk) + Л(βi) + Л(βj) + Л(βk));
			E.add(+ Л(0.5 * (PI - αi - αj - αk)));
			E.add(-0.5 * PI * (λjk + λki + λij));
			E.add(-initialEnergy.getInitialEnergy(f));
		}
		if (alpha != null) {
			alpha.setAlpha(eij, αk);
			alpha.setAlpha(ejk, αi);
			alpha.setAlpha(eki, αj);
		}
		return valid;
	}

	@Override
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> int[][] getNonZeroPattern(HDS hds) {
		int n = 0;
		for (V v : hds.getVertices()) {
			if (var.isVariable(v)) {
				n++;
			}
		}
		for (E e : hds.getPositiveEdges()) {
			if (var.isVariable(e)) {
				n++;
			}
		}
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
		double l2 = var.isVariable(e) ? u.get(ei) : lambda.getLambda(e) + u1 + u2;
		return 2 * MathUtility.arsinh( exp(l2 / 2) );
	}
	@Override
	public double getVertexU(V v, DomainValue u) {
		int i = var.getVarIndex(v);
		return var.isVariable(v) ? u.get(i) : 0.0; 
	};
	
}

