package de.varylab.discreteconformal.functional;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.varylab.discreteconformal.functional.Clausen.Л;
import static java.lang.Math.PI;
import static java.lang.Math.asin;
import static java.lang.Math.atan2;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import no.uib.cipr.matrix.DenseVector;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;
import de.jtem.numericalMethods.calculus.minimizing.DBrent;
import de.jtem.numericalMethods.calculus.minimizing.Info;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Lambda;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.unwrapper.numerics.MTJGradient;
import de.varylab.discreteconformal.unwrapper.numerics.SimpleEnergy;

public class SphericalFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> implements ConformalFunctional<V, E, F> {

	private Variable<V, E> 
		var = null;
	private Theta<V> 
		theta = null;
	private Lambda<E> 
		lambda = null;
	private Alpha<E> 
		alpha = null;
	private InitialEnergy<F> 
		initE = null;
	private boolean
		reduced = true;
	private double
		logScale = 0.0;
	
	public SphericalFunctional(
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
		if (isReduced()) {
			u = maximizeInNegativeDirection(hds, u);
		}
		evaluate_internal(hds, u, E, G, H);
//		if (G != null && G instanceof MTJGradient) {
//			System.out.print("|G|: " + ((MTJGradient)G).getG().norm(Norm.Two) + " - ");
//			if (E != null) {
//				System.out.println(E.get());
//			} else {
//				System.out.println(" ");
//			}
//		}
//		if (G != null && G instanceof TaoGradient) {
//			System.out.print("|G|: " + ((TaoGradient)G).getG().norm(NormType.NORM_FROBENIUS) + " - ");
//			if (E != null) {
//				System.out.println(E.get());
//			} else {
//				System.out.println(" ");
//			}
//		}
	};
	
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> void evaluate_internal(
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

	private class MaximizingFunctional implements RealFunctionOfOneVariable {

		private DomainValue
			baseU = null;
		private int	
			dim = 0;
		private HalfEdgeDataStructure<V, E, F>
			hds = null;
		
		private MaximizingFunctional(HalfEdgeDataStructure<V, E, F> hds, DomainValue baseU, int dim) {
			this.hds = hds;
			this.baseU = baseU;
			this.dim = dim;
		}

		@Override
		public double eval(double x) {
			DomainValue u = new MTJDomain(new DenseVector(dim));
			for (int i = 0; i < dim; i++) {
				u.add(i, baseU.get(i) + x);
			}
			SimpleEnergy e = new SimpleEnergy();
			evaluate_internal(hds, u, e, null, null);
			return -e.get();
		}

	}
	
	private class MaximizingDerivative implements RealFunctionOfOneVariable {

		private DomainValue
			baseU = null;
		private int	
			dim = 0;
		private HalfEdgeDataStructure<V, E, F>
			hds = null;
		
		private MaximizingDerivative(HalfEdgeDataStructure<V, E, F> hds, DomainValue baseU, int dim) {
			this.hds = hds;
			this.baseU = baseU;
			this.dim = dim;
		}
		
		@Override
		public double eval(double x) {
			DomainValue u = new MTJDomain(new DenseVector(dim));
			for (int i = 0; i < dim; i++) {
				u.add(i, baseU.get(i) + x);
			}
			Gradient G = new MTJGradient(new DenseVector(dim));
			evaluate_internal(hds, u, null, G, null);
			double sum = 0.0;
			for (int i = 0; i < dim; i++) {
				sum += G.get(i);
			}
			return -sum;
		}

	}	

	
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> DomainValue maximizeInNegativeDirection(HDS hds, DomainValue u) {
		double[] x = {0.0, 0.0};
		int dim = getDimension(hds);
		MaximizingFunctional f = new MaximizingFunctional(hds, u, dim);
		MaximizingDerivative df = new MaximizingDerivative(hds, u, dim);
		Info info = new Info(true);
		DBrent.search(-1E5, logScale, 1E5, x, f, df, 1E-8, info);
//		System.out.println("maximization iter: " + info.getCurrentIter());
		DomainValue newU = new MTJDomain(new DenseVector(dim));
		logScale = x[0];
		for (int i = 0; i < dim; i++) {
			newU.add(i, u.get(i) + logScale);
		}
		return newU;
	}
	
	
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
			for (final E e : hds.getPositiveEdges()) {
				if (!var.isVariable(e)) continue;
				int i = var.getVarIndex(e);
				double αk = alpha.getAlpha(e);
				double αl = alpha.getAlpha(e.getOppositeEdge());
				G.add(i, αk + αl - PI);
			}
		}
	}
	
	
	public void conformalHessian(
		// combinatorics
			final HalfEdgeDataStructure<V, E, F> hds,
		// input
			final DomainValue u,
		// output	
			final Hessian H
	) {
		H.setZero();
		// Face Energy
		for(final E e : hds.getPositiveEdges()) {
			double cot = 0.0;
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
					cot += 1 / tan(βk);
				}
			}
			if (e.getRightFace() != null) {
				E el = e.getOppositeEdge();
				boolean valid = triangleEnergyAndAlphas(u, el.getLeftFace(), null, initE);
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
			double lEuc = exp(λij / 2);
			lEuc = lEuc > 1 ? 1 : lEuc;
//			if (lEuc == 1) {
//				System.err.println("warning (hessian): length clamped to PI");
//			}
			final double
				lij = 2*asin(lEuc);
			double
				tan2 = tan(lij / 2);
			tan2 *= tan2;
			final double 
				Hii = 0.5 * cot * (1 - tan2),
				Hij = 0.5 * cot * (-tan2 - 1);
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
			λij = λk + (var.isVariable(eij) ? 0 : ui + uj),
			λjk = λi + (var.isVariable(ejk) ? 0 : uj + uk),
			λki = λj + (var.isVariable(eki) ? 0 : uk + ui);
		double lijEuc = exp(λij / 2);
		double ljkEuc = exp(λjk / 2);
		double lkiEuc = exp(λki / 2);
		lijEuc = lijEuc > 1 ? 1 : lijEuc;
		ljkEuc = ljkEuc > 1 ? 1 : ljkEuc;
		lkiEuc = lkiEuc > 1 ? 1 : lkiEuc;
//		if (lijEuc == 1 || ljkEuc == 1 || lkiEuc == 1) {
//			System.err.println("warning (energy): length clamped to PI");
//		}
		final double
			lij = 2 * asin(lijEuc),
			ljk = 2 * asin(ljkEuc),
			lki = 2 * asin(lkiEuc);
		final double
			Δij = - lij + ljk + lki,
			Δjk = + lij - ljk + lki,
			Δki = + lij + ljk - lki,
			Δijk = + lij + ljk + lki;
		final double
			si = sin(Δjk / 2),
			sj = sin(Δki / 2),
			sk = sin(Δij / 2),
			s = sin(Δijk / 2);
		double
			αi = 0.0,
			αj = 0.0,
			αk = 0.0;
		final boolean 
			valid = Δij > 0 && Δjk > 0 && Δki > 0 && Δijk < 2*PI;
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
		} else if (Δijk >= 2*PI) {
			αi = PI;
			αj = PI;
			αk = PI;
		}
		final double
			βi = 0.5 * (PI + αi - αj - αk),
			βj = 0.5 * (PI - αi + αj - αk),
			βk = 0.5 * (PI - αi - αj + αk);
		if (E != null) {
			E.add(- αi*ui - αj*uj - αk*uk);
			E.add(βi*λi + βj*λj + βk*λk);
			E.add(+ Л(αi) + Л(αj) + Л(αk) + Л(βi) + Л(βj) + Л(βk));
			E.add(+ Л(0.5 * (PI - αi - αj - αk)));
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
		int n = getDimension(hds);
		int[][] nz = new int[n][];
		for (V v : hds.getVertices()) {
			if (!var.isVariable(v)) {
				continue;
			}
			int i = var.getVarIndex(v);
			List<E> star = incomingEdges(v);
			Set<Integer> nzList = new HashSet<Integer>();
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
		double u1 = var.isVariable(v1) ? u.get(i1) : 0.0; 
		double u2 = var.isVariable(v2) ? u.get(i2) : 0.0;
		double l2 = var.isVariable(e) ? u.get(ei) : lambda.getLambda(e) + u1 + u2;
		double lEuc = getLength(l2);
		return 2 * asin(lEuc);
	}
	@Override
	public double getVertexU(V v, DomainValue u) {
		int i = var.getVarIndex(v);
		return var.isVariable(v) ? u.get(i) : 0.0; 
	};
	
	public void setReduced(boolean reduced) {
		this.reduced = reduced;
	}
	public boolean isReduced() {
		return reduced;
	}
	
	public double getLogScale() {
		return logScale;
	}
	
}

