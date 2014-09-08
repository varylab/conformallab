package de.varylab.discreteconformal.functional;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.jreality.math.Pn.innerProduct;
import static de.jtem.halfedgetools.functional.FunctionalUtils.addRowToHessian;
import static de.jtem.halfedgetools.functional.FunctionalUtils.addVectorToGradient;
import static java.lang.Math.log;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.unwrapper.numerics.MTJGradient;
import de.varylab.discreteconformal.unwrapper.numerics.MTJHessian;
import de.varylab.discreteconformal.unwrapper.numerics.SimpleEnergy;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.Optimizable;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;

public class MobiusCenteringFunctional implements Functional<CoVertex, CoEdge, CoFace> {

	private AdapterSet
		aSet = null;
	private double[]
	    g1 = new double[4],
		g2 = new double[4];
	
	public MobiusCenteringFunctional(AdapterSet a) {
		this.aSet = a;
	}
	
	
	public Optimizable getOptimizatble(final CoHDS hds) {
		return new Optimizable() {
			
			@Override
			public Matrix getHessianTemplate() {
				return new DenseMatrix(4, 4);
			}
			
			@Override
			public Integer getDomainDimension() {
				return 4;
			}
			
			@Override
			public Double evaluate(Vector x, Vector gradient, Matrix hessian) {
				MTJDomain u = new MTJDomain(x);
				MTJGradient G = new MTJGradient(gradient);
				MTJHessian H = new MTJHessian(hessian);
				SimpleEnergy E = new SimpleEnergy();
				MobiusCenteringFunctional.this.evaluate(hds, u, E, G, H);
				return E.get();
			}

			@Override
			public Double evaluate(Vector x, Vector gradient) {
				MTJDomain u = new MTJDomain(x);
				MTJGradient G = new MTJGradient(gradient);
				SimpleEnergy E = new SimpleEnergy();
				MobiusCenteringFunctional.this.evaluate(hds, u, E, G, null);
				return E.get();
			}

			@Override
			public Double evaluate(Vector x, Matrix hessian) {
				MTJDomain u = new MTJDomain(x);
				MTJHessian H = new MTJHessian(hessian);
				SimpleEnergy E = new SimpleEnergy();
				MobiusCenteringFunctional.this.evaluate(hds, u, E, null, H);
				return E.get();
			}

			@Override
			public Double evaluate(Vector x) {
				MTJDomain u = new MTJDomain(x);
				SimpleEnergy E = new SimpleEnergy();
				MobiusCenteringFunctional.this.evaluate(hds, u, E, null, null);
				return E.get();
			}
		}; 
	}

	
	public void normalizeVertices(CoHDS hds) {
		Optimizable opt = getOptimizatble(hds);
		Vector x = new DenseVector(new double[] {0,0,0,1});
		NewtonOptimizer min = new NewtonOptimizer();
		min.setMaxIterations(100);
		min.setError(1E-13);
		try {
			min.minimize(x, opt);
		} catch (NotConvergentException e) {
			e.printStackTrace();
			return;
		}
		double xp[] = {x.get(0), x.get(1), x.get(2), x.get(3)};
		double scale = 1 / Math.sqrt(-Pn.innerProduct(xp, xp, Pn.HYPERBOLIC));
		Rn.times(xp, scale, xp);
		System.out.println("sqrt(-<x,x>): " + Pn.norm(xp, Pn.HYPERBOLIC));
		x = new DenseVector(xp);
		
	}
	
	
	@Override
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> void evaluate(
		HDS hds, 
		DomainValue d, 
		Energy E, 
		Gradient G, 
		Hessian H
	) {
		double[] x = {d.get(0), d.get(1), d.get(2), d.get(3)};
		double xx = innerProduct(x, x, HYPERBOLIC);
		if (E != null) {
			E.setZero();
		}
		if (G != null) {
			G.setZero();
		}
		if (H != null) {
			H.setZero();
		}
		for (CoVertex v : hds.getVertices()) {
			double[] p = aSet.getD(Position4d.class, v);
			double xp = innerProduct(p, x, HYPERBOLIC);
			if (E != null) {
				E.add(-log(xp / xx));
				E.add(xx);
			}
			if (G != null) {
				Rn.times(g1, 2/xx, x);
				Rn.times(g2, -1/xp, p);
				g1[3] *= -1; g2[3] *= -1;
				addVectorToGradient(G, 0, g1);
				addVectorToGradient(G, 0, g2);
				Rn.times(g1, 2, x);
				g1[3] *= -1;
				addVectorToGradient(G, 0, g1);
			}
			if (H != null) {
				for (int i = 0; i < 4; i++) {
					double sign = i == 3 ? -1 : 1;
					Rn.times(g1, sign * p[i] / (xp*xp), p);
					Rn.times(g2, -sign * 4*x[i]/(xx*xx), x);
					g1[3] *= -1; g2[3] *= -1;
					addRowToHessian(H, i, g2);
					addRowToHessian(H, i, g1);
					H.add(i, i, sign * 2 / xx);
					H.add(i, i, sign * 2);
				}
			}
		}
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
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> int getDimension(HDS hds) {
		return 4;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> int[][] getNonZeroPattern(HDS hds) {
		int[][] nnz = new int[4][4];
		for (int i = 0; i < nnz.length; i++) {
			for (int j = 0; j < nnz[i].length; j++) {
				nnz[i][j] = j;
			}
		}
		return nnz;
	}

}
