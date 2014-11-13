package de.varylab.discreteconformal.functional;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.jreality.math.Pn.innerProduct;
import static de.jtem.halfedgetools.functional.FunctionalUtils.addRowToHessian;
import static de.jtem.halfedgetools.functional.FunctionalUtils.addVectorToGradient;
import static java.lang.Math.log;

import java.lang.annotation.Annotation;
import java.util.List;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.unwrapper.numerics.MTJGradient;
import de.varylab.discreteconformal.unwrapper.numerics.MTJHessian;
import de.varylab.discreteconformal.unwrapper.numerics.SimpleEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.TaoDomain;
import de.varylab.discreteconformal.unwrapper.numerics.TaoGradient;
import de.varylab.discreteconformal.unwrapper.numerics.TaoHessian;
import de.varylab.discreteconformal.util.SparseUtility;
import de.varylab.mtjoptimization.Optimizable;

public class MobiusCenteringFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	DATAGET extends Annotation
> implements Functional<V, E, F> {

	private AdapterSet
		a = null;
	private Class<DATAGET>
		pos = null;
	private List<V> 	
		include = null; 
	private double[]
	    g1 = new double[4],
		g2 = new double[4];

	public MobiusCenteringFunctional(Class<DATAGET> pos, AdapterSet a) {
		this(null, pos, a);
	}
	
	public MobiusCenteringFunctional(List<V> include, Class<DATAGET> pos, AdapterSet a) {
		this.a = a;
		this.pos = pos;
		this.include = include;
	}
	
	
	@Override
	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void evaluate(
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
		List<V> vertices = include;
		if (vertices == null) {
			vertices = hds.getVertices();
		}
		for (V v : vertices) {
			double[] p = a.getD(pos, v);
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
	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> int getDimension(HDS hds) {
		return 4;
	}

	@Override
	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> int[][] getNonZeroPattern(HDS hds) {
		int[][] nnz = new int[4][4];
		for (int i = 0; i < nnz.length; i++) {
			for (int j = 0; j < nnz[i].length; j++) {
				nnz[i][j] = j;
			}
		}
		return nnz;
	}
	
	
	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Optimizable getOptimizatble(HDS hds) {
		return new MoebiusCenteringOptimizable<HDS>(hds);
	}
	
	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> TaoApplication getTaoApplication(HDS hds) {
		return new MoebiusCenteringApplication<HDS>(hds);
	}
	
	
	public class MoebiusCenteringOptimizable <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> implements Optimizable { 
		
		private MobiusCenteringFunctional<V, E, F, DATAGET>
			mcf = MobiusCenteringFunctional.this;
		private SimpleEnergy 
			E = new SimpleEnergy();
		private HDS
			hds = null;
		
		public MoebiusCenteringOptimizable(HDS hds) {
			this.hds = hds;
		}
		
		@Override
		public Matrix getHessianTemplate() {
			return new DenseMatrix(4, 4);
		}
		
		@Override
		public Integer getDomainDimension() {
			return mcf.getDimension(hds);
		}
		
		@Override
		public Double evaluate(Vector x, Vector gradient, Matrix hessian) {
			MTJDomain u = new MTJDomain(x);
			MTJGradient G = new MTJGradient(gradient);
			MTJHessian H = new MTJHessian(hessian);
			mcf.evaluate(hds, u, E, G, H);
			return E.get();
		}

		@Override
		public Double evaluate(Vector x, Vector gradient) {
			MTJDomain u = new MTJDomain(x);
			MTJGradient G = new MTJGradient(gradient);
			mcf.evaluate(hds, u, E, G, null);
			return E.get();
		}

		@Override
		public Double evaluate(Vector x, Matrix hessian) {
			MTJDomain u = new MTJDomain(x);
			MTJHessian H = new MTJHessian(hessian);
			mcf.evaluate(hds, u, E, null, H);
			return E.get();
		}

		@Override
		public Double evaluate(Vector x) {
			MTJDomain u = new MTJDomain(x);
			mcf.evaluate(hds, u, E, null, null);
			return E.get();
		}
	};
	

	
	public class MoebiusCenteringApplication <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> extends TaoApplication implements TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess { 
		
		private SimpleEnergy 
			E = new SimpleEnergy();
		private HDS
			hds = null;
		
		public MoebiusCenteringApplication(HDS hds) {
			this.hds = hds;
		}

		@Override
		public double evaluateObjectiveAndGradient(Vec x, Vec g) {
			TaoDomain u = new TaoDomain(x);
			TaoGradient G = new TaoGradient(g);
			evaluate(hds, u, E, G, null);
			g.assemble();
			return E.get();
		}

		@Override
		public PreconditionerType evaluateHessian(Vec x, Mat H, Mat Hpre) {
			TaoDomain u = new TaoDomain(x);
			TaoHessian taoHess = new TaoHessian(H);
			evaluate(hds, u, null, null, taoHess);
			H.assemble();
			return PreconditionerType.SAME_NONZERO_PATTERN;
		}
		
		public int getDomainDimension() {
			return getDimension(hds);
		}
		
		public Mat getHessianTemplate() {
			int dim = getDomainDimension();
			int[][] sparceStructure = getNonZeroPattern(hds);
			int[] nonZeros = SparseUtility.getPETScNonZeros(sparceStructure);
			Mat H = Mat.createSeqAIJ(dim, dim, PETSc.PETSC_DEFAULT, nonZeros);
			H.assemble();
			return H;
		}
		
	}

}