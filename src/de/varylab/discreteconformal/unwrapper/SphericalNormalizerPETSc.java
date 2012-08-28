package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;

import java.lang.annotation.Annotation;
import java.util.List;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;

public class SphericalNormalizerPETSc {
	
	private static double
		tolerance = 1E-6;
	private static int
		maxIterations = 400;
	
	static {
		Tao.Initialize();		
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>,
		DATAGET extends Annotation, 
		DATASET extends Annotation
	> void normalize(HDS hds, AdapterSet a, Class<DATAGET> get, Class<DATASET> set) {
		normalize(hds, hds.getVertices(), a, get, set);
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>,
		DATAGET extends Annotation, 
		DATASET extends Annotation
	> void normalize(HDS hds, List<V> include, AdapterSet a, Class<DATAGET> get, Class<DATASET> set) {
		CNormalizerOptimizable<V, E, F, HDS, DATAGET> 
			opt = new CNormalizerOptimizable<V, E, F, HDS, DATAGET>(include, a, get);
		Vec center = new Vec(opt.getDomainDimension());
		opt.setInitialSolutionVec(center);
		Mat H = opt.getHessianTemplate();
		opt.setHessianMat(H, H);
		
		Tao optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(opt);
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setGradientTolerances(tolerance, tolerance, tolerance);
		optimizer.setMaximumIterates(maxIterations);
		optimizer.solve();
		
		if (lengthEuclid(center) == Double.NaN) {
			throw new RuntimeException("normalization did not succeed in PolyederNormalizer: NaN");
		}
		if (lengthEuclid(center) >= 1) {
			throw new RuntimeException("normalization did not succeed in PolyederNormalizer: |center| >= 1");
		}
		int_normalize(center, hds, a, get, set);
	}

	
	private static double lengthEuclid(Vec v){
		double result = 0.0;
		for (int i = 0; i < v.getSize(); i++)
			result += v.getValue(i)*v.getValue(i);
		return Math.sqrt(result);
	}
	

	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>,
		DATAGET extends Annotation, 
		DATASET extends Annotation
	> void int_normalize(Vec center, HDS hds, AdapterSet a, Class<DATAGET> get, Class<DATASET> set){
		Vector e1 = new DenseVector(new double[]{1,0,0,0});
		Vector e2 = new DenseVector(new double[]{0,1,0,0});
		Vector e3 = new DenseVector(new double[]{0,0,1,0});
		
		Vector a4 = new DenseVector(new double[]{-center.getValue(0), -center.getValue(1), -center.getValue(2), -1});
		a4.scale(1 / length2(a4));
		Vector a3 = new DenseVector(e3); 
		a3.add(ldot(e3, a4), a4).scale(1 / length(a3));
		Vector a2 = new DenseVector(e2);
		a2.add(ldot(e2, a4), a4).add(-ldot(e2, a3), a3).scale(1 / length(a2));
		Vector a1 = new DenseVector(e1); 
		a1.add(ldot(e1, a4), a4).add(-ldot(e1, a3), a3).add(-ldot(e1, a2), a2).scale(1 / length(a1));
		
		Matrix At = new DenseMatrix(new Vector[]{a1, a2, a3, a4}).transpose();
		Matrix I_l = new DenseMatrix(4,4);
		I_l.set(0,0,1);I_l.set(1,1,1);I_l.set(2,2,1);I_l.set(3,3,-1);
		Matrix A_inv = At.mult(I_l, new DenseMatrix(4,4));
		
		Vector test = new DenseVector(a4);
		test = A_inv.mult(test, new DenseVector(4));
		
		// transform
		for (V v : hds.getVertices()){
			double[] vt = a.getD(get, v);
			Vector v1 = new DenseVector(vt);
			Vector newV = A_inv.mult(v1, new DenseVector(4));
			double[] newVT = {newV.get(0), newV.get(1), newV.get(2), newV.get(3)};
			a.set(set, v, newVT);
		}
	}
	
	private static double length2(Vector x){
		return Math.sqrt(x.get(3)*x.get(3) - x.get(0)*x.get(0) - x.get(1)*x.get(1) - x.get(2)*x.get(2));
	}
	
	private static double length(Vector x){
		return Math.sqrt(-x.get(3)*x.get(3) + x.get(0)*x.get(0) + x.get(1)*x.get(1) + x.get(2)*x.get(2));
	}
	
	private static double ldot(Vector x, Vector y){
		return x.get(0)*y.get(0) + x.get(1)*y.get(1) + x.get(2)*y.get(2) - x.get(3)*y.get(3);
	}
	
	
	protected static class CNormalizerOptimizable <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>,
		DATAGET extends Annotation
	> extends TaoApplication implements TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {
		
		private List<V>
			vertices = null;
		private AdapterSet
			a = new AdapterSet();
		private Class<DATAGET>
			get = null;
		
		public CNormalizerOptimizable(List<V> vList, AdapterSet a, Class<DATAGET> get){
			this.vertices = vList;
			this.a = a;
			this.get = get;
		}
		
		@Override
		public double evaluateObjectiveAndGradient(Vec x, Vec g) {
			makeGradient(x, g);
			g.assemble();
			return evaluate(x);
		}
		
		@Override
		public PreconditionerType evaluateHessian(Vec x, Mat H, Mat Hpre) {
			makeHessian(x, H);
			H.assemble();
			return PreconditionerType.SAME_NONZERO_PATTERN;
		}

		public Double evaluate(Vec x)	{
			return evaluate(x, 1.0);
		}
		public Double evaluate(Vec x, double weight) {
			double result = 0;
			double l = myLength(x);
			for (V v : vertices) {
				double[] vt = a.getD(get, v);
				Pn.dehomogenize(vt, vt);
				result += log( dot(vt, x) / sqrt(l) );
			}
			return result;
		}

		
		private double dot(double[] p, Vec x){
			return -x.getValue(0)*p[0] - x.getValue(1)*p[1] - x.getValue(2)*p[2] + 1.0;
		}
		
		
		private double myLength(Vec x){
			return 1 - x.getValue(0)*x.getValue(0) - x.getValue(1)*x.getValue(1) - x.getValue(2)*x.getValue(2);
		}
		
		private void makeGradient(Vec x, Vec g){
			g.zeroEntries();
			List<V> vList = vertices; 
			for (int i = 0; i < 3; i++){
				for (V v : vList){
					double[] vt = a.getD(get, v);
					Pn.dehomogenize(vt, vt);
					double pi = 0;
					if (i == 0) pi = vt[0];
					if (i == 1) pi = vt[1];
					if (i == 2) pi = vt[2];
					double xi = x.getValue(i);
					double dot = dot(vt, x);
					double l = myLength(x);
					g.add(i, (-pi/dot + xi/l));
				}
			}
		}
		
		
		private void makeHessian(Vec x, Mat hess){
			hess.zeroEntries();
			List<V> vList = vertices; 
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					for (V v : vList){
						double[] vt = a.getD(get, v);
						Pn.dehomogenize(vt, vt);
						double xi = x.getValue(i);
						double xj = x.getValue(j);
						double pi = 0;
						if (i == 0) pi = vt[0];
						if (i == 1) pi = vt[1];
						if (i == 2) pi = vt[2];
						double pj = 0;
						if (j == 0) pj = vt[0];
						if (j == 1) pj = vt[1];
						if (j == 2) pj = vt[2];
						double d = dot(vt, x);
						double l = myLength(x);
						double diag = i == j ? 1 : 0;
						hess.add(i, j, diag/l + 2*xi*xj/(l*l) - pi*pj/(d*d));
					}
				}	
			}				
		}
		
		public int getDomainDimension() {
			return 3;
		}

		public Mat getHessianTemplate() {
			Mat H = Mat.createSeqDense(3, 3, null);
			H.assemble();
			return H;
		}
		
	}
	
}
