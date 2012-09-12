package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;

import java.util.List;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.util.NativePathUtility;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;


public class SphericalNormalizerPETScOld {
	
	private static double
		tolerance = 1E-6;
	private static int
		maxIterations = 400;
	
	static {
		NativePathUtility.set("native");
		NativePathUtility.set("../DiscreteConformalLab/native");
		Tao.Initialize();		
	}
	
	
	
	public static void normalize(CoHDS hds) throws Exception {
	normalize(hds, hds.getVertices());
	}
	
	public static void normalize(CoHDS hds, List<CoVertex> effectiveList) throws Exception {
		double[][] verts = new double[effectiveList.size()][];
		int i = 0;
		for (CoVertex v : effectiveList) {
			verts[i] = v.T;
			i++;
		}
		normalize(verts);
	}
	
	public static double[] normalize(double[][] verts) throws Exception {
		double[] weights = new double[verts.length];
		for (int i = 0; i < verts.length; i++) {
			Pn.dehomogenize(verts[i], verts[i]);
			verts[i][3] = 0.0;//convert to euclidean vector
			weights[i] = Rn.euclideanNorm(verts[i]);
			Rn.normalize(verts[i], verts[i]);
			// convert back to homogeneous point
			verts[i][3] = 1.0;
		}
		CNormalizerOptimizable opt = new CNormalizerOptimizable(verts, weights);
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
			throw new Exception("normalization did not succeed in PolyederNormalizer: NaN");
		}
		if (lengthEuclid(center) >= 1) {
			throw new Exception("normalization did not succeed in PolyederNormalizer: |center| >= 1");
		}
		
		int_normalize(center, verts);
		return weights;
	}

	
	private static double lengthEuclid(Vec v){
		double result = 0.0;
		for (int i = 0; i < v.getSize(); i++)
			result += v.getValue(i)*v.getValue(i);
		return Math.sqrt(result);
	}
	

	private static void int_normalize(Vec center, double[][] verts){
		Vector e1 = new DenseVector(new double[]{1,0,0,0});
		Vector e2 = new DenseVector(new double[]{0,1,0,0});
		Vector e3 = new DenseVector(new double[]{0,0,1,0});
		// Gram Schmidt w.r.t. Lorentz inner product
		Vector a4 = new DenseVector(new double[]{-center.getValue(0), -center.getValue(1), -center.getValue(2), -1});
		a4.scale(1 / length2(a4));
		Vector a3 = new DenseVector(e3); 
		a3.add(ldot(e3, a4), a4).scale(1 / length(a3));
		Vector a2 = new DenseVector(e2);
		a2.add(ldot(e2, a4), a4).add(-ldot(e2, a3), a3).scale(1 / length(a2));
		Vector a1 = new DenseVector(e1); 
		a1.add(ldot(e1, a4), a4).add(-ldot(e1, a3), a3).add(-ldot(e1, a2), a2).scale(1 / length(a1));
		//
		Matrix At = new DenseMatrix(new Vector[]{a1, a2, a3, a4}).transpose();
		Matrix I_l = new DenseMatrix(4,4);
		I_l.set(0,0,1);I_l.set(1,1,1);I_l.set(2,2,1);I_l.set(3,3,-1);
		Matrix A_inv = At.mult(I_l, new DenseMatrix(4,4));
		
		Vector test = new DenseVector(a4);
		test = A_inv.mult(test, new DenseVector(4));
		
		// transform
		for (double[] v : verts){
			Vector v1 = new DenseVector(v);
			Vector newV = A_inv.mult(v1, new DenseVector(4));
			v[0] = newV.get(0);
			v[1] = newV.get(1);
			v[2] = newV.get(2);
			v[3] = 0.0;
			Rn.normalize(v, v);
			v[3] = 1.0;
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
	
	
	protected static class CNormalizerOptimizable extends TaoApplication implements
		TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {
		
		private double[]
			weights = null;
		
		private double[][]
			verts = null;	
		
		public CNormalizerOptimizable(double[][] normalizedVerts, double[] weights){
			this.weights = weights;
			this.verts = normalizedVerts;
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
		

		public Double evaluate(Vec x) {
			double result = 0;
			double l = -lorentzQuadratic(x);
			for (int i = 0; i < verts.length; i++) {
				result += weights[i] * log( -lorentzDot(verts[i], x) / sqrt(l) );
			}
			return result;
		}

		private double lorentzDot(double[] p, Vec x){
			return x.getValue(0)*p[0] + x.getValue(1)*p[1] + x.getValue(2)*p[2] - 1.0;
		}
		
		private double lorentzQuadratic(Vec x){
			return x.getValue(0)*x.getValue(0) + x.getValue(1)*x.getValue(1) + x.getValue(2)*x.getValue(2) - 1.0;
		}
		
		private void makeGradient(Vec x, Vec g){
			g.zeroEntries();
			double[][] vList = verts;
			double l = -lorentzQuadratic(x);
			for (int i = 0; i < 3; i++){
				double xi = x.getValue(i);
				for (int j = 0; j < vList.length; j++) {
					double[] v = vList[j];
					double pi = v[i];
					double dot = -lorentzDot(v, x);
					g.add(i, weights[j] * (-pi/dot + xi/l));
				}
			}
		}
		
		private void makeHessian(Vec x, Mat hess){
			hess.zeroEntries();
			double[][] v = verts; 
			double l = -lorentzQuadratic(x);
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					for (int k = 0; k < v.length; k++) {
						double xi = x.getValue(i);
						double xj = x.getValue(j);
						double pi = v[k][i];
						double pj = v[k][j];
						double d = -lorentzDot(v[k], x);
						double diag = i == j ? 1 : 0;
						hess.add(i, j,weights[k]*(diag/l + 2*xi*xj/(l*l) - pi*pj/(d*d)));
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

