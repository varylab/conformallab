package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;

import java.util.List;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.Optimizable;
import de.varylab.mtjoptimization.Optimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class SphericalNormalizer {

	public static void normalize(CoHDS hds) throws NotConvergentException {
		normalize(hds, hds.getVertices());
	}
	
	public static void normalize(CoHDS hds, List<CoVertex> effectiveList) throws NotConvergentException {
		CNormalizerOptimizable opt = new CNormalizerOptimizable(effectiveList);
		Optimizer o = new NewtonOptimizer();
		ArmijoStepController stepController = new ArmijoStepController();
		o.setStepController(stepController);
		
		o.setMaxIterations(1000);
		o.setError(1E-4);
		Vector result = new DenseVector(opt.getDomainDimension());
		o.minimize(result, opt);
		if (lengthEuclid(result) == Double.NaN)
			throw new NotConvergentException("normalization did not succeed in PolyederNormalizer: NaN", -1.0);
		if (lengthEuclid(result) >= 1)
			throw new NotConvergentException("normalization did not succeed in PolyederNormalizer: |center| >= 1", -1.0);
		int_normalize(result, hds);
	}

	
	private static double lengthEuclid(Vector v){
		double result = 0.0;
		for (int i = 0; i < v.size(); i++)
			result += v.get(i)*v.get(i);
		return Math.sqrt(result);
	}
	

	private static void int_normalize(Vector center, CoHDS hds){
		Vector e1 = new DenseVector(new double[]{1,0,0,0});
		Vector e2 = new DenseVector(new double[]{0,1,0,0});
		Vector e3 = new DenseVector(new double[]{0,0,1,0});
		
		Vector a4 = new DenseVector(new double[]{-center.get(0), -center.get(1), -center.get(2), -1});
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
		for (CoVertex v : hds.getVertices()){
			Vector v1 = new DenseVector(v.T);
			Vector newV = A_inv.mult(v1, new DenseVector(4));
			v.T[0] = newV.get(0);
			v.T[1] = newV.get(1);
			v.T[2] = newV.get(2);
			v.T[3] = newV.get(3);
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
	
	
	
	
	protected static class CNormalizerOptimizable implements Optimizable{
		
		private List<CoVertex>
			vertices = null;
		
		public CNormalizerOptimizable(List<CoVertex> vList){
			this.vertices = vList;
		}
		
		@Override
		public Double evaluate(Vector x, Vector gradient, Matrix hessian) {
			makeGradient(x, gradient);
			makeHessian(x, hessian);
			return evaluate(x);
		}

		@Override
		public Double evaluate(Vector x, Vector gradient) {
			makeGradient(x, gradient);
			return evaluate(x);
		}

		@Override
		public Double evaluate(Vector x, Matrix hessian) {
			makeHessian(x, hessian);
			return evaluate(x);
		}

		@Override
		public Double evaluate(Vector x) {
			double result = 0;
			double l = myLength(x);
			for (CoVertex v : vertices) {
				Pn.dehomogenize(v.T, v.T);
				result += log( dot(v.T, x) / sqrt(l) );
			}
			return result;
		}

		
		private double dot(double[] p, Vector x){
			return -x.get(0)*p[0] - x.get(1)*p[1] - x.get(2)*p[2] + 1.0;
		}
		
		
		private double myLength(Vector x){
			return 1 - x.get(0)*x.get(0) - x.get(1)*x.get(1) - x.get(2)*x.get(2);
		}
		
		private void makeGradient(Vector x, Vector g){
			g.zero();
			List<CoVertex> vList = vertices; 
			for (int i = 0; i < 3; i++){
				for (CoVertex v : vList){
					Pn.dehomogenize(v.T, v.T);
					double pi = 0;
					if (i == 0) pi = v.T[0];
					if (i == 1) pi = v.T[1];
					if (i == 2) pi = v.T[2];
					double xi = x.get(i);
					double dot = dot(v.T, x);
					double l = myLength(x);
					
					g.add(i, (-pi/dot + xi/l));
				}
			}
		}
		
		
		private void makeHessian(Vector x, Matrix hess){
			hess.zero();
			List<CoVertex> vList = vertices; 
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					for (CoVertex v : vList){
						Pn.dehomogenize(v.T, v.T);
						double xi = x.get(i);
						double xj = x.get(j);
						double pi = 0;
						if (i == 0) pi = v.T[0];
						if (i == 1) pi = v.T[1];
						if (i == 2) pi = v.T[2];
						double pj = 0;
						if (j == 0) pj = v.T[0];
						if (j == 1) pj = v.T[1];
						if (j == 2) pj = v.T[2];
						
						double d = dot(v.T, x);
						double l = myLength(x);
						double diag = i == j ? 1 : 0;
						hess.add(i, j, diag/l + 2*xi*xj/(l*l) - pi*pj/(d*d));
					}
				}	
			}				
		}
		
		
		@Override
		public Integer getDomainDimension() {
			return 3;
		}

		@Override
		public Matrix getHessianTemplate() {
			return new DenseMatrix(3, 3);
		}
		
		
	}
	
}
