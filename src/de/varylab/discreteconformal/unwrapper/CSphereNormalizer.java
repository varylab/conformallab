package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import geom3d.Point;

import java.util.List;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.math.optimization.NotConvergentException;
import de.varylab.discreteconformal.math.optimization.Optimizable;
import de.varylab.discreteconformal.math.optimization.Optimizer;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer;
import de.varylab.discreteconformal.math.optimization.stepcontrol.ArmijoStepController;

public class CSphereNormalizer {

	
	public static void normalize(CHDS hds) throws NotConvergentException{
		CNormalizerOptimizable opt = new CNormalizerOptimizable(hds);
		Optimizer o = new NewtonOptimizer();
		ArmijoStepController stepController = new ArmijoStepController();
		o.setStepController(stepController);
		
		o.setMaxIterations(20);
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
	

	private static void int_normalize(Vector center, CHDS hds){
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
		for (CVertex v : hds.getVertices()){
			Point p = v.getTextureCoord();
			Vector v1 = new DenseVector(new double[]{p.x(), p.y(), p.z(), 1.0});
			Vector newV = A_inv.mult(v1, new DenseVector(4));
			p.set(0, newV.get(0) / newV.get(3));
			p.set(1, newV.get(1) / newV.get(3));
			p.set(2, newV.get(2) / newV.get(3));
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
		
		private CHDS
			hds = null;
		
		public CNormalizerOptimizable(CHDS hds){
			this.hds = hds;
		}
		
		public Double evaluate(Vector x, Vector gradient, Matrix hessian) {
			makeGradient(x, gradient);
			makeHessian(x, hessian);
			return evaluate(x);
		}

		public Double evaluate(Vector x, Vector gradient) {
			makeGradient(x, gradient);
			return evaluate(x);
		}

		public Double evaluate(Vector x, Matrix hessian) {
			makeHessian(x, hessian);
			return evaluate(x);
		}

		public Double evaluate(Vector x) {
			double result = 0;
			double l = myLength(x);
			for (CVertex v : hds.getVertices())
				result += log( dot(v.getTextureCoord(), x) / sqrt(l) );
			return result;
		}

		
		private double dot(Point p, Vector x){
			return -x.get(0)*p.x() - x.get(1)*p.y() - x.get(2)*p.z() + 1.0;
		}
		
		
		private double myLength(Vector x){
			return 1 - x.get(0)*x.get(0) - x.get(1)*x.get(1) - x.get(2)*x.get(2);
		}
		
		private void makeGradient(Vector x, Vector g){
			g.zero();
			List<CVertex> vList = hds.getVertices(); 
			for (int i = 0; i < 3; i++){
				for (CVertex v : vList){
					Point p = v.getTextureCoord();
					double pi = 0;
					if (i == 0) pi = p.x();
					if (i == 1) pi = p.y();
					if (i == 2) pi = p.z();
					double xi = x.get(i);
					
					double dot = dot(v.getTextureCoord(), x);
					double l = myLength(x);
					
					g.add(i, (-pi/dot + xi/l));
				}
			}
		}
		
		
		private void makeHessian(Vector x, Matrix hess){
			hess.zero();
			List<CVertex> vList = hds.getVertices(); 
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					for (CVertex v : vList){
						Point p = v.getTextureCoord();
						double xi = x.get(i);
						double xj = x.get(j);
						double pi = 0;
						if (i == 0) pi = p.x();
						if (i == 1) pi = p.y();
						if (i == 2) pi = p.z();
						double pj = 0;
						if (j == 0) pj = p.x();
						if (j == 1) pj = p.y();
						if (j == 2) pj = p.z();
						
						double d = dot(v.getTextureCoord(), x);
						double l = myLength(x);
						double diag = i == j ? 1 : 0;
						hess.add(i, j, diag/l + 2*xi*xj/(l*l) - pi*pj/(d*d));
					}
				}	
			}				
		}
		
		
		public Integer getDomainDimension() {
			return 3;
		}
		
		
	}
	
}
