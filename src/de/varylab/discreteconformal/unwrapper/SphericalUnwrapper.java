package de.varylab.discreteconformal.unwrapper;

import static no.uib.cipr.matrix.Vector.Norm.Two;

import java.util.Arrays;
import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;
import de.jtem.numericalMethods.calculus.minimizing.DBrent;
import de.jtem.numericalMethods.calculus.minimizing.Info;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CSphericalOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.Optimizable;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.StepController;

public class SphericalUnwrapper implements Unwrapper {

	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	private CoVertex
		layoutRoot = null;
	
	@Override
	public void unwrap(CoHDS hds, int g, AdapterSet a) throws Exception {
		double maxLength = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = a.get(Length.class, e, Double.class);
			maxLength = maxLength < l ? l : maxLength;
		}
		double scale = 1.0/2.0/maxLength;
		
		CSphericalOptimizable opt = new CSphericalOptimizable(hds);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), hds, a, zeroU, scale);

		// optimization
		DenseVector u = calculateConformalFactors(opt);

		layoutRoot = hds.getVertex(0);
		SphericalLayout.doLayout(hds, layoutRoot, opt.getFunctional(), u);
	}
	
	private class MaximizingFunctional implements RealFunctionOfOneVariable {

		private Optimizable
			f = null;
		private Vector
			baseX = null;
		
		private MaximizingFunctional(Optimizable f, Vector baseX) {
			this.f = f;
			this.baseX = baseX;
		}

		@Override
		public double eval(double x) {
			double[] dxArr = new double[f.getDomainDimension()];
			Arrays.fill(dxArr, x);
			Vector xEval = new DenseVector(dxArr);
			xEval.add(baseX);
			return -f.evaluate(xEval);
		}

	}
	
	private class MaximizingDerivative implements RealFunctionOfOneVariable {

		private Optimizable
			f = null;
		private Vector
			baseX = null;
		
		private MaximizingDerivative(Optimizable f, Vector baseX) {
			this.f = f;
			this.baseX = baseX;
		}
		
		@Override
		public double eval(double x) {
			double[] dxArr = new double[f.getDomainDimension()];
			Arrays.fill(dxArr, x);
			Vector xEval = new DenseVector(dxArr);
			xEval.add(baseX);
			Vector gEval = new DenseVector(f.getDomainDimension());
			f.evaluate(xEval, gEval);
			double sum = 0.0;
			for (int i = 0; i < f.getDomainDimension(); i++) {
				sum += gEval.get(i);
			}
			return -sum;
		}

	}	
	
	private class MaximizingStepController implements StepController {

		@Override
		public Double step(Vector x, Double value, Vector dx, Optimizable func, Vector grad, Matrix hess) {
			Vector oldX = new DenseVector(x);
			Vector oldGrad = new DenseVector(grad);
			Double oldGradLength = oldGrad.norm(Two);
			
			x.add(dx);
			Double result = func.evaluate(x, grad, hess);
			Double gradLength = grad.norm(Two);
			int counter = 0;
			boolean success = true;
			while (oldGradLength <= gradLength || gradLength.equals(Double.NaN)) {
				dx.scale(0.5);
				x.set(oldX).add(dx);
				result = func.evaluate(x, grad, hess);
				gradLength = grad.norm(Two);
				counter++;
				if (counter == 100) {
					success = false;
//					throw new RuntimeException("No valid step in step controller!");
				}
			}
			if (!success) {
				result = maximize(x, func, grad, hess);
			}
			return result;
		}

		protected Double maximize(Vector x, Optimizable func, Vector grad, Matrix hess) {
			Double result;
			// maximize in negative direction
			Info info = new Info();
			info.setDebug(true);
			double[] xm = {0.0, 0.0};
			MaximizingFunctional f = new MaximizingFunctional(func, x);
			MaximizingDerivative df = new MaximizingDerivative(func, x);
			DBrent.search(-1E5, 0, 1E5, xm, f, df, 1E-8, info);
			System.out.println("dbrent iterations: " + info.getCurrentIter());
			
			double[] dxArr = new double[func.getDomainDimension()];
			Arrays.fill(dxArr, xm[0]);
			Vector dxMinimize = new DenseVector(dxArr);
			x.add(dxMinimize);
			result = func.evaluate(x, grad, hess);
			
			System.out.println("grad length: " + grad.norm(Norm.Two));
			return result;
		}
		
	}

	DenseVector calculateConformalFactors(CSphericalOptimizable opt) throws UnwrapException {
		int n = opt.getDomainDimension();
		DenseVector u = new DenseVector(n);
		Matrix H = opt.getHessianTemplate();
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		MaximizingStepController stepController = new MaximizingStepController();
		optimizer.setStepController(stepController);
		optimizer.setSolver(Solver.CGS);
		optimizer.setError(gradTolerance);
		optimizer.setMaxIterations(maxIterations);
//		Vector G = new DenseVector(opt.getDomainDimension());
//		stepController.maximize(u, opt, G, H);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		System.out.println("u: " + u);
		return u;
	}

	@Override
	public void setGradientTolerance(double tol) {
		gradTolerance = tol;
	}
	@Override
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}

	@Override
	public void setCutRoot(CoVertex root) {
	}

	@Override
	public CuttingInfo<CoVertex, CoEdge, CoFace> getCutInfo() {
		return null;
	}

	@Override
	public Map<CoEdge, Double> getlengthMap() {
		return null;
	}

	@Override
	public CoVertex getLayoutRoot() {
		return layoutRoot;
	}
	
}
