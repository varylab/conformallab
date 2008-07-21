package de.varylab.discreteconformal.math.optimization.newton;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.QMR;
import de.varylab.discreteconformal.math.optimization.FunctionNotDefinedException;
import de.varylab.discreteconformal.math.optimization.IterationMonitor;
import de.varylab.discreteconformal.math.optimization.Linearizable;
import de.varylab.discreteconformal.math.optimization.NotConvergentException;
import de.varylab.discreteconformal.math.optimization.Solver;



/**
 * This class implements a Newton solver
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 */
public class NewtonSolver implements Solver{

	private Integer
		maxIterations = 20;
	private Double
		error = 1E-8;
	private IterationMonitor
		iterMon = null;

	
	public void solve(Linearizable fun, Vector x, Vector b) throws FunctionNotDefinedException, NotConvergentException{
		Integer iteration = 0;
		Double actError = Double.MAX_VALUE;
		Vector tmp = new DenseVector(fun.getDomainDimension());
		Vector h = new DenseVector(fun.getDomainDimension());
		Vector fx = new DenseVector(fun.getCoDomainDimension());
		Vector nextFx = new DenseVector(fun.getCoDomainDimension());
		Vector offset = new DenseVector(fun.getCoDomainDimension());
		while(actError > error && iteration < maxIterations){
			// find linearized solution
			Matrix jacobian = new DenseMatrix(fun.getDomainDimension(), fun.getCoDomainDimension());
			fun.evaluate(x, fx, offset, jacobian);
			Vector rightSide = fx.copy();
			rightSide.scale(-1);
			rightSide.add(b);
			QMR solver = new QMR(tmp);
			try {
				solver.solve(jacobian, rightSide, h);
			} catch (IterativeSolverNotConvergedException e) {
				throw new NotConvergentException(e.getMessage(), actError);
			}
			
			Vector test = new DenseVector(fun.getCoDomainDimension());
			jacobian.mult(h, test);
			
			// find stepwidth
			Vector nextX = x.copy().add(h);
			fun.evaluate(nextX, nextFx, offset);
			Vector remainder = fx.copy().add(-1, b);
			Vector nextRemainder = nextFx.copy().add(-1, b);
			while (hasNANEntry(nextFx) || Math.abs(nextRemainder.norm(Norm.Two)) > Math.abs(remainder.norm(Norm.Two))){
				h.scale(0.1);
				nextX = x.copy().add(h);
				fun.evaluate(nextX, nextFx, offset);
				nextRemainder = nextFx.copy().add(-1, b);	
			}
			// do step
			x.set(nextX);
			fx.set(nextFx);
			actError = nextRemainder.norm(Norm.Two);
			if (iterMon != null)
				iterMon.setIteration(iteration, actError);
			iteration++;
		}
	}

	
	private boolean hasNANEntry(Vector v){
		for (int i = 0; i < v.size(); i++)
			if (Double.isNaN(v.get(i)))
				return true;
		return false;
	}
	
	
	public void setIterationMonitor(IterationMonitor monitor) {
		this.iterMon = monitor;
	}

	public void setMaxIterations(Integer maxIterations) {
		this.maxIterations = maxIterations;
	}

	public Integer getMaxIterations() {
		return maxIterations;
	}

	public void setError(Double error) {
		this.error = error;
	}

	public Double getError() {
		return error;
	}

}
