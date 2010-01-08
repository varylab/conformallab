package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.halfedgetools.functional.Energy;
import de.varylab.discreteconformal.functional.ConformalAdapters.Alpha;
import de.varylab.discreteconformal.functional.ConformalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.ConformalAdapters.Lambda;
import de.varylab.discreteconformal.functional.ConformalAdapters.Theta;
import de.varylab.discreteconformal.functional.ConformalAdapters.Variable;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

public abstract class Adapters {


	public static class CVariable implements Variable<CoVertex> {
		public boolean isVariable(CoVertex v) {
			return getVarIndex(v) >= 0;
		}
		public int getVarIndex(CoVertex v) {
			return v.getSolverIndex();
		}
	}
	
	public static class CAlpha implements Alpha<CoEdge> {
		public double getAlpha(CoEdge e) {
			return e.getAlpha();
		}
		public void setAlpha(CoEdge e, double alpha) {
			e.setAlpha(alpha);
		}
	}
	
	public static class CLambda implements Lambda<CoEdge> {
		public double getLambda(CoEdge e) {
			return e.getLambda();
		}
		@Override
		public void setLambda(CoEdge e, double lambda) {
			e.setLambda(lambda);
		}
	}
	
	public static class CTheta implements Theta<CoVertex> {
		public double getTheta(CoVertex v) {
			return v.getTheta();
		}
	}
	
	public static class CInitialEnergy implements InitialEnergy<CoFace> {
		public double getInitialEnergy(CoFace f) {
			return f.getInitialEnergy();
		}
	}
	
	
	public static class ConformalEnergy implements Energy {

		public double 
			E = 0.0;
		
		public double get() {
			return E;
		}
		
		@Override
		public void add(double E) {
			this.E += E;
		}

		@Override
		public void set(double E) {
			this.E = E;
		}

		@Override
		public void setZero() {
			E = 0.0;
		}
		
	}
	
}
