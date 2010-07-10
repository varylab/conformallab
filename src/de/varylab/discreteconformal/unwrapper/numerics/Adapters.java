package de.varylab.discreteconformal.unwrapper.numerics;

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
		@Override
		public boolean isVariable(CoVertex v) {
			return getVarIndex(v) >= 0;
		}
		@Override
		public int getVarIndex(CoVertex v) {
			return v.getSolverIndex();
		}
	}
	
	public static class CAlpha implements Alpha<CoEdge> {
		@Override
		public double getAlpha(CoEdge e) {
			return e.getAlpha();
		}
		@Override
		public void setAlpha(CoEdge e, double alpha) {
			e.setAlpha(alpha);
		}
	}
	
	public static class CLambda implements Lambda<CoEdge> {
		@Override
		public double getLambda(CoEdge e) {
			return e.getLambda();
		}
		@Override
		public void setLambda(CoEdge e, double lambda) {
			e.setLambda(lambda);
		}
	}
	
	public static class CTheta implements Theta<CoVertex> {
		@Override
		public double getTheta(CoVertex v) {
			return v.getTheta();
		}
	}
	
	public static class CInitialEnergy implements InitialEnergy<CoFace> {
		@Override
		public double getInitialEnergy(CoFace f) {
			return f.getInitialEnergy();
		}
	}
	
}
