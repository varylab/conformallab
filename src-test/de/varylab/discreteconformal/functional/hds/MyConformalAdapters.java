package de.varylab.discreteconformal.functional.hds;

import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Lambda;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;


public abstract class MyConformalAdapters {


	public static class CVariable implements Variable<MyConformalVertex, MyConformalEdge> {
		@Override
		public boolean isVariable(MyConformalVertex v) {
			return getVarIndex(v) >= 0;
		}
		@Override
		public int getVarIndex(MyConformalVertex v) {
			return v.getSolverIndex();
		}
		
		@Override
		public boolean isVariable(MyConformalEdge e) {
			return getVarIndex(e) >= 0;
		}
		@Override
		public int getVarIndex(MyConformalEdge e) {
			return e.getSolverIndex();
		}
	}
	
	public static class CAlpha implements Alpha<MyConformalEdge> {
		@Override
		public double getAlpha(MyConformalEdge e) {
			return e.getAlpha();
		}
		@Override
		public void setAlpha(MyConformalEdge e, double alpha) {
			e.setAlpha(alpha);
		}
	}
	
	public static class CLambda implements Lambda<MyConformalEdge> {
		@Override
		public double getLambda(MyConformalEdge e) {
			return e.getLambda();
		}
		@Override
		public void setLambda(MyConformalEdge e, double lambda) {
			e.setLambda(lambda);
		}
	}
	
	public static class CTheta implements Theta<MyConformalVertex> {
		@Override
		public double getTheta(MyConformalVertex v) {
			return v.getTheta();
		}
	}
	
	public static class CInitialEnergy implements InitialEnergy<MyConformalFace> {
		@Override
		public double getInitialEnergy(MyConformalFace f) {
			return f.getInitialEnergy();
		}
	}
	
}
