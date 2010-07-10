package de.jtem.halfedgetools.functional.conformal.hds;

import de.varylab.discreteconformal.functional.ConformalAdapters.Alpha;
import de.varylab.discreteconformal.functional.ConformalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.ConformalAdapters.Lambda;
import de.varylab.discreteconformal.functional.ConformalAdapters.Theta;
import de.varylab.discreteconformal.functional.ConformalAdapters.Variable;


public abstract class MyConformalAdapters {


	public static class CVariable implements Variable<MyConformalVertex> {
		@Override
		public boolean isVariable(MyConformalVertex v) {
			return getVarIndex(v) >= 0;
		}
		@Override
		public int getVarIndex(MyConformalVertex v) {
			return v.getSolverIndex();
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
