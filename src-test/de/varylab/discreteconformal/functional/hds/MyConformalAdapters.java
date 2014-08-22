package de.varylab.discreteconformal.functional.hds;

import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Lambda;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Phi;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;
import de.varylab.discreteconformal.heds.CoEdge;


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
	
	public static class CTheta implements Theta<MyConformalVertex, MyConformalEdge> {
		@Override
		public double getTheta(MyConformalVertex v) {
			return v.getTheta();
		}
		@Override
		public double getTheta(MyConformalEdge e) {
			return e.getTheta();
		}
	}
	
	public static class CPhi implements Phi<CoEdge> {
		@Override
		public double getPhi(CoEdge e) {
			return e.getPhi();
		}
	}
	
	public static class CInitialEnergy implements InitialEnergy<MyConformalFace> {
		@Override
		public double getInitialEnergy(MyConformalFace f) {
			return f.getInitialEnergy();
		}
	}
	
}
