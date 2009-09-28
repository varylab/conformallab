package de.jtem.halfedgetools.functional.conformal.hds;

import de.jtem.halfedgetools.functional.conformal.ConformalAdapters.Alpha;
import de.jtem.halfedgetools.functional.conformal.ConformalAdapters.InitialEnergy;
import de.jtem.halfedgetools.functional.conformal.ConformalAdapters.Lambda;
import de.jtem.halfedgetools.functional.conformal.ConformalAdapters.Theta;
import de.jtem.halfedgetools.functional.conformal.ConformalAdapters.Variable;


public abstract class MyConformalAdapters {


	public static class CVariable implements Variable<MyConformalVertex> {
		public boolean isVariable(MyConformalVertex v) {
			return getVarIndex(v) >= 0;
		}
		public int getVarIndex(MyConformalVertex v) {
			return v.getSolverIndex();
		}
	}
	
	public static class CAlpha implements Alpha<MyConformalEdge> {
		public double getAlpha(MyConformalEdge e) {
			return e.getAlpha();
		}
		public void setAlpha(MyConformalEdge e, double alpha) {
			e.setAlpha(alpha);
		}
	}
	
	public static class CLambda implements Lambda<MyConformalEdge> {
		public double getLambda(MyConformalEdge e) {
			return e.getLambda();
		}
		@Override
		public void setLambda(MyConformalEdge e, double lambda) {
			e.setLambda(lambda);
		}
	}
	
	public static class CTheta implements Theta<MyConformalVertex> {
		public double getTheta(MyConformalVertex v) {
			return v.getTheta();
		}
	}
	
	public static class CInitialEnergy implements InitialEnergy<MyConformalFace> {
		public double getInitialEnergy(MyConformalFace f) {
			return f.getInitialEnergy();
		}
	}
	
}
