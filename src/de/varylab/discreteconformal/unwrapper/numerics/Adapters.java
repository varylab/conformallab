package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.halfedge.functional.Energy;
import de.jtem.halfedge.functional.conformal.CAdapters.Alpha;
import de.jtem.halfedge.functional.conformal.CAdapters.InitialEnergy;
import de.jtem.halfedge.functional.conformal.CAdapters.Lambda;
import de.jtem.halfedge.functional.conformal.CAdapters.Theta;
import de.jtem.halfedge.functional.conformal.CAdapters.Variable;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CFace;
import de.varylab.discreteconformal.heds.CVertex;

public abstract class Adapters {


	public static class CVariable implements Variable<CVertex> {
		public boolean isVariable(CVertex v) {
			return getVarIndex(v) >= 0;
		}
		public int getVarIndex(CVertex v) {
			return v.getSolverIndex();
		}
	}
	
	public static class CAlpha implements Alpha<CEdge> {
		public double getAlpha(CEdge e) {
			return e.getAlpha();
		}
		public void setAlpha(CEdge e, double alpha) {
			e.setAlpha(alpha);
		}
	}
	
	public static class CLambda implements Lambda<CEdge> {
		public double getLambda(CEdge e) {
			return e.getLambda();
		}
	}
	
	public static class CTheta implements Theta<CVertex> {
		public double getTheta(CVertex v) {
			return v.getTheta();
		}
	}
	
	public static class CInitialEnergy implements InitialEnergy<CFace> {
		public double getInitialEnergy(CFace f) {
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
