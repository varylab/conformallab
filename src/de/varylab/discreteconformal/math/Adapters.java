package de.varylab.discreteconformal.math;

import de.varylab.discreteconformal.functional.CAdapters.Alpha;
import de.varylab.discreteconformal.functional.CAdapters.Energy;
import de.varylab.discreteconformal.functional.CAdapters.Lambda;
import de.varylab.discreteconformal.functional.CAdapters.Theta;
import de.varylab.discreteconformal.functional.CAdapters.Variable;
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
	
	public static class CEnergy implements Energy<CFace> {
		public double getEnergy(CFace f) {
			return f.getEnergy();
		}
	}
	
}
