package de.varylab.discreteconformal.unwrapper.numerics;

import de.varylab.discreteconformal.functional.CEuclideanFuctional.Alpha;
import de.varylab.discreteconformal.functional.CEuclideanFuctional.InitialEnergy;
import de.varylab.discreteconformal.functional.CEuclideanFuctional.Lambda;
import de.varylab.discreteconformal.functional.CEuclideanFuctional.Theta;
import de.varylab.discreteconformal.functional.CEuclideanFuctional.Variable;
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
	
}
