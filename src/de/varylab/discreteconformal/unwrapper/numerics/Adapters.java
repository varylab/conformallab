package de.varylab.discreteconformal.unwrapper.numerics;

import de.jreality.math.Pn;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Alpha;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Beta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Lambda;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Phi;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Position;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Theta;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

public abstract class Adapters {


	public static class CVariable implements Variable<CoVertex, CoEdge> {
		@Override
		public boolean isVariable(CoVertex v) {
			return getVarIndex(v) >= 0;
		}
		@Override
		public int getVarIndex(CoVertex v) {
			return v.getSolverIndex();
		}
		
		@Override
		public boolean isVariable(CoEdge e) {
			return getVarIndex(e) >= 0;
		}
		@Override
		public int getVarIndex(CoEdge e) {
			return e.getSolverIndex();
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
	
	public static class CBeta implements Beta<CoEdge> {
		@Override
		public double getBeta(CoEdge e) {
			return e.getBeta();
		}
		@Override
		public void setBeta(CoEdge e, double beta) {
			e.setBeta(beta);
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
	
	public static class CTheta implements Theta<CoVertex, CoEdge> {
		@Override
		public double getTheta(CoVertex v) {
			return v.getTheta();
		}
		@Override
		public double getTheta(CoEdge e) {
			return e.getTheta();
		}
	}
	
	public static class CPhi implements Phi<CoEdge> {
		@Override
		public double getPhi(CoEdge e) {
			return e.getPhi();
		}
	}
	
	public static class CInitialEnergy implements InitialEnergy<CoFace> {
		@Override
		public double getInitialEnergy(CoFace f) {
			return f.getInitialEnergy();
		}
	}
	
	public static class CPosition implements Position<CoVertex> {
		@Override
		public double[] getPosition(CoVertex v) {
			return Pn.dehomogenize(null, v.P);
		}
	}
	
}
