package de.varylab.discreteconformal.functional;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Vertex;

public interface CAdapters {

	public static interface U <V extends Vertex<V, ?, ?>> {
		public double getU(V v);
		public void setU(V v, double u);
	}
	
	public static interface Variable <V extends Vertex<V, ?, ?>> {
		public boolean isVariable(V v);
		public int getVarIndex(V v);
	}
	
	public static interface Alpha <E extends Edge<?, E, ?>> {
		public double getAlpha(E e);
		public void setAlpha(E e, double alpha);
	}
	
	public static interface Lambda <E extends Edge<?, E, ?>> {
		public double getLambda(E e);
	}
	
	public static interface Theta <V extends Vertex<V, ?, ?>> {
		public double getTheta(V v);
	}
	
	public static interface Energy <F extends Face<?, ?, F>> {
		public double getEnergy(F f);
	}
	
	public static interface Gradient {
		public void addGradient(int i, double value);
		public void setZero();
	}

	public static interface Hessian {
		public void addHessian(int i, int j, double value);
		public void setZero();
	}
	
}
