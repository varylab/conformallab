package de.varylab.discreteconformal.heds.dec;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.varylab.discreteconformal.heds.dec.AbstractDECObject.Type;

/**
 * Factory to create DEC objects on a half edge data structure.
 * 
 * @author knoeppel
 * 
 */
public class DEC <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>	   
	> {
	
	protected HalfEdgeDataStructure<V,E,F> hds;

	public DEC(HalfEdgeDataStructure<V,E,F> hds) {
		this.hds= hds;
	}
	
	public Chain<V> make0Chain() {
		return new Chain<V>(Type.VERTEX);
	}
	
	public Chain<E> make1Chain() {
		return new Chain<E>(Type.EDGE);
	}
	
	public Chain<F> make2Chain() {
		return new Chain<F>(Type.FACE);
	}
	
	public DualChain<F> makeDual0Chain() {
		return new DualChain<F>(Type.FACE);
	}
	
	public DualChain<E> makeDual1Chain() {
		return new DualChain<E>(Type.EDGE);
	}
	
	public DualChain<V> makeDual2Chain() {
		return new DualChain<V>(Type.VERTEX);
	}
	
	public Form<V> make0Form(){
		return new Form<V>(Type.VERTEX); 
	}
	
	public Form<E> make1Form(){
		return new Form<E>(Type.EDGE); 
	}
	
	public Form<F> make2Form(){
		return new Form<F>(Type.FACE); 
	}
	
	public DualForm<F> makeDual0Form(){
		return new DualForm<F>(Type.FACE);
	}
	
	public DualForm<E> makeDual1Form(){
		return new DualForm<E>(Type.EDGE);
	}
	
	public DualForm<V> makeDual2Form(){
		return new DualForm<V>(Type.VERTEX);
	}
}
