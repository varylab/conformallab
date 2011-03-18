package de.varylab.discreteconformal.heds.dec;

import java.util.HashMap;
import java.util.Map;

import de.jtem.halfedge.HalfEdgeDataStructure;

/**
 * Class represents an object used in DEC, in general this means a chain or
 * co-chain with the specified class W of weights. Type class contained,
 * specifying where the object acts on.
 * 
 * @author knoeppel
 * 
 * @param <W>
 */
public abstract class AbstractDECObject<W extends Number> {

	// the half edge data structure the objects operate on
	protected HalfEdgeDataStructure<?, ?, ?> hds;
	protected Map<Integer, W> weights = new HashMap<Integer, W>();

	// the type says where the object lives on
	public static enum Type {
		VERTEX, EDGE, FACE
	};

	protected Type type;
	protected int dimension;

	/**
	 * Constructs a DEC object of given type. Type class is contained in
	 * AbstractDECObject.
	 * 
	 * @param type
	 */
	public AbstractDECObject(Type type) {
		this.type = type;
	}

	/**
	 * Returns the type of the object.
	 * 
	 * @return
	 */
	public Type getType() {
		return type;
	}

	/**
	 * Returns the dimension of the space of the object.
	 * 
	 * @return
	 */
	public int getDimension() {
		return dimension;
	}

	/**
	 * Sets the dimension of the object space. TODO: i think this should be done
	 * different. Think about later!
	 * 
	 * @return
	 */
	protected void setDimension(int dimension) {
		this.dimension = dimension;
	}

	/**
	 * Get weight with specified index.
	 * 
	 * @param index
	 * @return weight
	 */
	public W getWeight(int index) {
		return this.weights.get(index);
	}

	/**
	 * Sets the weight at specified index.
	 * 
	 * @param index
	 * @param weight
	 */
	public void setWeight(int index, W weight) {
		this.weights.put(index, weight);
	}

	/**
	 * Returns an array of all weights.
	 * @return
	 */
	public Map<Integer, W> getWeights() {
		return weights;
	}
}
