package de.varylab.discreteconformal.heds.dec;

import de.jtem.halfedge.Node;

@SuppressWarnings("rawtypes")
public class Chain<N extends Node> extends AbstractDECObject<Integer> {

	protected int dimension= 0;
	
	public Chain(Type type) {
		super(type);
	}
	
}
