package de.varylab.discreteconformal.heds.dec;

/**
 * DEC - Operator interface: Shall be implemented to define an operator mapping
 * from a class of DEC object A to another one. TODO: Perhaps better to
 * implement this as an abstract class?
 * 
 * @author knoeppel
 * 
 * @param <A>
 * @param <B>
 */
@SuppressWarnings("rawtypes")
public abstract class AbstractDECOperator<A extends AbstractDECObject, B extends AbstractDECObject>
extends AbstractDECObject<Number>
{

	public AbstractDECOperator(Type type) {
		super(type);
	}

	public abstract B apply(A obj);
	
}
