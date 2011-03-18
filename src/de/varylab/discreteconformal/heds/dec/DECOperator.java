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
public interface DECOperator<A extends AbstractDECObject, B extends AbstractDECObject> {

	public B apply(A obj);

}
