package de.varylab.discreteconformal.heds.dec;

/**
 * DEC - Pairing interface: Shall be implemented to define a pairing between
 * one class of DEC objects A to another one. TODO: Perhaps better to
 * implement this as an abstract class?
 * 
 * @author knoeppel
 * 
 * @param <A>
 * @param <B>
 */
@SuppressWarnings("rawtypes")
public interface DECPairing<A extends AbstractDECObject, B extends AbstractDECObject> {

	public Number pair(A obj1, B obj2);

}
