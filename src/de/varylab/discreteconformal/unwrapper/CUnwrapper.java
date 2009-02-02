package de.varylab.discreteconformal.unwrapper;

import de.varylab.discreteconformal.heds.CoHDS;

public interface CUnwrapper {

	public void unwrap(CoHDS hds) throws UnwrapException;
	
}
