package de.varylab.discreteconformal.unwrapper;

import org.eclipse.core.runtime.IProgressMonitor;

import de.varylab.discreteconformal.heds.CoHDS;

public interface CUnwrapper {

	public void unwrap(CoHDS hds, IProgressMonitor mon) throws UnwrapException;
	
}
