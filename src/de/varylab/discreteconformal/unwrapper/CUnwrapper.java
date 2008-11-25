package de.varylab.discreteconformal.unwrapper;

import org.eclipse.core.runtime.IProgressMonitor;

import de.varylab.discreteconformal.heds.CHDS;

public interface CUnwrapper {

	public void unwrap(CHDS hds, IProgressMonitor mon) throws UnwrapException;
	
}
