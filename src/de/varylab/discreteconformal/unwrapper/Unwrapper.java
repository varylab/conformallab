package de.varylab.discreteconformal.unwrapper;

import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.heds.CoHDS;

public interface Unwrapper {

	public Vector unwrap(CoHDS surface) throws Exception;
	
}
