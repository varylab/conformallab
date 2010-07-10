package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import de.varylab.discreteconformal.unwrapper.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.unwrapper.UnwrapUtility.QuantizationMode;

public class CustomVertexInfo {

	public boolean	
		useCustomTheta = false;
	public double
		theta = 2 * PI;
	public BoundaryMode
		boundaryMode = BoundaryMode.Isometric;
	public QuantizationMode
		quantizationMode = QuantizationMode.AllAngles;
	
}
