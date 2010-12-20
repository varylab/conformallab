package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import de.varylab.discreteconformal.util.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.util.UnwrapUtility.QuantizationMode;

public class CustomVertexInfo {

	public boolean	
		useCustomTheta = false;
	public double
		theta = 2 * PI;
	public BoundaryMode
		boundaryMode = BoundaryMode.Isometric;
	public QuantizationMode
		quantizationMode = QuantizationMode.AllAngles;
	

	public CustomVertexInfo() {
	}

	public CustomVertexInfo(CustomVertexInfo info) {
		this.useCustomTheta = info.useCustomTheta;
		this.theta = info.theta;
		this.boundaryMode = info.boundaryMode;
		this.quantizationMode = info.quantizationMode;
	}
	
}
