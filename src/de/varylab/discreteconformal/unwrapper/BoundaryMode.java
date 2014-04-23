package de.varylab.discreteconformal.unwrapper;

public enum BoundaryMode {
	
	Isometric,
	QuantizedAngles("Quantized Angles"),
	ConformalCurvature("Conformal Curvature"),
	Circle,
	ReadIsometricAngles("Read Isometric Angles"),
	QuantizeAnglePeriods("Quantized Angle Periods");
	
	private String name;
	
	private BoundaryMode() {
		this.name = super.toString();
	}
	private BoundaryMode(String name) {
		this.name = name;
	}
	
	@Override
	public String toString() {
		return name;
	}
	
}