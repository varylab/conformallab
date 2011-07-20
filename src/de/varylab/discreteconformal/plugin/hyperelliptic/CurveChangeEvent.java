package de.varylab.discreteconformal.plugin.hyperelliptic;

public class CurveChangeEvent {
	public enum EventType {
		DIVISOR_CHANGED, RATIONALDATA_CHANGED, DOUBLEVALUE_CHANGED, BRANCHPOINTS_CHANGED, CURVE_CHANGED, GENUS_CHANGED, CURVETYPE_CHANGED, NEW_CURVE_SET
	}

	public Curve curve;
	public Object source;
	public EventType type;

	public CurveChangeEvent(Curve curve, Object source,
			EventType type) {
		this.curve = curve;
		this.source = source;
		this.type = type;
	}
}
