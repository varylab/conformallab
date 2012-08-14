package de.varylab.discreteconformal.plugin.schottky;

import de.jtem.mfc.field.Complex;

class SchottkyCircle {
	
	private boolean 
		orientation = true;
	private Complex 
		c = new Complex();
	private double 
		r = 1.0;
	
	public SchottkyCircle(SchottkyCircle c) {
		this.c = new Complex(c.c);
		this.orientation = c.orientation;
		this.r = c.r;
	}
	
	public SchottkyCircle(Complex c, double r, boolean o) {
		super();
		this.c = c;
		this.r = r;
		this.orientation = o;
	}
	
	public boolean isInside(Complex z, double tol) {
		Complex d = z.minus(c);
		if (orientation) {
			return d.abs() < r + tol;
		} else {
			return d.abs() > r - tol;
		}
	}
	
	public Complex getCenter() {
		return c;
	}
	public void setCenter(Complex c) {
		this.c = c;
	}
	
	public double getRadius() {
		return r;
	}
	public void setRadius(double r) {
		this.r = r;
	}
	
	public boolean getOrientation() {
		return orientation;
	}
	public void setOrientation(boolean orientation) {
		this.orientation = orientation;
	}
	
	@Override
	public String toString() {
		return "Circle: " + c + " r=" + r + " o=" + orientation;
	}
	
}