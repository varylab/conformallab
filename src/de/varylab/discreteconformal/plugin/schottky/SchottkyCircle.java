package de.varylab.discreteconformal.plugin.schottky;

import de.jtem.mfc.field.Complex;

public class SchottkyCircle {
	
	private boolean 
		orientation = true;
	private Complex 
		center = new Complex();
	private double 
		radius = 1.0;
	
	public SchottkyCircle() {
	}
	
	public SchottkyCircle(SchottkyCircle c) {
		this.center = new Complex(c.center);
		this.orientation = c.orientation;
		this.radius = c.radius;
	}
	
	public SchottkyCircle(Complex c, double r, boolean o) {
		super();
		this.center = c;
		this.radius = r;
		this.orientation = o;
	}
	
	public boolean isInside(Complex z, double tol) {
		Complex d = z.minus(center);
		if (orientation) {
			return d.abs() < radius + tol;
		} else {
			return d.abs() > radius - tol;
		}
	}
	
	public Complex getCenter() {
		return center;
	}
	public void setCenter(Complex c) {
		this.center = c;
	}
	
	public double getRadius() {
		return radius;
	}
	public void setRadius(double r) {
		this.radius = r;
	}
	
	public boolean getOrientation() {
		return orientation;
	}
	public void setOrientation(boolean orientation) {
		this.orientation = orientation;
	}
	
	@Override
	public String toString() {
		return "Circle: " + center + " r=" + radius + " o=" + orientation;
	}
	
}