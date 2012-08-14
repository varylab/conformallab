package de.varylab.discreteconformal.plugin.schottky;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;

class SchottkyGenerator {
	
	private Complex
		A = new Complex(-1),
		B = new Complex(1),
		m = new Complex(0.05);
	private SchottkyCircle
		cycle = new SchottkyCircle(new Complex(-1), 0.5, true);
	
	public SchottkyGenerator() {
	}
	
	public SchottkyGenerator(SchottkyGenerator G) {
		this(G.A, G.B, G.m, G.cycle);
	}
		
	public SchottkyGenerator(Complex A, Complex B, Complex m, SchottkyCircle cycle) {
		setFixPoints(A, B);
		setMu(m);
	}

	public Moebius getMoebius() {
		return new Moebius(A, B, m);
	}

	public void setFixPoints(Complex A, Complex B) {
		this.A.assign(A);
		this.B.assign(B);
	}
	
	public void setA(Complex a) {
		A = a;
	}
	public Complex getA() {
		return A;
	}
	public void setB(Complex b) {
		B = b;
	}
	public Complex getB() {
		return B;
	}

	public Complex getMu() {
		return m;
	}
	public void setMu(Complex m) {
		this.m.assign(m);
	}
	
	public SchottkyCircle getCycle() {
		return cycle;
	}
	public void setCycle(SchottkyCircle cycle) {
		this.cycle = cycle;
	}
	
	public Complex mapPoint(Complex z) {
		return getMoebius().applyTo(z);
	}
	
	public SchottkyCircle mapCircle(SchottkyCircle c) {
		Moebius T = getMoebius();
		Complex center = new Complex().invert(); 
		double R = T.getRadiusOfMappedCircle(c.getCenter(), c.getRadius(), center);
		double bDist = center.minus(B).abs();
		return new SchottkyCircle(center, R, bDist < R);
	}
	
	public Complex unmapPoint(Complex z) {
		return getMoebius().applyTo(z);
	}
	
	public SchottkyCircle unmapCircle(SchottkyCircle c) {
		Moebius T = getMoebius().invert();
		Complex center = new Complex(); 
		double R = T.getRadiusOfMappedCircle(c.getCenter(), c.getRadius(), center);
		double aDist = center.minus(A).abs();
		return new SchottkyCircle(center, R, aDist < R);
	}
	
}