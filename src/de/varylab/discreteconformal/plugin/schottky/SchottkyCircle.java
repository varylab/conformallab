package de.varylab.discreteconformal.plugin.schottky;

import static java.lang.Math.sqrt;

import java.awt.BasicStroke;
import java.awt.Color;

import de.jtem.java2d.DragListener;
import de.jtem.java2d.SceneComponent;
import de.jtem.java2d.TransformedMouseEvent;
import de.jtem.java2d.Viewer2D;
import de.jtem.java2dx.Ellipse2DDouble;
import de.jtem.mfc.field.Complex;

class SchottkyCircle implements DragListener {
	
	private boolean 
		orientation = true;
	private Complex 
		c = new Complex();
	private double 
		r = 1.0;
	private SceneComponent	
		editor = new SceneComponent();
	
	public SchottkyCircle(Complex c, double r, boolean o) {
		super();
		this.c = c;
		this.r = r;
		this.orientation = o;
		updateEditor();
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
		updateEditor();
	}
	
	public double getRadius() {
		return r;
	}
	public void setRadius(double r) {
		this.r = r;
		updateEditor();
	}
	
	public boolean getOrientation() {
		return orientation;
	}
	public void setOrientation(boolean orientation) {
		this.orientation = orientation;
		updateEditor();
	}
	
	@Override
	public String toString() {
		return "Circle: " + c + " r=" + r;
	}
	
	
	
	double lastX, lastY, rOffset;
	
	@Override
	public void dragStart(TransformedMouseEvent e) {
		lastX = e.getX();
		lastY = e.getY();
		double cdx = c.re - e.getX();
		double cdy = c.im - e.getY();
		double radius = sqrt(cdx*cdx + cdy*cdy);
		rOffset = radius - r;
	}
	
	@Override
	public void dragEnd(TransformedMouseEvent e) {
	}
	
	@Override
	public void drag(TransformedMouseEvent e) {
		switch (e.getIndexOfHitPart()) {
			case -1:
				double dx = e.getX() - lastX;
				double dy = e.getY() - lastY;
				c.setRe(c.re + dx);
				c.setIm(c.im + dy);
				break;
			case -2:
				double cdx = c.re - e.getX();
				double cdy = c.im - e.getY();
				r = sqrt(cdx*cdx + cdy*cdy) - rOffset;
				break;
		}
		lastX = e.getX();
		lastY = e.getY();
		editor.setShape(new Ellipse2DDouble(c.re-r,c.im-r,2*r,2*r));
		((Viewer2D)e.getMouseEvent().getSource()).repaint();
	}
	
	public void updateEditor() {
		editor.setShape(new Ellipse2DDouble(c.re-r,c.im-r,2*r,2*r));
		if (orientation) {
			editor.setPaint(Color.WHITE);
		} else {
			editor.setPaint(new Color(1,1,0,0.5f));
		}
		editor.setStroke(new BasicStroke(3));
		editor.setOutlineDragEnabled(true);
		editor.setAreaDragEnabled(true);
		editor.removeDragListener(this);
		editor.addDragListener(this);
	}
	
	public SceneComponent getEditor() {
		return editor;
	}
}