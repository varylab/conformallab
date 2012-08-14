package de.varylab.discreteconformal.plugin.schottky;

import java.awt.BasicStroke;
import java.awt.Color;
import java.util.Random;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import de.jtem.java2d.SceneComponent;
import de.jtem.java2dx.Ellipse2DDouble;
import de.jtem.java2dx.modelling.DraggableCircle2D;
import de.jtem.java2dx.modelling.DraggablePoint2D;
import de.jtem.java2dx.modelling.ModellingTool2D;
import de.jtem.mfc.field.Complex;

public class SchottkyGeneratorTool extends ModellingTool2D implements ChangeListener {

	private static final long 
		serialVersionUID = 1L;
	private DraggablePoint2D
		fixPointA = new DraggablePoint2D(),
		fixPointB = new DraggablePoint2D();
	private DraggableCircle2D
		sourceCircle = new DraggableCircle2D(),
		targetCircle = new DraggableCircle2D();
	private SchottkyGenerator
		generator = new SchottkyGenerator();
	private Random
		rnd = new Random();
	
	
	public SchottkyGeneratorTool(SchottkyGenerator G) {
		this();
		this.generator = G;
	}
	
	public SchottkyGeneratorTool(SchottkyGeneratorTool tool) {
		this();
		SchottkyGenerator G = tool.getModel();
		generator = new SchottkyGenerator(G);
	}
	
	public SchottkyGeneratorTool() {
		setName("Generator");
		setModelType("schottky generator");
		SceneComponent view = getViewScene();
		view.setStroke(new BasicStroke(3));
//		view.setPaint(new Color(0.3f, 0.6f, 0.0f, 0.5f));
		view.setPaint(new Color(rnd.nextFloat(), rnd.nextFloat(), 0.0f, 0.5f));
		view.addChild(sourceCircle.getViewScene());
		view.addChild(targetCircle.getViewScene());
		view.addChild(fixPointA.getViewScene());
		view.addChild(fixPointB.getViewScene());
		fixPointA.setText("A");
		fixPointB.setText("B");
		fixPointA.addChangeListener(this);
		fixPointB.addChangeListener(this);
		sourceCircle.addChangeListener(this);
		targetCircle.addChangeListener(this);
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		super.stateChanged(e);
		updateModelFromController(e.getSource() == sourceCircle || e.getSource() == fixPointB);
		updateFromModel();
	}
	
	@Override
	public SchottkyGenerator getModel() {
		return generator;
	}

	@Override
	public void storeModelState() {
	}

	@Override
	public void recallModelState() {
	}

	@Override
	public SceneComponent getControlScene() {
		return null;
	}

	@Override
	public void translate(double x, double y) {
		generator.getA().re += x;
		generator.getA().im += y;
		generator.getB().re += x;
		generator.getB().im += y;
		generator.getCycle().getCenter().re += x;
		generator.getCycle().getCenter().im += y;
	}

	@Override
	protected void updateFromModel() {
		Complex A = generator.getA();
		Complex B = generator.getB();
		fixPointA.getModel().x = A.re;
		fixPointA.getModel().y = A.im;
		fixPointA.updateFromModel();
		fixPointB.getModel().x = B.re;
		fixPointB.getModel().y = B.im;
		fixPointB.updateFromModel();
		SchottkyCircle sCircle = generator.getCycle();
		SchottkyCircle tCircle = generator.mapCircle(sCircle);
		Complex sCenter = sCircle.getCenter();
		Complex tCenter = tCircle.getCenter();
		sourceCircle.getModel().x = sCenter.re - sCircle.getRadius();
		sourceCircle.getModel().y = sCenter.im - sCircle.getRadius();
		sourceCircle.getModel().width = sCircle.getRadius() * 2;
		sourceCircle.getModel().height = sCircle.getRadius() * 2;
		sourceCircle.updateFromModel();
		targetCircle.getModel().x = tCenter.re - tCircle.getRadius();
		targetCircle.getModel().y = tCenter.im - tCircle.getRadius();
		targetCircle.getModel().width = tCircle.getRadius() * 2;
		targetCircle.getModel().height = tCircle.getRadius() * 2;
		targetCircle.updateFromModel();
		fixPointA.getViewScene().fireAppearanceChange();
		fixPointB.getViewScene().fireAppearanceChange();
		sourceCircle.getViewScene().fireAppearanceChange();
		targetCircle.getViewScene().fireAppearanceChange();
		getViewScene().fireAppearanceChange();
	}
	
	protected void updateModelFromController(boolean source) {
		generator.getA().re = fixPointA.getModel().x;
		generator.getA().im = fixPointA.getModel().y;
		generator.getB().re = fixPointB.getModel().x;
		generator.getB().im = fixPointB.getModel().y;
		if (source) {
			SchottkyCircle sCircle = generator.getCycle();
			Ellipse2DDouble sEllipse = sourceCircle.getModel();
			double sRadius = sEllipse.width / 2;
			Complex sCenter = new Complex(sEllipse.x + sRadius, sEllipse.y + sRadius);
			sCircle.setCenter(sCenter);
			sCircle.setRadius(sRadius);
			Complex A = generator.getA();
			sCircle.setOrientation(sEllipse.contains(A.re, A.im));
		} else {
			Ellipse2DDouble tEllipse = targetCircle.getModel();
			double tRadius = tEllipse.width / 2;
			Complex tCenter = new Complex(tEllipse.x + tRadius, tEllipse.y + tRadius);
			SchottkyCircle bCircle = new SchottkyCircle(tCenter, tRadius, true); 
			SchottkyCircle aCircle = generator.unmapCircle(bCircle);
			generator.setCycle(aCircle);
		}
		updateFromModel();
	}

	@Override
	public ModellingTool2D clone() {
		return new SchottkyGeneratorTool(this);
	}
	

}
