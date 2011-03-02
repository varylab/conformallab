package de.varylab.discreteconformal.plugin.schottky;

import static de.jtem.java2d.Annotation.SOUTHWEST;
import de.jtem.java2d.Annotation;
import de.jtem.java2d.DragListener;
import de.jtem.java2d.SceneComponent;
import de.jtem.java2d.TransformedMouseEvent;
import de.jtem.java2d.Viewer2D;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;

class SchottkyGenerator {
	
	private Moebius 
		s = null;
	private SchottkyCircle 
		source,
		target;
	private SceneComponent
		editor = new SceneComponent();
	private String 
		label = "1";
		
	public SchottkyGenerator(Moebius s, SchottkyCircle source, SchottkyCircle target, String label) {
		super();
		this.s = s;
		this.source = source;
		this.target = target;
		updateEditor();
	}

	public Moebius getS() {
		return s;
	}
	public void setS(Moebius s) {
		this.s = s;
	}

	public SchottkyCircle getSource() {
		return source;
	}
	public void setSource(SchottkyCircle source) {
		this.source = source;
	}

	public SchottkyCircle getTarget() {
		return target;
	}
	public void setTarget(SchottkyCircle target) {
		this.target = target;
	}
	
	public void updateEditor() {
		editor.removeAllChildren();
		SceneComponent sEditor = source.getEditor();
		SceneComponent tEditor = target.getEditor();
		sEditor.addDragListener(new DragListener() {
			@Override
			public void dragStart(TransformedMouseEvent e) {
			}
			@Override
			public void dragEnd(TransformedMouseEvent e) {
			}
			@Override
			public void drag(TransformedMouseEvent e) {
				SchottkyCircle circle = new SchottkyCircle(source.getCenter(), source.getRadius(), source.getOrientation());
				Complex centerOfMappedCircle = new Complex();
				double sr = s.getRadiusOfMappedCircle(circle.getCenter(), circle.getRadius(), centerOfMappedCircle);
				target.setCenter(centerOfMappedCircle);
				target.setRadius(sr);
				target.updateEditor();
				((Viewer2D)e.getMouseEvent().getSource()).repaint();
			}
		});
		tEditor.addDragListener(new DragListener() {
			@Override
			public void dragStart(TransformedMouseEvent e) {
			}
			@Override
			public void dragEnd(TransformedMouseEvent e) {
			}
			@Override
			public void drag(TransformedMouseEvent e) {
				SchottkyCircle circle = new SchottkyCircle(target.getCenter(), target.getRadius(), target.getOrientation());
				Complex centerOfMappedCircle = new Complex();
				Moebius t = s.invert();
				double sr = t.getRadiusOfMappedCircle(circle.getCenter(), circle.getRadius(), centerOfMappedCircle);
				source.setCenter(centerOfMappedCircle);
				source.setRadius(sr);
				source.updateEditor();
				((Viewer2D)e.getMouseEvent().getSource()).repaint();
			}
		});
		sEditor.getAnnotations().clear();
		sEditor.getAnnotations().add(new Annotation("S" + label, 0, 0, SOUTHWEST));
		tEditor.getAnnotations().clear();
		tEditor.getAnnotations().add(new Annotation("T" + label, 0, 0, SOUTHWEST));
		editor.addChild(sEditor);
		editor.addChild(tEditor);
	}
	
	
	public SceneComponent getEditor() {
		return editor;
	}
	
}