package de.varylab.discreteconformal.plugin.hyperelliptic;

import java.util.List;
import java.util.Vector;

import de.jtem.java2d.CoordinateGrid;
import de.jtem.java2d.SceneComponent;
import de.jtem.java2d.Viewer2D;
import de.jtem.java2d.ViewportChangeEvent;
import de.jtem.java2d.ViewportChangeListener;

@SuppressWarnings("serial")
public abstract class Editor extends Viewer2D{
	
	protected List<EditChangeListener> listeners = new Vector<EditChangeListener>();
	protected EditChangeEvent event;
	protected SceneComponent scene;
	
	protected Editor(){
		super();
		
		final CoordinateGrid grid = new CoordinateGrid();

		addViewportChangeListener(new ViewportChangeListener() {
			public void viewportChange(ViewportChangeEvent event) {
				grid.setRectangle(event.getViewport());
				grid.fireAppearanceChange();
			}
		});

		getBackdrop().addChild(grid);

		event= new EditChangeEvent(this);
		
		initScene();
		getRoot().addChild(scene);
		repaint();
		
	}

	public void addEditChangeListener(EditChangeListener listener) {
		listeners.add(listener);
	}

	public void removeEditChangeListener(EditChangeListener listener) {
		listeners.remove(listener);
	}

	public void setEditChangeListeners(
			List<EditChangeListener> listeners) {
		this.listeners = listeners;
	}

	public List<EditChangeListener> getEditChangeListeners() {
		return listeners;
	}

	public void fireEditChangeEvent() {
		for (EditChangeListener l : getEditChangeListeners()) {
			l.editChange(event);
		}
	}
	
	protected abstract void initScene();
	public abstract void update();

}
