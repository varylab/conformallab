package de.varylab.discreteconformal.frontend.widget;

import org.eclipse.swt.SWT;
import org.eclipse.swt.events.MouseEvent;
import org.eclipse.swt.events.MouseListener;
import org.eclipse.swt.events.MouseMoveListener;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.events.SelectionListener;
import org.eclipse.swt.graphics.Point;
import org.eclipse.swt.graphics.Rectangle;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Event;
import org.eclipse.swt.widgets.Layout;
import org.eclipse.swt.widgets.Listener;

public class ShrinkPanelContainer extends Composite implements SelectionListener, MouseMoveListener, Listener, MouseListener{

	private boolean 
		dragging = false;
	private int 	
		scroll = 0,
		lastMouseY = 0,
		width = 200;
	
	
	public ShrinkPanelContainer(Composite parent, int style) {
		super(parent, SWT.V_SCROLL | SWT.DOUBLE_BUFFERED | style);
		setLayout(new ShrinkLayout());
		getVerticalBar().setVisible(false);
		getVerticalBar().setMaximum(0);
		getVerticalBar().setMinimum(0);
		getVerticalBar().addSelectionListener(this);
	}

	public ShrinkPanelContainer(Composite parent){
		this(parent,SWT.NONE);
	}
	private class ShrinkLayout extends Layout{

		@Override
		protected Point computeSize(Composite comp, int xh, int yh, boolean ch) {
			Control[] children = comp.getChildren();
			int height = 0;
			for (Control c : children){
				Point cSize = c.computeSize(SWT.DEFAULT, SWT.DEFAULT, ch);
				height += cSize.y;
			}
			return new Point(width, height);
		}

		@Override
		protected void layout(Composite comp, boolean arg1) { 
			Control[] children = comp.getChildren();
			int height = 0;
			for (Control c : children)
				height += c.computeSize(SWT.DEFAULT, SWT.DEFAULT).y;
			
			Rectangle clientArea = getClientArea();
			if (clientArea.height < height){
				getVerticalBar().setMaximum(height - clientArea.height);
			} else {
				getVerticalBar().setMaximum(0);
				scroll = 0;
			}
			
			int yLoc = 0;
			for (Control c : children){
				Point cSize = c.computeSize(SWT.DEFAULT, SWT.DEFAULT);
				c.setBounds(0, yLoc - scroll, comp.getSize().x, cSize.y);
				yLoc += cSize.y;
			}
		}
	}

	public void widgetSelected(SelectionEvent arg0) {
		scroll = getVerticalBar().getSelection();
		layout();
	}


	public void widgetDefaultSelected(SelectionEvent arg0) {}


	
	public void mouseMove(MouseEvent me) {
		Control c = (Control)me.widget;
		int mouseY = c.toDisplay(me.x, me.y).y;
		if (dragging){
			scroll += lastMouseY - mouseY;
			if (scroll < 0) scroll = 0;
			if (scroll > getVerticalBar().getMaximum())
				scroll = getVerticalBar().getMaximum();
			lastMouseY = mouseY;
			layout();
		}
	}

	public void handleEvent(Event e) {
		Control c = (Control)e.widget;
		lastMouseY = c.toDisplay(e.x, e.y).y;
		dragging = true;
	}

	public void mouseDoubleClick(MouseEvent arg0) {}
	public void mouseDown(MouseEvent arg0) {}

	public void mouseUp(MouseEvent arg0) {
		dragging = false;
	}

	public int getWidth() {
		return width;
	}

	public void setWidth(int width) {
		this.width = width;
	}
	
	protected boolean isDragging() {
		return dragging;
	}

}