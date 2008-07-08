package de.varylab.discreteconformal.frontend.widget;

import java.util.LinkedList;

import org.eclipse.swt.SWT;
import org.eclipse.swt.events.PaintEvent;
import org.eclipse.swt.events.PaintListener;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.events.SelectionListener;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.GC;
import org.eclipse.swt.graphics.Point;
import org.eclipse.swt.graphics.RGB;
import org.eclipse.swt.layout.FillLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.ColorDialog;
import org.eclipse.swt.widgets.Composite;

import de.varylab.discreteconformal.ConformalLab;

public class ColorButton extends Composite implements SelectionListener, PaintListener {

	private Button
		colorButton = null;
	private RGB
		color = new RGB(1.0f, 1.0f, 1.0f);
	private LinkedList<ColorChangedListener>
		listeners = new LinkedList<ColorChangedListener>();
	private ColorDialog 
		colorDialog = new ColorDialog(ConformalLab.getApplicationWindow().getShell());
	
	public static interface ColorChangedListener {
		
		public void colorChanged(ColorButton btn);
		
	}
	
	public ColorButton(Composite parent, int style, java.awt.Color color) {
		super(parent, style);
		setColor(color);
		setLayout(new FillLayout());
		colorButton = new Button(this, SWT.PUSH);
		colorButton.addSelectionListener(this);
		colorButton.addPaintListener(this);
	}

	
	public void setText(String text) {
		colorDialog.setText(text);
	}
	
	public boolean addColorChangedListener(ColorChangedListener l){
		return listeners.add(l);
	}
	
	public boolean removeColorChangedListener(ColorChangedListener l){
		return listeners.remove(l);
	}
	
	private void fireColorChangedEvent() {
		for (ColorChangedListener l : listeners)
			l.colorChanged(this);
	}
	
	
	public void setColor(java.awt.Color color) {
		this.color = new RGB(color.getRed(), color.getGreen(), color.getBlue());
	}
	
	public java.awt.Color getColor() {
		java.awt.Color c = new java.awt.Color(color.red, color.green, color.blue);
		return c;
	}


	public void widgetDefaultSelected(SelectionEvent arg0) {
		
	}


	public void widgetSelected(SelectionEvent arg0) {
		colorDialog.setRGB(color);
		RGB color = colorDialog.open();
		if (color == null)
			return;
		this.color = colorDialog.getRGB();
		fireColorChangedEvent();
	}


	public void paintControl(PaintEvent e) {
		Color c = new Color(colorButton.getDisplay(), color.red, color.green, color.blue);
		GC g = e.gc;
		g.setBackground(c);
		Point dim = colorButton.getSize();
		g.fillRectangle(6, 5, dim.x-13, dim.y-11);
		c.dispose();
	}
	
}
