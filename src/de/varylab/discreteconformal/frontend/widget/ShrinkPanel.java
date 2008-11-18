package de.varylab.discreteconformal.frontend.widget;

import org.eclipse.swt.SWT;
import org.eclipse.swt.events.DisposeEvent;
import org.eclipse.swt.events.DisposeListener;
import org.eclipse.swt.events.MouseEvent;
import org.eclipse.swt.events.MouseListener;
import org.eclipse.swt.events.MouseMoveListener;
import org.eclipse.swt.events.PaintEvent;
import org.eclipse.swt.events.PaintListener;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.Font;
import org.eclipse.swt.graphics.GC;
import org.eclipse.swt.graphics.Point;
import org.eclipse.swt.layout.FillLayout;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Layout;
import org.eclipse.swt.widgets.Shell;


public class ShrinkPanel extends Composite implements PaintListener, MouseListener, MouseMoveListener, DisposeListener{
	private static final long 
		serialVersionUID = 1L;
	private String 
		name = "";
	private int 
		inset = 3,
		name_space = 6,
		name_box_height = 15,
		name_box_middle = inset + name_box_height / 2;
	private boolean
		shrinked = false,
		floated = false;
	private static Color
		headerColor = new Color(Display.getCurrent(), 200, 200, 200),
		blackColor = new Color(Display.getCurrent(), 0,0,0),
		redColor = new Color(Display.getCurrent(), 255, 0, 0);
	private static Font
		titleFont = new Font(Display.getCurrent(), "Arial", 8, SWT.NORMAL);
	protected ShrinkPanelContainer
		parentContainer = null;
	private Shell
		floatingShell = null;
	
	private Composite
		container = null;
	
	
	public ShrinkPanel(ShrinkPanelContainer parent, String name) {
		super(parent, SWT.NONE);
		container = new Composite(parent, SWT.NORMAL);
		setParent(container);
		container.setLayout(new ShrinkPanelLayout());
		this.name = name;
		this.parentContainer = parent;
		container.addPaintListener(this);
		container.addMouseListener(this);
		container.pack(); 
		
		container.addListener(SWT.DragDetect, parent);
		container.addMouseMoveListener(parent);
		container.addMouseListener(parent);
		addListener(SWT.DragDetect, parent);
		addMouseMoveListener(parent);
		addMouseListener(parent);
	}

	
	private class ShrinkPanelLayout extends Layout{

		@Override
		protected Point computeSize(Composite comp, int xh, int yh, boolean ch) {
			if (isShrinked())
				return getShrinkedSize();
			else
				return getUnshrinkedSize();
		}

		@Override
		protected void layout(Composite comp, boolean flush) {
			Point compSize = comp.getSize();
			if (isShrinked()){
				setBounds(inset + 4, inset * 2 + name_box_height - 1, compSize.x - inset * 2 - 7, 0);
			} else {
				Point size = ShrinkPanel.this.computeSize(SWT.DEFAULT, SWT.DEFAULT);
				setBounds(inset + 4, inset * 2 + name_box_height + 2, compSize.x - inset * 2 - 7, size.y);
			}
		}
		
	}
	
	
	public void paintControl(PaintEvent e) {
		GC g = e.gc;
		Point dim = container.getSize();
		// the outline
		g.setForeground(blackColor);
		g.drawRectangle(inset, name_box_middle, dim.x - inset * 2, dim.y - name_box_middle - inset);
		g.setBackground(headerColor);
		g.fillRectangle(name_space, inset, dim.x - name_space * 2, name_box_height);
		g.setForeground(blackColor);
		g.drawRectangle(name_space, inset, dim.x - name_space * 2, name_box_height);

		//float handle
		g.setBackground(redColor);
		if (floated)
			g.fillRectangle(dim.x - name_space - name_box_height, inset + 4, name_box_height - 8, name_box_height - 8);	
		g.setForeground(blackColor);
		g.drawRectangle(dim.x - name_space - name_box_height, inset + 4, name_box_height - 8, name_box_height - 8);

		//title
		g.setFont(titleFont);		
		if (shrinked){
			g.drawString("+", name_space + inset, inset + 1, true);
		} else {
			g.drawString("-", name_space + inset, inset + 1, true);			
		}
		int headerTextWidth = 0;
		for (int i = 0; i < name.length(); i++)
			headerTextWidth += g.getAdvanceWidth(name.charAt(i));
		int offset = dim.x / 2 - headerTextWidth / 2;
		g.drawString(name, offset, inset + 1, true);
	}

	public Point getUnshrinkedSize(){
		Point innerSize = computeSize(SWT.DEFAULT, SWT.DEFAULT);
		return new Point(container.getSize().x, innerSize.y + inset * 3 + name_box_height + 6);
	}
	
	public Point getShrinkedSize(){
		return new Point(container.getSize().x, name_box_height + inset * 2);
	}
	
	
	public void widgetDisposed(DisposeEvent e) {
		setFloating(false);
	}
	
	
	public void setFloating(boolean floating){
		if (floated == floating)
			return;
		if (floating) {
			if (shrinked)
				setShrink(false);
			floated = true;
			floatingShell = new Shell(parentContainer.getShell(), SWT.DIALOG_TRIM);
			floatingShell.addDisposeListener(this);
			setParent(floatingShell);
			floatingShell.setText(name);
			floatingShell.setLayout(new FillLayout());
			floatingShell.setLocation(container.toDisplay(new Point(0,0)));
			floatingShell.setSize(container.getSize());
			floatingShell.open();
			setShrink(true);
			floatingShell.layout();
		} else {
			floated = false;
			setParent(container);
			setShrink(false);
			container.layout();
			if (!floatingShell.isDisposed())
				floatingShell.close();
		}
	}
	
	
	public void setShrink(boolean shrink){
		if (shrinked == shrink)
			return;
		Point new_size = container.getSize();
		if (shrink){
			shrinked = true;
			new_size.y = inset * 2 + name_box_height;
		} else {
			shrinked = false;
			new_size = getUnshrinkedSize();
			if (floated) setFloating(false);
		}
		container.setSize(new_size);
		parentContainer.layout();
	}

	
	public boolean isShrinked(){
		return shrinked;
	}
	
	
//	public static void main(String[] args) {
//		Display display = null;
//		if (Display.getCurrent() != null)
//			display = Display.getCurrent();
//		else
//			display = new Display();
//		Shell shell = new Shell(display);
//		shell.setText("ShrinkPanel Test");
//		shell.setSize(210, 300);
//		shell.setLayout(new FillLayout());
//		ShrinkPanelContainer spc = new ShrinkPanelContainer(shell,SWT.NONE);
//		new ShrinkPanel(spc, "Test Panel 1");
//		new ShrinkPanel(spc, "Test Panel 2");
//		new ShrinkPanel(spc, "Test Panel 3");
//		new ShrinkPanel(spc, "Test Panel 4");
//		ShrinkPanel panel2 = new ShrinkPanel(spc, "Test Panel 5");
//		GridLayout layout = new GridLayout(2, true);
//		GridData gData = new GridData(GRAB_HORIZONTAL | FILL_BOTH);
//		panel2.setLayout(layout);
//		Label l1 = new Label(panel2, SWT.NORMAL);
//		l1.setText("Test Text:");
//		Text t1 = new Text(panel2, SWT.SINGLE);
//		t1.setLayoutData(gData);
//		Label l2 = new Label(panel2, SWT.NORMAL);
//		l2.setText("Longer Text:");
//		Text t2 = new Text(panel2, SWT.SINGLE);
//		t2.setLayoutData(gData);
//		
//		shell.open();
//		while (!shell.isDisposed()){
//			if (!display.readAndDispatch ()) 
//				display.sleep ();
//		}
//		display.dispose();
//	}

	public void mouseDoubleClick(MouseEvent arg0) {}

	public void mouseDown(MouseEvent arg0) {}

	public void mouseUp(MouseEvent m) {
		if (parentContainer.isDragging())
			return;
		Point dim = container.getSize();
		if (!floated && dim.x - name_space - name_box_height <= m.x && 
		        m.x <= dim.x - name_space - 8 &&
				inset + 4 <= m.y && m.y <= inset + name_box_height - 4){
		    setFloating(true);
		} else
		if (name_space <= m.x && m.x <= dim.x - name_space &&
			inset <= m.y && m.y <= inset + name_box_height)
			setShrink(!shrinked);
	}

	
	
	public void mouseMove(MouseEvent arg0) {
		parentContainer.mouseMove(arg0);
	}

}