package de.varylab.discreteconformal.frontend.action;

import org.eclipse.jface.action.Action;
import org.eclipse.jface.resource.ImageDescriptor;
import org.eclipse.swt.SWT;

import de.varylab.discreteconformal.frontend.MainWindow;
import de.varylab.discreteconformal.image.ImageHook;

public class QuitProgramAction extends Action{
	
	private MainWindow
		remeshingTool = null;
	
	public QuitProgramAction(MainWindow tool){
		remeshingTool = tool;
		setToolTipText("Quit Curvature Remesher");
		setText("E&xit");
		setAccelerator(SWT.ALT | 'X');
		setImageDescriptor(ImageDescriptor.createFromImage(ImageHook.getImage("close.png")));
	}
	
	@Override
	public void run() {
		remeshingTool.close();
	}
	
}