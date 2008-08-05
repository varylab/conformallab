package de.varylab.discreteconformal.frontend.action;

import org.eclipse.jface.action.Action;
import org.eclipse.swt.SWT;

import de.varylab.discreteconformal.frontend.MainWindow;

public class QuitProgramAction extends Action{
	
	private MainWindow
		remeshingTool = null;
	
	public QuitProgramAction(MainWindow tool){
		remeshingTool = tool;
		setToolTipText("Quit Curvature Remesher");
		setText("E&xit");
		setAccelerator(SWT.ALT | 'X');
	}
	
	@Override
	public void run() {
		remeshingTool.close();
	}
	
}