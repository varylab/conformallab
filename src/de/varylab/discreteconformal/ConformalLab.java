package de.varylab.discreteconformal;

import static org.eclipse.jface.dialogs.MessageDialog.openError;

import org.eclipse.jface.window.ApplicationWindow;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Shell;

import de.varylab.discreteconformal.frontend.MainWindow;
import de.varylab.discreteconformal.frontend.action.OpenMeshAction;
import de.varylab.discreteconformal.frontend.controller.GeometryController;
import de.varylab.discreteconformal.frontend.controller.UIController;

public class ConformalLab {

	private static ApplicationWindow	
		applicationWindow = new MainWindow();
	private static GeometryController
		geometryController = new GeometryController();
	private static UIController
		uiController = new UIController();
	
	public static void main(String[] args) {
		OpenMeshAction.openMesh("data/mann.obj");
		applicationWindow.setBlockOnOpen(true);
		applicationWindow.open();
		Display.getCurrent().dispose();
		System.exit(0);
	}
	
	public static ApplicationWindow getApplicationWindow() {
		return applicationWindow;
	}
	
	public static Shell getMainShell() {
		return applicationWindow.getShell();
	}
	
	public static Display getDisplay() {
		return applicationWindow.getShell().getDisplay();
	}
	
	
	public static GeometryController getGeometryController() {
		return geometryController;
	}
	
	public static UIController getUIController() {
		return uiController;
	}
	
	
	public static void setStatus(final String status) {
		invokeOnSWT(new Runnable(){
			public void run() {
				applicationWindow.setStatus(status);				
			}
		});
	}
	
	
	public static void errorMessage(final String title, final String message) {
		invokeOnSWT(new Runnable() {
			public void run() {
				openError(getMainShell(), title, message);						
			}
		});
	}
	
	public static void invokeOnSWT(Runnable r) {
		applicationWindow.getShell().getDisplay().asyncExec(r);
	}
	
}
