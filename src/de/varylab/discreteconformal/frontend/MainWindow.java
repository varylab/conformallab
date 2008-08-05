package de.varylab.discreteconformal.frontend;

import static org.eclipse.swt.SWT.BORDER;
import static org.eclipse.swt.SWT.FLAT;
import static org.eclipse.swt.SWT.NONE;

import javax.swing.JComponent;
import javax.swing.JPanel;

import org.eclipse.jface.action.MenuManager;
import org.eclipse.jface.action.Separator;
import org.eclipse.jface.action.StatusLineManager;
import org.eclipse.jface.dialogs.MessageDialog;
import org.eclipse.jface.window.ApplicationWindow;
import org.eclipse.jface.window.Window.IExceptionHandler;
import org.eclipse.swt.SWT;
import org.eclipse.swt.custom.CTabFolder;
import org.eclipse.swt.custom.CTabItem;
import org.eclipse.swt.graphics.Point;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Layout;

import swingintegration.example.EmbeddedSwingComposite;
import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.frontend.action.ExportU3DAction;
import de.varylab.discreteconformal.frontend.action.OpenMeshAction;
import de.varylab.discreteconformal.frontend.action.QuitProgramAction;
import de.varylab.discreteconformal.frontend.shrinkpanels.AppearanceShrinker;
import de.varylab.discreteconformal.frontend.shrinkpanels.ConformalUnwrapShrinker;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanelContainer;
import de.varylab.discreteconformal.image.ImageHook;

public class MainWindow extends ApplicationWindow implements IExceptionHandler{
	
	private EmbeddedSwingComposite 
		sourceViewer = null;
	
	public MainWindow(){
		super(null);
		addStatusLine();
		addMenuBar();
		addCoolBar(SWT.FLAT | SWT.WRAP);
		setExceptionHandler(this);
	}
	
	private class ControlsLayout extends Layout {

		private Control
			controls = null,
			mainView = null;
		private int 
			controlsWidth = 180;
		
		public ControlsLayout(Control controls, Control mainView, int width) {
			this.controls = controls;
			this.mainView = mainView;
			this.controlsWidth = width;
		}
		
		protected Point computeSize(Composite comp, int xh, int yh, boolean ch) {
			return comp.computeSize(SWT.DEFAULT, SWT.DEFAULT, ch);
		}

		protected void layout(Composite comp, boolean arg1) {
			Point pSize = comp.getSize();
			for (Control c : comp.getChildren()) {
				if (c == controls){
					c.setBounds(0, 0, controlsWidth, pSize.y);
				}
				if (c == mainView){
					c.setBounds(controlsWidth + 1, 0, pSize.x - controlsWidth - 2, pSize.y);
				}
			}
		}
	}
	
	
	protected Control createContents(Composite parent) {
		setStatus("Welcome");
		getShell().setText("Discrete Conformal Parametrization");
		getShell().setSize(1000,740);
		
		CTabFolder mainFolder = new CTabFolder(parent, BORDER | FLAT);


		// object content -----------------------
		Composite objectContent = new Composite(mainFolder, NONE);
		sourceViewer = new EmbeddedSwingComposite(objectContent, BORDER) {
			protected JComponent createSwingComponent() {
				JPanel panel = new JPanel();
				panel.setLayout(new java.awt.GridLayout(1, 1));
				panel.add(ConformalLab.getUIController().getViewerAppSource().getViewingComponent());
				return panel;
			}
		};
		sourceViewer.populate();
		
		ShrinkPanelContainer objectShrinkContainer = new ShrinkPanelContainer(objectContent);
		new ConformalUnwrapShrinker(objectShrinkContainer);
		new AppearanceShrinker(objectShrinkContainer);

		ControlsLayout layout1 = new ControlsLayout(objectShrinkContainer, sourceViewer, 180);
		objectContent.setLayout(layout1);

		CTabItem mainPage = new CTabItem(mainFolder, NONE);
		mainPage.setText("Object");
		mainPage.setImage(ImageHook.getImage("Standard_24i_box.png"));
		mainPage.setControl(objectContent);
		mainFolder.setSelection(mainPage);
		
		// texture content -------------------------
		Composite textureContent = new Composite(mainFolder, NONE);
		Composite textureViewer = new Composite(objectContent, BORDER); // dummy
		
		ShrinkPanelContainer textureShrinkContainer = new ShrinkPanelContainer(textureContent);
		ControlsLayout layout2 = new ControlsLayout(textureShrinkContainer, textureViewer, 180);
		textureContent.setLayout(layout2);	
		
		CTabItem texturePage = new CTabItem(mainFolder, NONE);
		texturePage.setText("Texture");
		texturePage.setImage(ImageHook.getImage("texture.png"));
		texturePage.setControl(textureContent);
		
 		return parent;
    }	
	
	@Override
	protected MenuManager createMenuManager() {
		MenuManager manager = new MenuManager();
		MenuManager fileMenu = new MenuManager("File");
		fileMenu.add(new OpenMeshAction());
		fileMenu.add(new Separator());
		fileMenu.add(new ExportU3DAction());
		fileMenu.add(new Separator());
		fileMenu.add(new QuitProgramAction(this));
		manager.add(fileMenu);
		return manager;
	}
	
//	@Override
//	protected CoolBarManager createCoolBarManager(int style) {
//		CoolBarManager cm = new CoolBarManager(style);
//		ToolBarManager fileTools = new ToolBarManager(style);
//		fileTools.add(new OpenMeshAction());
//		cm.add(fileTools);
//		return cm;
//	}
	
	@Override
	protected StatusLineManager createStatusLineManager() {
		StatusLineManager manager = new StatusLineManager();
		manager.setCancelEnabled(false);
		return manager;
	}
	

	public void handleException(Throwable e) {
		String message = e.getLocalizedMessage();
		if (message == null || message.equals(""))
			message = "No error message provided";
		MessageDialog.openError(getShell(), "Error", message);
		e.printStackTrace();
	}
	
	@Override
	protected void handleShellCloseEvent() {
		super.handleShellCloseEvent();
		ConformalLab.getUIController().getViewerAppSource().dispose();
		sourceViewer.dispose();
	}
	
}

