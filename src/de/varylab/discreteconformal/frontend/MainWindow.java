package de.varylab.discreteconformal.frontend;

import static org.eclipse.swt.SWT.BORDER;
import static org.eclipse.swt.SWT.CLOSE;
import static org.eclipse.swt.SWT.FLAT;
import static org.eclipse.swt.SWT.NONE;
import static org.eclipse.swt.SWT.TOP;

import javax.swing.JComponent;
import javax.swing.JPanel;

import org.eclipse.jface.action.CoolBarManager;
import org.eclipse.jface.action.MenuManager;
import org.eclipse.jface.action.Separator;
import org.eclipse.jface.action.StatusLineManager;
import org.eclipse.jface.action.ToolBarManager;
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
import de.varylab.discreteconformal.frontend.widget.ShrinkPanelContainer;

public class MainWindow extends ApplicationWindow implements IExceptionHandler{
	
	private int 
		folderWidth = 180;
	private EmbeddedSwingComposite 
		sourceViewer = null;
	private CTabFolder
		leftFolder = null,
		rightFolder = null;
	
	public MainWindow(){
		super(null);
		addStatusLine();
		addMenuBar();
		addCoolBar(SWT.FLAT | SWT.WRAP);
		setExceptionHandler(this);
	}
	
	private class MainLayout extends Layout {

		protected Point computeSize(Composite comp, int xh, int yh, boolean ch) {
			return comp.computeSize(SWT.DEFAULT, SWT.DEFAULT, ch);
		}

		protected void layout(Composite comp, boolean arg1) {
			Point pSize = comp.getSize();
			for (Control c : comp.getChildren()) {
				if (c == rightFolder){
					c.setBounds(pSize.x - folderWidth, 0, folderWidth, pSize.y);
				}
				if (c == leftFolder){
					c.setBounds(0, 0, folderWidth, pSize.y);
				}
				if (c == sourceViewer){
					c.setBounds(folderWidth + 1, 0, pSize.x - 2*folderWidth - 2, pSize.y);
				}
			}
		}
	}
	
	
	protected Control createContents(Composite parent) {
		setStatus("Welcome");
		getShell().setText("Discrete Conformal Parametrization");
		getShell().setSize(1000,740);
		
		Composite content = new Composite(parent, NONE);
		content.setLayout(new MainLayout());

		leftFolder = new CTabFolder(content, TOP | BORDER | CLOSE | FLAT);
		rightFolder = new CTabFolder(content, TOP | BORDER | CLOSE | FLAT);
		CTabItem mainPage = new CTabItem(leftFolder, NONE);
		CTabItem testsPage = new CTabItem(rightFolder, NONE);
		testsPage.setText("Tests");
		mainPage.setText("Controls");
		leftFolder.setSelection(mainPage);
		rightFolder.setSelection(testsPage);
		
		ShrinkPanelContainer testsShrinkContainer = new ShrinkPanelContainer(rightFolder);
//		new StableSpread(testsShrinkContainer);
//		new SpreadTests(testsShrinkContainer);
//		new QuadTestShrinker(testsShrinkContainer);
//		new CircleIntersectionTest(testsShrinkContainer);
//		new KdTreeTest(testsShrinkContainer);
//		new PrincipalTests(testsShrinkContainer);
		ShrinkPanelContainer mainShrinkContainer = new ShrinkPanelContainer(leftFolder);
		new AppearanceShrinker(mainShrinkContainer);
//		new CircleCoverTests(mainShrinkContainer);
//		new QuadMeshShrinker(mainShrinkContainer);
//		new StatusShrinker(mainShrinkContainer);
//		new PrincipalShrinker(mainShrinkContainer);
		
		
		testsPage.setControl(testsShrinkContainer);
		mainPage.setControl(mainShrinkContainer);
		
		sourceViewer = new EmbeddedSwingComposite(content, BORDER) {
			protected JComponent createSwingComponent() {
				JPanel panel = new JPanel();
				panel.setLayout(new java.awt.GridLayout(1, 1));
				panel.add(ConformalLab.getUIController().getViewerAppSource().getViewingComponent());
				return panel;
			}
		};
		sourceViewer.populate();
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
	
	@Override
	protected CoolBarManager createCoolBarManager(int style) {
		CoolBarManager cm = new CoolBarManager(style);
		ToolBarManager fileTools = new ToolBarManager(style);
		fileTools.add(new OpenMeshAction());
		cm.add(fileTools);
		return cm;
	}
	
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

