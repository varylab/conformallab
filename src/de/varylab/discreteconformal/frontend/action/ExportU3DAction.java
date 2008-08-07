package de.varylab.discreteconformal.frontend.action;

import static de.varylab.discreteconformal.ConformalLab.getApplicationWindow;
import static org.eclipse.jface.resource.ImageDescriptor.createFromImage;

import java.io.File;
import java.io.FileOutputStream;

import org.eclipse.jface.action.Action;
import org.eclipse.jface.dialogs.MessageDialog;
import org.eclipse.swt.SWT;
import org.eclipse.swt.widgets.FileDialog;

import de.jreality.writer.u3d.WriterU3D;
import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.image.ImageHook;

public class ExportU3DAction extends Action {
	
	private FileDialog 
		dialog = null;
	
	public ExportU3DAction(){
		super("Export U3D...");
		setToolTipText("Export current scene as U3D file");
		setAccelerator(SWT.CTRL | 'E');
		setImageDescriptor(createFromImage(ImageHook.getSWTImage("save.png")));
	}
	
	
	protected FileDialog getDialog() {
		if (dialog == null) {
			dialog = new FileDialog (getApplicationWindow().getShell(), SWT.SAVE);
			dialog.setFilterNames (new String [] {"U3D Files", "All Files (*.*)"});
			dialog.setFilterExtensions (new String [] {"*.u3d", "*.*"});
			dialog.setFilterPath(System.getProperty("user.dir") + "/data");
		}
		return dialog;
	}
	
	
	@Override
	public void run() {
		String filename = getDialog().open();
		if (filename == null) return;
		File file = new File(filename);
		if (!file.getName().toLowerCase().endsWith(".u3d"))
			file = new File(file.getAbsolutePath() + ".u3d");
		if (file.exists()) {
			boolean overwrite = MessageDialog.openConfirm(ConformalLab.getMainShell(), 
					"File exists", "Do you want to overwrite the file: " + file + "?");
			if (!overwrite) return;
		}
		try {
			FileOutputStream fos = new FileOutputStream(file);
			WriterU3D writer = new WriterU3D();
			writer.writeScene(ConformalLab.getUIController().getViewerAppSource().getJrScene(), fos);
			fos.close();
		} catch (Exception e1) {
			e1.printStackTrace();
			MessageDialog.openError(ConformalLab.getMainShell(), e1.getClass().getSimpleName(), e1.getMessage());
		}
	}
	
}
