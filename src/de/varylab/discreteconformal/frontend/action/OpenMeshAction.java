package de.varylab.discreteconformal.frontend.action;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.eclipse.jface.action.Action;
import org.eclipse.jface.dialogs.MessageDialog;
import org.eclipse.jface.resource.ImageDescriptor;
import org.eclipse.swt.SWT;
import org.eclipse.swt.widgets.FileDialog;

import de.jreality.geometry.IndexedFaceSetFactory;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.frontend.io.OFFReader;
import de.varylab.discreteconformal.image.ImageHook;

public class OpenMeshAction extends Action {
	public OpenMeshAction(){
		super("Open");
		setToolTipText("Open a new mesh");
		setAccelerator(SWT.CTRL | 'O');
		setImageDescriptor(ImageDescriptor.createFromImage(ImageHook.getImage("open.gif")));
	}
	
	
	public static void openMesh(String filename) {
		try {
			File f = new File(filename);
			if (filename.toLowerCase().endsWith(".obj")) {
				ReaderOBJ reader = new ReaderOBJ();
				SceneGraphComponent c =reader.read(f);
				ConformalLab.getGeometryController().setGeometry((IndexedFaceSet)c.getChildComponent(0).getGeometry());
			} else 
			if (filename.toLowerCase().endsWith(".off")) {
				OFFReader reader = new OFFReader(filename);
				IndexedFaceSetFactory ifsf = new IndexedFaceSetFactory();
				ifsf.setVertexCount(reader.getVertexCount());
				ifsf.setFaceCount(reader.getFaceCount());
				ifsf.setVertexCoordinates(reader.getVertices());
				ifsf.setFaceIndices(reader.getFaces());
				ifsf.update();
				ConformalLab.getGeometryController().setGeometry(ifsf.getIndexedFaceSet());
			} else {
				MessageDialog.openError(ConformalLab.getApplicationWindow().getShell(), "Unknown file extension", "File format unknown");
			}
		} catch (FileNotFoundException e) {
			MessageDialog.openError(ConformalLab.getApplicationWindow().getShell(), "File not found", e.getMessage());
		} catch (IOException e) {
			MessageDialog.openError(ConformalLab.getApplicationWindow().getShell(), "IO Error", e.getMessage());
		} catch (Exception e) {
			MessageDialog.openError(ConformalLab.getApplicationWindow().getShell(), "Error", e.getMessage());			
			e.printStackTrace();
		}
	}
	
	
	
	@Override
	public void run() {
		FileDialog dialog = new FileDialog (ConformalLab.getApplicationWindow().getShell(), SWT.OPEN);
		dialog.setFilterNames (new String [] {"OBJ Files", "OFF Files", "All Files (*.*)"});
		dialog.setFilterExtensions (new String [] {"*.obj", "*.off", "*.*"});
		dialog.setFilterPath(System.getProperty("user.dir") + "/data");
		String filename = dialog.open();
		if (filename != null) {
			openMesh(filename);
		}
		super.run();
	}
	
}
