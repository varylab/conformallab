package de.varylab.discreteconformal.frontend.action;

import java.io.FileInputStream;
import java.io.InputStream;

import org.eclipse.jface.action.Action;
import org.eclipse.jface.dialogs.MessageDialog;
import org.eclipse.jface.resource.ImageDescriptor;
import org.eclipse.swt.SWT;
import org.eclipse.swt.widgets.FileDialog;

import de.jreality.geometry.IndexedFaceSetFactory;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.frontend.io.OFFReader;
import de.varylab.discreteconformal.image.ImageHook;

public class OpenMeshAction extends Action {
	
	private FileDialog 
		dialog = null; 
	
	public OpenMeshAction(){
		super("Open");
		setToolTipText("Open a new mesh");
		setAccelerator(SWT.CTRL | 'O');
		setImageDescriptor(ImageDescriptor.createFromImage(ImageHook.getSWTImage("open.gif")));
	}
	
	public static void openMesh(String filename) {
		try {
			FileInputStream  fin = new FileInputStream(filename);
			if (filename.toLowerCase().endsWith(".obj")) {
				openOBJMesh(fin);
			} else if (filename.toLowerCase().endsWith(".off")) {
				openOBJMesh(fin);
			} else {
				MessageDialog.openError(ConformalLab.getApplicationWindow().getShell(), "Unknown file extension", "File format unknown");
			}
		} catch (Exception e) {
			ConformalLab.handleException(e);
			return;
		} 
	}
	
	public static void openOBJMesh(InputStream in) throws Exception{
		Input input = new Input("OBJ Input", in);
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c =reader.read(input);
		ConformalLab.getGeometryController().setGeometry((IndexedFaceSet)c.getChildComponent(0).getGeometry());
	}
	
	
	public static void openOFFMesh(InputStream in) throws Exception {
		OFFReader reader = new OFFReader(in);
		IndexedFaceSetFactory ifsf = new IndexedFaceSetFactory();
		ifsf.setVertexCount(reader.getVertexCount());
		ifsf.setFaceCount(reader.getFaceCount());
		ifsf.setVertexCoordinates(reader.getVertices());
		ifsf.setFaceIndices(reader.getFaces());
		ifsf.update();
		ConformalLab.getGeometryController().setGeometry(ifsf.getIndexedFaceSet());
		ConformalLab.getUIController().encompass();
	}
	
	
	private FileDialog getFileDialog() {
		if (dialog == null) {
			dialog = new FileDialog (ConformalLab.getApplicationWindow().getShell(), SWT.OPEN);
			dialog.setFilterNames (new String [] {"OBJ Files", "OFF Files", "All Files (*.*)"});
			dialog.setFilterExtensions (new String [] {"*.obj", "*.off", "*.*"});
			dialog.setFilterPath(System.getProperty("user.dir") + "/data");
		}
		return dialog;
	}
	
	
	@Override
	public void run() {
		String filename = getFileDialog().open();
		if (filename != null) {
			openMesh(filename);
		}
		super.run();
	}
	
}
