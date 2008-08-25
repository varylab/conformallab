package de.varylab.discreteconformal.frontend.shrinkpanels;

import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.TextureUtility.createTexture;
import static de.jreality.shader.TextureUtility.removeTexture;
import static java.lang.Math.PI;
import static org.eclipse.jface.layout.GridDataFactory.fillDefaults;
import static org.eclipse.swt.SWT.BORDER;
import static org.eclipse.swt.SWT.NONE;
import static org.eclipse.swt.SWT.PUSH;
import static org.eclipse.swt.layout.GridData.BEGINNING;
import static org.eclipse.swt.layout.GridData.CENTER;

import java.awt.Image;
import java.io.FileInputStream;

import javax.imageio.ImageIO;

import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.events.SelectionListener;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.FileDialog;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Spinner;

import de.jreality.math.MatrixBuilder;
import de.jreality.scene.Appearance;
import de.jreality.shader.ImageData;
import de.jreality.shader.Texture2D;
import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.frontend.controller.GeometryController.GeometryChangedListener;
import de.varylab.discreteconformal.frontend.widget.ColorButton;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanel;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanelContainer;
import de.varylab.discreteconformal.frontend.widget.ColorButton.ColorChangedListener;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.util.MeshUtility;
import de.varylab.discreteconformal.image.ImageHook;

public class TextureShrinker extends ShrinkPanel implements SelectionListener, ColorChangedListener, GeometryChangedListener{

	private double
		meanEdgeLength = 1.0;
	private Image
		texImage = ImageHook.getAWTImage("checker.jpg");
	private org.eclipse.swt.graphics.Image
		swtTexImage = ImageHook.getSWTImage("checker.jpg");
	private Texture2D 
		jrTex2D = null;
	private boolean
		jrTextureIsInvalid = false;
	private Label
		texLabel = null;
	private double
	    uTexScale = 1.0,
	    vTexScale = 1.0,
	    uTexOffset = 0.0,
	    vTexOffset = 0.0;
	private double 
		texAngle = 0.0;
	private Spinner
		uScaleSpinner = null,
		vScaleSpinner = null,
		uOffsetSpinner = null,
		vOffsetSpinner = null,
		texAngleSpinner = null;
	private Button
		loadImageBtn = null;
	private FileDialog
		dialog = null;
	
	public TextureShrinker(ShrinkPanelContainer parent) {
		super(parent, "Texture");
		createLayout();
		
		updateStates();
		ConformalLab.getGeometryController().addChangeListener(this);
		
		uScaleSpinner.addSelectionListener(this);
		vScaleSpinner.addSelectionListener(this);
		uOffsetSpinner.addSelectionListener(this);
		vOffsetSpinner.addSelectionListener(this);
		texAngleSpinner.addSelectionListener(this);
		loadImageBtn.addSelectionListener(this);
	}

	private void createLayout() {
		setLayout(new GridLayout(1, true));
		
		Group transformGroup = new Group(this, NONE);
		transformGroup.setText("Transformation");
		transformGroup.setLayout(new GridLayout(4, false));
		fillDefaults().grab(true, false).applyTo(transformGroup);
		
		Label uScaleLabel = new Label(transformGroup, NONE);
		uScaleLabel.setText("U");
		fillDefaults().grab(false, false).align(BEGINNING, CENTER).applyTo(uScaleLabel);
		uScaleSpinner = new Spinner(transformGroup, BORDER);
		uScaleSpinner.setValues((int)(uTexScale * 100), -10000, 10000, 2, 1, 10);
		fillDefaults().grab(true, false).applyTo(uScaleSpinner);
		Label vScaleLabel = new Label(transformGroup, NONE);
		vScaleLabel.setText("V");
		fillDefaults().grab(false, false).align(BEGINNING, CENTER).applyTo(vScaleLabel);
		vScaleSpinner = new Spinner(transformGroup, BORDER);
		vScaleSpinner.setValues((int)(vTexScale * 100), -10000, 10000, 2, 1, 10);
		fillDefaults().grab(true, false).applyTo(vScaleSpinner);
		
		Label uOffsetLabel = new Label(transformGroup, NONE);
		uOffsetLabel.setText("U");
		fillDefaults().grab(false, false).align(BEGINNING, CENTER).applyTo(uOffsetLabel);
		uOffsetSpinner = new Spinner(transformGroup, BORDER);
		uOffsetSpinner.setValues((int)(uTexOffset * 100), -10000, 10000, 2, 1, 10);
		fillDefaults().grab(true, false).applyTo(uOffsetSpinner);
		Label vOffsetLabel = new Label(transformGroup, NONE);
		vOffsetLabel.setText("V");
		fillDefaults().grab(false, false).align(BEGINNING, CENTER).applyTo(vOffsetLabel);
		vOffsetSpinner = new Spinner(transformGroup, BORDER);
		vOffsetSpinner.setValues((int)(vTexOffset * 100), -10000, 10000, 2, 1, 10);
		fillDefaults().grab(true, false).applyTo(vOffsetSpinner);
		
		Label texAngleLabel = new Label(transformGroup, NONE);
		texAngleLabel.setText("Angle");
		fillDefaults().span(2, 1).grab(false, false).align(BEGINNING, CENTER).applyTo(texAngleLabel);
		texAngleSpinner = new Spinner(transformGroup, BORDER);
		texAngleSpinner.setValues((int)(texAngle * 10), -3600, 3600, 1, 10, 10);
		fillDefaults().span(2, 1).grab(true, false).applyTo(texAngleSpinner);
		
		loadImageBtn = new Button(this, PUSH);
		loadImageBtn.setText("Load Image...");
		fillDefaults().grab(true, false).applyTo(loadImageBtn);
		
		texLabel = new Label(this, BORDER);
		texLabel.setImage(swtTexImage);
		fillDefaults().grab(true, false).applyTo(texLabel);
	}

	
	private void updateStates() {
		Appearance meshApp = ConformalLab.getUIController().getMeshAppearance();
		CHDS hds = ConformalLab.getGeometryController().getCHDS();
		if (hds.isTexCoordinatesValid()) {
			if (jrTex2D == null || jrTextureIsInvalid) {
				jrTex2D = createTexture(meshApp, POLYGON_SHADER, new ImageData(texImage));
				jrTextureIsInvalid = false;
			}
			MatrixBuilder texBuilder = MatrixBuilder.euclidean();
			texBuilder.scale(uTexScale / meanEdgeLength, vTexScale / meanEdgeLength, 1.0);
			texBuilder.rotate(texAngle / 180 * PI, 0.0, 0.0, 1.0);
			texBuilder.translate(uTexOffset / meanEdgeLength, vTexOffset / meanEdgeLength, 0.0);
			jrTex2D.setTextureMatrix(texBuilder.getMatrix());
		} else {
			jrTex2D = null;
			removeTexture(meshApp, POLYGON_SHADER);
		}
	}
	
	
		public void colorChanged(ColorButton btn) {
		updateStates();
	}

	public void widgetDefaultSelected(SelectionEvent arg0) {
		
	}

	public void widgetSelected(SelectionEvent e) {
		Object s = e.getSource();
		if (uScaleSpinner == s) {
			uTexScale = 100.0 / uScaleSpinner.getSelection();
		}
		if (vScaleSpinner == s) {
			vTexScale = 100.0 / vScaleSpinner.getSelection();
		}
		if (uOffsetSpinner == s) {
			uTexOffset = uOffsetSpinner.getSelection() / 100.0;
		}
		if (vOffsetSpinner == s) {
			vTexOffset = vOffsetSpinner.getSelection() / 100.0;
		}
		if (texAngleSpinner == s) {
			texAngle = texAngleSpinner.getSelection() / 10.0;
		}
		if (loadImageBtn == s) {
			String fileName = getFileDialog().open();
			if (fileName == null) return;
			try {
				FileInputStream fin = new FileInputStream(fileName);
				texImage = ImageIO.read(fin);
				if (!swtTexImage.isDisposed())
					swtTexImage.dispose();
				fin = new FileInputStream(fileName);
				swtTexImage = new org.eclipse.swt.graphics.Image(Display.getCurrent(), fin);
				texLabel.setImage(swtTexImage);
				jrTextureIsInvalid = true;
			} catch (Exception e1) {
				ConformalLab.handleException(e1);
				return;
			}
		}
		updateStates();
	}

	public void geometryChanged(final CHDS heds) {
		ConformalLab.invokeOnSWT(new Runnable() {
			public void run() {
				meanEdgeLength = MeshUtility.meanEdgeLength(heds);
				updateStates();
			}
		});
	}


	private FileDialog getFileDialog() {
		if (dialog == null) {
			dialog = new FileDialog (ConformalLab.getApplicationWindow().getShell(), SWT.OPEN);
			dialog.setFilterNames (new String [] {"Image Files(jpg, png, gif, bmp, ico, tif)", "All Files"});
			dialog.setFilterExtensions (new String [] {"*.jpg; *.jpeg; *.png; *.gif; *.bmp; *.ico; *.tif, *.tiff", "*.*"});
			dialog.setFilterPath(System.getProperty("user.dir") + "/data");
		}
		return dialog;
	}
	
	
}