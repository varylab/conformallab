package de.varylab.discreteconformal.frontend.shrinkpanels;

import static de.jreality.scene.data.Attribute.COLORS;
import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.LINE_SHADER;
import static de.jreality.shader.CommonAttributes.LINE_WIDTH;
import static de.jreality.shader.CommonAttributes.POINT_RADIUS;
import static de.jreality.shader.CommonAttributes.POINT_SHADER;
import static de.jreality.shader.CommonAttributes.POINT_SIZE;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.SMOOTH_SHADING;
import static de.jreality.shader.CommonAttributes.SPHERES_DRAW;
import static de.jreality.shader.CommonAttributes.TRANSPARENCY;
import static de.jreality.shader.CommonAttributes.TRANSPARENCY_ENABLED;
import static de.jreality.shader.CommonAttributes.TUBES_DRAW;
import static de.jreality.shader.CommonAttributes.TUBE_RADIUS;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static java.awt.Color.DARK_GRAY;
import static java.awt.Color.GRAY;
import static java.util.Collections.max;
import static java.util.Collections.min;
import static org.eclipse.jface.layout.GridDataFactory.fillDefaults;
import static org.eclipse.swt.SWT.CHECK;
import static org.eclipse.swt.SWT.NONE;
import static org.eclipse.swt.SWT.RADIO;
import static org.eclipse.swt.SWT.SHADOW_ETCHED_IN;
import static org.eclipse.swt.layout.GridData.BEGINNING;
import static org.eclipse.swt.layout.GridData.CENTER;

import java.awt.Color;
import java.util.Arrays;

import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.events.SelectionListener;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Scale;

import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.data.DoubleArrayArray;
import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.frontend.controller.GeometryController.GeometryChangedListener;
import de.varylab.discreteconformal.frontend.widget.ColorButton;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanel;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanelContainer;
import de.varylab.discreteconformal.frontend.widget.ColorButton.ColorChangedListener;
import de.varylab.discreteconformal.heds.CFace;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.heds.HDS;
import de.varylab.discreteconformal.heds.bsp.KdTree;
import de.varylab.discreteconformal.heds.util.MeshUtility;

public class AppearanceShrinker extends ShrinkPanel implements SelectionListener, ColorChangedListener, GeometryChangedListener{

	private final Color
		DEFAULT_VERTEX_COLOR = DARK_GRAY,
		DEFAULT_EDGE_COLOR = DARK_GRAY,
		DEFAULT_FACE_COLOR = GRAY;
	private final double
		DEFAULT_TUBE_RADIUS = 0.02;
	
	private ColorButton
		faceColorButton = null,
		edgeColorButton = null,
		vertexColorButton = null;
	private Button
		showEdges = null,
		showFaces = null,
		showVertices = null,
		vColorCurvature = null,
		vColorNone = null,
		transparentButton = null;
	private Scale
		meshWidthScale = null;
	private double[][]
	    vertexColors = null;
	private double
		meanEdgeLength = 1.0;
	
	public AppearanceShrinker(ShrinkPanelContainer parent) {
		super(parent, "Appearance");
		createLayout();
		
		showVertices.addSelectionListener(this);
		showEdges.addSelectionListener(this);
		showFaces.addSelectionListener(this);
		vertexColorButton.addColorChangedListener(this);
		edgeColorButton.addColorChangedListener(this);
		faceColorButton.addColorChangedListener(this);
		meshWidthScale.addSelectionListener(this);
		vColorNone.addSelectionListener(this);
		vColorCurvature.addSelectionListener(this);
		transparentButton.addSelectionListener(this);
		
		updateStates();
		ConformalLab.getGeometryController().addChangeListener(this);
	}

	private void createLayout() {
		setLayout(new GridLayout(2, true));
		
		Group viewGroup = new Group(this, SHADOW_ETCHED_IN);
		viewGroup.setText("View Settings");
		viewGroup.setLayout(new GridLayout(3, true));
		fillDefaults().span(2, 1).grab(true, false).applyTo(viewGroup);
	
		showVertices = new Button(viewGroup, CHECK);
		showVertices.setText("V");
		fillDefaults().grab(true, false).applyTo(showVertices);
		
		showEdges = new Button(viewGroup, CHECK);
		showEdges.setText("E");
		fillDefaults().grab(true, false).applyTo(showEdges);
		
		showFaces = new Button(viewGroup, CHECK);
		showFaces.setText("F");
		showFaces.setSelection(true);
		fillDefaults().grab(true, false).applyTo(showFaces);
		
		transparentButton = new Button(this, CHECK);
		transparentButton.setText("Transparent Faces");
		transparentButton.setSelection(false);
		fillDefaults().span(2, 1).grab(true, false).applyTo(transparentButton);
		
		Group colorGroup = new Group(this, SHADOW_ETCHED_IN);
		colorGroup.setText("Color Settings");
		colorGroup.setLayout(new GridLayout(2, true));
		fillDefaults().span(2, 1).grab(true, false).applyTo(colorGroup);
		
		Label vertexColorLabel = new Label(colorGroup, SWT.NONE);
		vertexColorLabel.setText("Vertex");
		fillDefaults().grab(true, false).align(BEGINNING, CENTER).applyTo(vertexColorLabel);
		vertexColorButton = new ColorButton(colorGroup, NONE, DEFAULT_VERTEX_COLOR);
		fillDefaults().grab(true, false).applyTo(vertexColorButton);
		
		Label edgeColorLabel = new Label(colorGroup, SWT.NONE);
		edgeColorLabel.setText("Edge");
		fillDefaults().grab(true, false).align(BEGINNING, CENTER).applyTo(edgeColorLabel);
		edgeColorButton = new ColorButton(colorGroup, NONE, DEFAULT_EDGE_COLOR);
		fillDefaults().grab(true, false).applyTo(edgeColorButton);
		
		Label meshColorLabel = new Label(colorGroup, SWT.NONE);
		meshColorLabel.setText("Face");
		fillDefaults().grab(true, false).align(BEGINNING, CENTER).applyTo(meshColorLabel);
		faceColorButton = new ColorButton(colorGroup, NONE, DEFAULT_FACE_COLOR);
		fillDefaults().grab(true, false).applyTo(faceColorButton);
		
		meshWidthScale = new Scale(this, NONE);
		meshWidthScale.setMinimum(0);
		meshWidthScale.setMaximum(100);
		meshWidthScale.setSelection((int)(DEFAULT_TUBE_RADIUS * 100));
		fillDefaults().span(2, 1).grab(true, false).applyTo(meshWidthScale);
		
		Group vColorGroup = new Group(this, SHADOW_ETCHED_IN);
		vColorGroup.setText("Face Color Settings");
		vColorGroup.setLayout(new GridLayout(2, true));
		fillDefaults().span(2, 1).grab(true, false).applyTo(vColorGroup);
		
		vColorNone = new Button(vColorGroup, RADIO);
		vColorNone.setText("Material");
		vColorNone.setSelection(true);
		fillDefaults().grab(true, false).applyTo(vColorNone);
		
		vColorCurvature = new Button(vColorGroup, RADIO);
		vColorCurvature.setText("Curvature");
		fillDefaults().grab(true, false).applyTo(vColorCurvature);
	}

	
	private void updateStates() {
		Appearance meshApp = ConformalLab.getUIController().getMeshAppearance();
	
		meshApp.setAttribute(FACE_DRAW, showFaces.getSelection());
		meshApp.setAttribute(EDGE_DRAW, showEdges.getSelection());
		meshApp.setAttribute(VERTEX_DRAW, showVertices.getSelection());
		
		meshApp.setAttribute(POLYGON_SHADER + "." + DIFFUSE_COLOR, faceColorButton.getColor());
		meshApp.setAttribute(TRANSPARENCY_ENABLED, transparentButton.getSelection());
		meshApp.setAttribute(TRANSPARENCY, transparentButton.getSelection() ? 0.8 : 0.0);
		
		meshApp.setAttribute(LINE_SHADER + "." + POLYGON_SHADER + "." + DIFFUSE_COLOR, edgeColorButton.getColor());
		meshApp.setAttribute(LINE_SHADER + "." + POLYGON_SHADER + "." + SMOOTH_SHADING, true);
		meshApp.setAttribute(LINE_SHADER + "." + TUBES_DRAW, true);
		meshApp.setAttribute(LINE_SHADER + "." + TUBE_RADIUS, meanEdgeLength * meshWidthScale.getSelection() / 100.0);
		meshApp.setAttribute(LINE_SHADER + "." + LINE_WIDTH, meanEdgeLength / 10);
		
		meshApp.setAttribute(POINT_SHADER + "." + POLYGON_SHADER + "." + DIFFUSE_COLOR, vertexColorButton.getColor());
		meshApp.setAttribute(POINT_SHADER + "." + POLYGON_SHADER + "." + SMOOTH_SHADING, true);
		meshApp.setAttribute(POINT_SHADER + "." + SPHERES_DRAW, true);
		meshApp.setAttribute(POINT_SHADER + "." + POINT_RADIUS, meanEdgeLength * meshWidthScale.getSelection() / 100.0);
		meshApp.setAttribute(POINT_SHADER + "." + POINT_SIZE, meanEdgeLength / 10);
	}
	
	
	private void updateVertexColors() {
		IndexedFaceSet ifs = ConformalLab.getGeometryController().getIndexedFaceSet();
		if (vColorCurvature.getSelection()) {
			if (vertexColors == null ) {
				HDS mesh = ConformalLab.getGeometryController().getHDS();
				double scale = meanEdgeLength * 5;
				KdTree<CVertex> kd = ConformalLab.getGeometryController().getKdTree();
				Double[] K    = new Double[mesh.numFaces()];
				for (CFace f : mesh.getFaces())
					K[f.getIndex()] = MeshUtility.absoluteCurvatureAt(f.toTriangle().getBaryCenter(), scale, kd);				
				double max = max(Arrays.asList(K));
				double min = min(Arrays.asList(K));
				double range = max - min;
				vertexColors = new double[mesh.numFaces()][];
				for (CFace f : mesh.getFaces()) {
					double k = K[f.getIndex()];
					float val = (float)((k + min) / range);
					vertexColors[f.getIndex()] = new double[] {val, 1 - val, 0};
				}
			}
			ifs.setFaceAttributes(COLORS, new DoubleArrayArray.Array(vertexColors));
		} else {
			ifs.setFaceAttributes(COLORS, null);
		}
	}
	

	public void colorChanged(ColorButton btn) {
		updateStates();
	}

	public void widgetDefaultSelected(SelectionEvent arg0) {
		
	}

	public void widgetSelected(SelectionEvent e) {
		updateStates();
		if (vColorCurvature == e.getSource() || vColorNone == e.getSource())
			updateVertexColors();
	}

	public void geometryChanged(HDS heds) {
		updateStates();
		vertexColors = null;
		vColorNone.setSelection(true);
		vColorCurvature.setSelection(false);
		meanEdgeLength = MeshUtility.meanEdgeLength(heds);
	}


}
