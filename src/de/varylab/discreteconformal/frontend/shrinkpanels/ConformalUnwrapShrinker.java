package de.varylab.discreteconformal.frontend.shrinkpanels;

import static de.varylab.discreteconformal.ConformalLab.getGeometryController;

import java.io.File;
import java.io.IOException;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;

import org.eclipse.jface.layout.GridDataFactory;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.events.SelectionListener;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jtem.halfedge.jReality.converter.ConverterJR2Heds;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanel;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanelContainer;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CFace;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CLayout;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.math.optimization.NotConvergentException;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer.Solver;
import de.varylab.discreteconformal.math.optimization.stepcontrol.ArmijoStepController;

public class ConformalUnwrapShrinker extends ShrinkPanel implements SelectionListener{

	private Button
		computeEnergyBtn = null;
	
	
	public ConformalUnwrapShrinker(ShrinkPanelContainer parent) {
		super(parent, "Conformal Unwrap");
		createLayout();
		computeEnergyBtn.addSelectionListener(this);
	}

	
	private void createLayout() {
		setLayout(new GridLayout(1, true));
		computeEnergyBtn = new Button(this, SWT.PUSH);
		computeEnergyBtn.setText("Compute Conformal Energy");
		GridDataFactory.fillDefaults().grab(true, false).applyTo(computeEnergyBtn);
	}


	public void widgetDefaultSelected(SelectionEvent e) {
		
	}


	public void widgetSelected(SelectionEvent e) {
		Object s = e.getSource();
		if (computeEnergyBtn == s) {
			CHDS hds = getGeometryController().getCHDS();
			computeEnergy(hds);
			getGeometryController().setGeometry(hds);
		}
	}
	
	
	private static void computeEnergy(CHDS hds) {
		hds.prepareInvariantData();
		int n = hds.getDomainDimension();
		
		DenseVector u = new DenseVector(n);
		double[] E = {0.0};
		DenseVector G = new DenseVector(n);
		Matrix H = new DenseMatrix(n, n);
		
		hds.conformalEnergy(u, E, G, H);
		System.err.println("Dimension: " + n);
		System.err.println("Energy before: " + E[0]);
//		System.err.println("Gradient before: \n" + G);
//		System.err.println("Hessian before: \n" + H);
		
		// optimization
		NewtonOptimizer optimizer = new NewtonOptimizer();
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(1E-5);
		try {
			optimizer.minimize(u, hds);
		} catch (NotConvergentException e) {
			e.printStackTrace();
		}

		
		hds.conformalEnergy(u, E, G, H);
//		System.err.println("Solution: " + u);
		System.err.println("Energy after: " + E[0]);
//		System.err.println("Gradient after: \n" + G);
//		System.err.println("Hessian after: \n" + H);
		
		
		// layout
		CLayout.doLayout(hds, u);
		System.err.println("Texture Coordinates:");
		for (CVertex v : hds.getVertices()) {
			v.setPosition(v.getTextureCoord());
			System.err.println(v.getIndex() + ": " + v.getTextureCoord());
		}
	}
	
	
	
	public static void main(String[] args) throws IOException{
		File file = new File("data/tetra.obj");
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c =reader.read(file);
		IndexedFaceSet ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
		ConverterJR2Heds<CVertex, CEdge, CFace> converter = new ConverterJR2Heds<CVertex, CEdge, CFace>(CVertex.class, CEdge.class, CFace.class);
		CHDS heds = new CHDS();
		converter.ifs2heds(ifs, heds, new PositionAdapter());
		computeEnergy(heds);
	}
	
	
}
