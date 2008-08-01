package de.varylab.discreteconformal.frontend.shrinkpanels;

import static de.varylab.discreteconformal.ConformalLab.getGeometryController;
import static java.lang.Math.PI;
import static org.eclipse.jface.layout.GridDataFactory.fillDefaults;
import static org.eclipse.swt.SWT.BORDER;
import static org.eclipse.swt.SWT.NONE;
import static org.eclipse.swt.SWT.SHADOW_ETCHED_IN;
import static org.eclipse.swt.layout.GridData.BEGINNING;
import static org.eclipse.swt.layout.GridData.CENTER;

import java.util.Collection;
import java.util.Map;

import no.uib.cipr.matrix.DenseVector;

import org.eclipse.jface.layout.GridDataFactory;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.events.SelectionListener;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Spinner;

import de.varylab.discreteconformal.frontend.widget.ShrinkPanel;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanelContainer;
import de.varylab.discreteconformal.heds.CCones;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CLayout;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.math.optimization.NotConvergentException;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer.Solver;
import de.varylab.discreteconformal.math.optimization.stepcontrol.ArmijoStepController;

public class ConformalUnwrapShrinker extends ShrinkPanel implements SelectionListener{

	private Button
		computeEnergyBtn = null;
	private Group
		coneConfigGroup = null;
	private Spinner
		numConesSpinner = null;
	
	private int
		numCones = 0; 
	
	
	public ConformalUnwrapShrinker(ShrinkPanelContainer parent) {
		super(parent, "Conformal Unwrap");
		createLayout();
		computeEnergyBtn.addSelectionListener(this);
		numConesSpinner.addSelectionListener(this);
	}

	
	private void createLayout() {
		setLayout(new GridLayout(1, true));
		
		coneConfigGroup = new Group(this, SHADOW_ETCHED_IN);
		coneConfigGroup.setText("Cone Singularities");
		coneConfigGroup.setLayout(new GridLayout(2, true));
		GridDataFactory.fillDefaults().grab(true, false).applyTo(coneConfigGroup);
		
		Label numConesLabel = new Label(coneConfigGroup, NONE);
		numConesLabel.setText("Cones");
		fillDefaults().grab(true, false).align(BEGINNING, CENTER).applyTo(numConesLabel);
		numConesSpinner = new Spinner(coneConfigGroup, BORDER);
		numConesSpinner.setValues(numCones, 0, 100, 0, 1, 1);
		fillDefaults().grab(true, false).applyTo(numConesSpinner);
		
		
		computeEnergyBtn = new Button(this, SWT.PUSH);
		computeEnergyBtn.setText("Unwrap");
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
		if (numConesSpinner == s) {
			numCones = numConesSpinner.getSelection();
		}
	}
	
	
	private void computeEnergy(CHDS hds) {
		System.out.println("Unwrapping with " + numCones + " cones");

		hds.prepareInvariantData();
		Collection<CVertex> cones = CCones.setUpMesh(hds, numCones);
		System.out.println("Cones: " + cones);
		
		
		int n = hds.getDomainDimension();
		
		DenseVector u = new DenseVector(n);
//		double[] E = {0.0};
//		DenseVector G = new DenseVector(n);
//		Matrix H = new DenseMatrix(n, n);
		
//		hds.conformalEnergy(u, E, G, H);
//		System.err.println("Dimension: " + n);
//		System.err.println("Energy before: " + E[0]);
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

//		hds.conformalEnergy(u, E, G, H);
//		System.err.println("Solution: " + u);
//		System.err.println("Energy after: " + E[0]);
//		System.err.println("Gradient after: \n" + G);
//		System.err.println("Hessian after: \n" + H);
		
		
		System.out.println("cone angles are:");
		Map<CEdge, Double> aMap = hds.calculateAlphas(u);
		for (CVertex c : cones) {
			System.out.println(c + ": " + 180 / PI * CLayout.getAngleSum(c, aMap));
		}
		
		
		// layout
		CCones.cutMesh(hds, cones);
		CLayout.doLayout(hds, u, aMap);
		for (CVertex v : hds.getVertices()) {
			v.setPosition(v.getTextureCoord());
		}
	}
	
	
}
