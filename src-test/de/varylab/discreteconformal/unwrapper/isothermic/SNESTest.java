package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import junit.framework.Assert;

import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.MatStructure;
import de.jtem.jpetsc.SNES;
import de.jtem.jpetsc.SNES.FunctionEvaluator;
import de.jtem.jpetsc.SNES.JacobianEvaluator;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;

public class SNESTest implements FunctionEvaluator, JacobianEvaluator {
	
	static {
		NativePathUtility.set("native");
		String[] args = {
				"-help",
				"-snes_view",
				"-snes_type", "ls",
				"-snes_test_display",
				"-ksp_converged_reason",
				"-pc_factor_shift_nonzero", "1.0e-10"
		};
		Tao.Initialize("Sinus Condition Test", args, false);
	}
	
	@Test
	public void testSNES() {
		SNES snes = SNES.create();
		snes.setFromOptions();
		Vec fTpl = new Vec(2);
		Mat JTpl = new Mat(2,2);
		JTpl.assemble();
		snes.setFunction(this, fTpl);
		snes.setJacobian(this, JTpl, JTpl);
		snes.getKSP().setInitialGuessNonzero(false);
		
		Vec sol = new Vec(2);
		sol.set(0.7);
		snes.solve(null, sol);
		evaluateFunction(sol, fTpl);
		System.out.println("Solution: " + sol + " -> " + fTpl);
		Assert.assertEquals(0, fTpl.getValue(0), 1E-10);
	}

	@Override
	public MatStructure evaluateJacobian(Vec x, Mat J, Mat Jpre) {
		J.setValue(0, 0, -2*x.getValue(0), INSERT_VALUES);
		J.setValue(0, 1, 2.0, INSERT_VALUES);
		J.setValue(1, 0, 0.0, INSERT_VALUES);
		J.setValue(1, 1, 0.0, INSERT_VALUES);		
		J.assemble();
		return MatStructure.SAME_NONZERO_PATTERN;
	}

	@Override
	public void evaluateFunction(Vec x, Vec f) {
		f.setValue(0, -x.getValue(0)*x.getValue(0) + x.getValue(1) * 2 + 1, INSERT_VALUES);
		f.setValue(1, 0.0, INSERT_VALUES);
	}
	
}
