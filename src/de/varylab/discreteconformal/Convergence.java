package de.varylab.discreteconformal;

import de.jreality.util.NativePathUtility;
import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class Convergence {

	public static void main(String[] args) {
		NativePathUtility.set("native");
		try {
			ConvergenceSeries.performConvergenceSeries(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
