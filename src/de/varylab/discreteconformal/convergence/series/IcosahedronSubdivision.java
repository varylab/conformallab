package de.varylab.discreteconformal.convergence.series;

import org.junit.Test;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class IcosahedronSubdivision {

	
	@Test
	public void ico14711() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Subdivision",
			"--base", "data/convergence/subdivision",
			"--name", "icosahedron1-4-7-11",
			"--pin", "data/convergence/icosahedron.obj",
			"--bpi", "1,4,7,11",
			"--max", "6"
		);
	}
	
	
}
