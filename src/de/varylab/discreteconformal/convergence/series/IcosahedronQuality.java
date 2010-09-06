package de.varylab.discreteconformal.convergence.series;

import org.junit.Test;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class IcosahedronQuality {

	
	@Test
	public void ico14711() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Quality",
			"--base", "data/convergence/quality",
			"--name", "icosahedron1-4-7-11-extra20",
			"--pin", "data/convergence/icosahedron.obj",
			"--bpi", "1,4,7,11",
			"--num", "200",
			"--extra", "20",
			"--nopt", "0"
		);
	}
	
	
}
