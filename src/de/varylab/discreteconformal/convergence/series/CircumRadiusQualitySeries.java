package de.varylab.discreteconformal.convergence.series;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class CircumRadiusQualitySeries {
	
	public void randomCircumCircle() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Quality",
			"--base", "data/convergence/qualityCircumCircle",
			"--name", "random02-extra100-noopt",
			"--pin", "data/convergence/random02.obj",
			"--num", "10000",
			"--extra", "100",
			"--nopt", "0",
			"--QM", "MaxCircumCircle"
		);
	}
	
}
