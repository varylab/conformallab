package de.varylab.discreteconformal.convergence.series;

import org.junit.Test;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class QualitySeries {

	
//	@Test
//	public void icoMultiRatioExp1() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp1",
//			"--name", "icosahedron1-3-7-11-extra200",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "1,3,7,11",
//			"--num", "200",
//			"--extra", "200",
//			"--nopt", "0",
//			"--QM", "MeanMultiRatio"
//		);
//	}
//	@Test
//	public void icoMultiRatioExp2() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp2",
//			"--name", "icosahedron1-3-7-11-extra200",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "1,3,7,11",
//			"--num", "200",
//			"--extra", "200",
//			"--nopt", "0",
//			"--exp", "2.0",
//			"--QM", "MeanMultiRatio"
//		);
//	}
//	@Test
//	public void icoMultiRatioExp10() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp10",
//			"--name", "icosahedron1-3-7-11-extra200",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "1,3,7,11",
//			"--num", "200",
//			"--extra", "200",
//			"--nopt", "0",
//			"--exp", "10.0",
//			"--QM", "MeanMultiRatio"
//		);
//	}
	
	
	@Test
	public void icoCrossRatioExp1() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Quality",
			"--base", "data/convergence/qualityCrossRatioExp1",
			"--name", "icosahedron1-3-7-11-extra200",
			"--pin", "data/convergence/icosahedron.obj",
			"--bpi", "1,3,7,11",
			"--num", "200",
			"--extra", "200",
			"--nopt", "0",
			"--QM", "MeanCrossRatio"
		);
	}
	@Test
	public void icoCrossRatioExp2() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Quality",
			"--base", "data/convergence/qualityCrossRatioExp2",
			"--name", "icosahedron1-3-7-11-extra200",
			"--pin", "data/convergence/icosahedron.obj",
			"--bpi", "1,3,7,11",
			"--num", "200",
			"--extra", "200",
			"--nopt", "0",
			"--exp", "2.0",
			"--QM", "MeanCrossRatio"
		);
	}
	@Test
	public void icoCrossRatioExp10() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Quality",
			"--base", "data/convergence/qualityCrossRatioExp10",
			"--name", "icosahedron1-3-7-11-extra200",
			"--pin", "data/convergence/icosahedron.obj",
			"--bpi", "1,3,7,11",
			"--num", "200",
			"--extra", "200",
			"--nopt", "0",
			"--exp", "10.0",
			"--QM", "MeanCrossRatio"
		);
	}
	
}
