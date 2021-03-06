package de.varylab.discreteconformal.convergence.series;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class QualitySeries {

	
//	@Test
//	public void icoMultiRatioExp1() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp1",
//			"--name", "icosahedron0-1-10-11-extra50-randopt",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--num", "100000",
//			"--extra", "50",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "1.0",
//			"--QM", "MeanMultiRatio"
//		);
//	}
//	@Test
//	public void icoMultiRatioExp2() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp2",
//			"--name", "icosahedron0-1-10-11-extra50-randopt",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--num", "100000",
//			"--extra", "50",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "2.0",
//			"--QM", "MeanMultiRatio"
//		);
//	}
//	@Test
//	public void icoMultiRatioExp10() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp10",
//			"--name", "icosahedron0-1-10-11-extra50-randopt",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--num", "100000",
//			"--extra", "50",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "1.0",
//			"--QM", "MeanMultiRatio"
//		);
//	}
//	@Test
//	public void icoMultiRatioExp05() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp05",
//			"--name", "icosahedron0-1-10-11-extra50-randopt",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--num", "100000",
//			"--extra", "50",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "0.5",
//			"--QM", "MeanMultiRatio"
//		);
//	}
//	
	
//	@Test
//	public void icoCrossRatioExp1() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityCrossRatioExp1",
//			"--name", "icosahedron0-1-10-11-extra50-randopt",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--num", "100000",
//			"--extra", "50",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "1.0",
//			"--QM", "MeanCrossRatio"
//		);
//	}
//	@Test
//	public void icoCrossRatioExp2() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityCrossRatioExp2",
//			"--name", "icosahedron0-1-10-11-extra50-randopt",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--num", "100000",
//			"--extra", "50",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "2.0",
//			"--QM", "MeanCrossRatio"
//		);
//	}
//	@Test
//	public void icoCrossRatioExp10() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityCrossRatioExp10",
//			"--name", "icosahedron0-1-10-11-extra50-randopt",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--num", "100000",
//			"--extra", "50",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "10.0",
//			"--QM", "MeanCrossRatio"
//		);
//	}
//	
//	public void icoCrossRatioExp05() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityCrossRatioExp05",
//			"--name", "icosahedron0-1-10-11-extra50-randopt",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--num", "100000",
//			"--extra", "50",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "0.5",
//			"--QM", "MeanCrossRatio"
//		);
//	}
	
/** random configuration **/
	
	public void randomMultiRatioExp1RandOpt() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Quality",
			"--base", "data/convergence/qualityMultiRatioExp1",
			"--name", "random02-extra100-randopt",
			"--pin", "data/convergence/random02.obj",
			"--num", "10000",
			"--extra", "100",
			"--nopt", "50",
			"--noptrand",
			"--exp", "1.0",
			"--QM", "MeanMultiRatio"
		);
	}
	public void randomMultiRatioExp1NoOpt() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Quality",
			"--base", "data/convergence/qualityMultiRatioExp1",
			"--name", "random02-extra100-noopt",
			"--pin", "data/convergence/random02.obj",
			"--num", "10000",
			"--extra", "100",
			"--nopt", "0",
			"--exp", "1.0",
			"--QM", "MeanMultiRatio"
		);
	}	
	
//	@Test
//	public void randomMultiRatioExp2() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp2",
//			"--name", "random02-extra50-randopt",
//			"--pin", "data/convergence/random02.obj",
//			"--num", "10000",
//			"--extra", "40",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "2.0",
//			"--QM", "MeanMultiRatio"
//		);
//	}
//	@Test
//	public void randomMultiRatioExp10() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp10",
//			"--name", "random02-extra50-randopt",
//			"--pin", "data/convergence/random02.obj",
//			"--num", "10000",
//			"--extra", "40",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "1.0",
//			"--QM", "MeanMultiRatio"
//		);
//	}
//	@Test
//	public void randomMultiRatioExp05() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Quality",
//			"--base", "data/convergence/qualityMultiRatioExp05",
//			"--name", "random02-extra50-randopt",
//			"--pin", "data/convergence/random02.obj",
//			"--num", "10000",
//			"--extra", "40",
//			"--nopt", "20",
//			"--noptrand",
//			"--exp", "0.5",
//			"--QM", "MeanMultiRatio"
//		);
//	}
	
}
