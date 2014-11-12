package de.varylab.discreteconformal.functional;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import org.junit.Ignore;
import org.junit.Test;

import de.jreality.math.Rn;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class ElectrostaticSphereFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	private ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace>	
		functional = new ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace>();
	private Random
		rnd = new Random();
	
	@Override
	public void init() {
		CoHDS hds = new CoHDS();
		hds.addNewVertices(10);
		rnd.setSeed(0);
		Vector x = new DenseVector(30);
		for (Integer i = 0; i < x.size() / 3; i++) {
			double[] p = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			Rn.setToLength(p, p, 1.0);
			x.set(i*3 + 0, p[0]);
			x.set(i*3 + 1, p[0]);
			x.set(i*3 + 2, p[0]);
		}
		MyDomainValue u = new MyDomainValue(x);
		
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
	}
	
	
	@Override@Test@Ignore
	public void testHessian() throws Exception {
		super.testHessian();
	}
	
	
}
