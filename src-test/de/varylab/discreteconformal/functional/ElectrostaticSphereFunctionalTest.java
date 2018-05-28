package de.varylab.discreteconformal.functional;

import java.util.Random;

import org.junit.Ignore;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Position;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CPosition;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

public class ElectrostaticSphereFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	private Variable<CoVertex, CoEdge>
		var = new CVariable();
	private Position<CoVertex>
		pos = new CPosition();
	private ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace>	
		functional = new ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace>(var, pos);
	private Random
		rnd = new Random();
	
	@Override
	public void init() {
		CoHDS hds = new CoHDS();
		hds.addNewVertices(10);
		rnd.setSeed(0);
		
		Vector x = new DenseVector(3 * 9);
		int index = 0;
		for (CoVertex v : hds.getVertices()) {
			double[] p = {
				rnd.nextGaussian(),
				rnd.nextGaussian(),
				rnd.nextGaussian()
			};
			v.P = Pn.homogenize(null, p);
			if (v.getIndex() == 0) {
				// fixed
				v.setSolverIndex(-1);
			} else {
				// variable
				x.set(index * 3 + 0, p[0]);
				x.set(index * 3 + 1, p[1]);
				x.set(index * 3 + 2, p[2]);
				v.setSolverIndex(index++);
			}
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
