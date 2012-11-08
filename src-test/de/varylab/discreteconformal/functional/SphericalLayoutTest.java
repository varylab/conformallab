package de.varylab.discreteconformal.functional;

import static java.lang.Math.PI;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import junit.framework.Assert;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import org.junit.Test;

import de.jreality.math.Rn;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition3d;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.SphericalLayout;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

public class SphericalLayoutTest {

	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CLambda
		lambda = new CLambda();
	private CAlpha
		alpha = new CAlpha();
	private CInitialEnergy
		energy = new CInitialEnergy();
	private SphericalFunctional<CoVertex, CoEdge, CoFace>
		functional = new SphericalFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	
	@Test
	public void testSphericalLayout() {
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		FunctionalTest.createOctahedron(hds, aSet);
		
		for (CoEdge e : hds.getEdges()) {
			e.setAlpha(PI / 2);
			e.setLambda(2 * log(sqrt(2)/2));
		}
		
		Vector u = new DenseVector(6);
		SphericalLayout.doLayout(hds, hds.getVertex(0), functional, u);
		
		for (CoEdge e : hds.getPositiveEdges()) {
			double[] t1 = aSet.getD(TexturePosition3d.class, e.getStartVertex());
			double[] t2 = aSet.getD(TexturePosition3d.class, e.getTargetVertex());
			double l = Rn.euclideanDistance(t1, t2);
			Assert.assertEquals(sqrt(2), l, 1E-8);
		}
		
	}
	
	
}

