package de.varylab.discreteconformal.functional;

import java.util.Random;

import org.junit.Ignore;
import org.junit.Test;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomEdgeInfo;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class HyperbolicCircularHolesFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	public static final Double
		eps = 1E-5,
		error = 1E-4;
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
	private HyperbolicCircularHolesFunctional<CoVertex, CoEdge, CoFace>
		functional = new HyperbolicCircularHolesFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	
	
	@Override
	public void init() {
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		createCube(hds, aSet);
		hds.removeFace(hds.getFace(0));
		
//		one triangle of edges is circular
		for (CoFace f : hds.getFaces()) {
			if (!HalfEdgeUtils.isInteriorFace(f)) continue;
			CoEdge e1 = f.getBoundaryEdge();
			CoEdge e2 = e1.getNextEdge();
			CoEdge e3 = e2.getNextEdge();
			CustomEdgeInfo info = new CustomEdgeInfo();
			info.circularHoleEdge = true;
			e1.info = info;
			e2.info = info;
			e3.info = info;
			e1.getOppositeEdge().info = info;
			e2.getOppositeEdge().info = info;
			e3.getOppositeEdge().info = info;
			break;
		}
		
		UnwrapUtility.prepareInvariantDataHyperbolic(functional, hds, aSet);
		
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, rnd.nextDouble() - 0.5);
		}
		MyDomainValue u = new MyDomainValue(x);
		
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
		setXHessian(u);
		setEps(eps);
		setError(error);
	}
	
	@Override@Test@Ignore
	public void testHessian() throws Exception {
		super.testHessian();
	}
	
	
}
