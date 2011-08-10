package de.varylab.discreteconformal.functional;

import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import static de.varylab.discreteconformal.util.SparseUtility.makeNonZeros;

import java.io.IOException;
import java.util.Random;

import junit.framework.Assert;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomEdgeInfo;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class EuclideanConvergenceTest  {

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
	public EuclideanFunctional<CoVertex, CoEdge, CoFace>
		functional = new EuclideanFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	
	@Test
	public void testEuclideanConvergence() {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		CoHDS hds = new CoHDS(); 
		try {
			Input in = new Input("Obj File", EuclideanConvergenceTest.class.getResourceAsStream("cathead.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds converter = new ConverterJR2Heds();
			hds = new CoHDS();
			converter.ifs2heds(ifs, hds, a, null);
			hds.normalizeCoordinates();
		} catch (IOException e) {
			e.printStackTrace();
		}
		int n = UnwrapUtility.prepareInvariantDataEuclidean(hds, a);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		
//		 one edge is circular
		for (CoEdge e : hds.getPositiveEdges()) {
			CoVertex s = e.getStartVertex();
			CoVertex t = e.getTargetVertex();
			if (isBoundaryVertex(s) || isBoundaryVertex(t)) {
				continue;
			}
			e.info = new CustomEdgeInfo();
			e.info.circularHoleEdge = true;
			break;
		}
		
		CEuclideanOptimizable opt = new CEuclideanOptimizable(hds);
		// optimization
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n,makeNonZeros(hds));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CGS);
		optimizer.setError(1E-11);
		optimizer.setMaxIterations(500);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			Assert.fail(e.getLocalizedMessage());
		}
	}
	
	
}
