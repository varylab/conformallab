package de.varylab.discreteconformal.functional;

import static java.lang.Math.PI;

import java.util.List;
import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CBeta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

public class HyperIdealFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	public static final Double
		eps = 1E-5,
		error = 1E-4;
	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CAlpha
		alpha = new CAlpha();
	private CBeta
		beta = new CBeta();
	private HyperIdealCirclePatternFunctional<CoVertex, CoEdge, CoFace>
		functional = new HyperIdealCirclePatternFunctional<>(variable, theta, alpha, beta);
	
	
	@Override
	public void init() {
		CoHDS hds = new CoHDS(); 
		CoVertex A = hds.addNewVertex();
		CoVertex B = hds.addNewVertex();
		CoVertex C = hds.addNewVertex();
		CoVertex D = hds.addNewVertex();
		hds.addNewEdges(24);
		hds.getEdge(0).linkOppositeEdge(hds.getEdge(6));
		hds.getEdge(1).linkOppositeEdge(hds.getEdge(21));
		hds.getEdge(2).linkOppositeEdge(hds.getEdge(4));
		hds.getEdge(3).linkOppositeEdge(hds.getEdge(23));
		hds.getEdge(5).linkOppositeEdge(hds.getEdge(11));
		hds.getEdge(7).linkOppositeEdge(hds.getEdge(9));
		hds.getEdge(8).linkOppositeEdge(hds.getEdge(14));
		hds.getEdge(10).linkOppositeEdge(hds.getEdge(12));
		hds.getEdge(13).linkOppositeEdge(hds.getEdge(19));
		hds.getEdge(15).linkOppositeEdge(hds.getEdge(17));
		hds.getEdge(16).linkOppositeEdge(hds.getEdge(22));
		hds.getEdge(18).linkOppositeEdge(hds.getEdge(20));
		hds.getEdge(0).setTargetVertex(A);
		hds.getEdge(1).setTargetVertex(B);
		hds.getEdge(2).setTargetVertex(D);
		hds.getEdge(3).setTargetVertex(C);
		hds.getEdge(4).setTargetVertex(B);
		hds.getEdge(5).setTargetVertex(A);
		hds.getEdge(6).setTargetVertex(C);
		hds.getEdge(7).setTargetVertex(D);
		hds.getEdge(8).setTargetVertex(D);
		hds.getEdge(9).setTargetVertex(C);
		hds.getEdge(10).setTargetVertex(A);
		hds.getEdge(11).setTargetVertex(B);
		hds.getEdge(12).setTargetVertex(C);
		hds.getEdge(13).setTargetVertex(D);
		hds.getEdge(14).setTargetVertex(B);
		hds.getEdge(15).setTargetVertex(A);
		hds.getEdge(16).setTargetVertex(A);
		hds.getEdge(17).setTargetVertex(B);
		hds.getEdge(18).setTargetVertex(D);
		hds.getEdge(19).setTargetVertex(C);
		hds.getEdge(20).setTargetVertex(B);
		hds.getEdge(21).setTargetVertex(A);
		hds.getEdge(22).setTargetVertex(C);
		hds.getEdge(23).setTargetVertex(D);
		
		hds.getEdge(0).linkNextEdge(hds.getEdge(1));
		hds.getEdge(1).linkNextEdge(hds.getEdge(2));
		hds.getEdge(2).linkNextEdge(hds.getEdge(3));
		hds.getEdge(3).linkNextEdge(hds.getEdge(0));
		
		hds.getEdge(4).linkNextEdge(hds.getEdge(5));
		hds.getEdge(5).linkNextEdge(hds.getEdge(6));
		hds.getEdge(6).linkNextEdge(hds.getEdge(7));
		hds.getEdge(7).linkNextEdge(hds.getEdge(4));
		
		hds.getEdge(8).linkNextEdge(hds.getEdge(9));
		hds.getEdge(9).linkNextEdge(hds.getEdge(10));
		hds.getEdge(10).linkNextEdge(hds.getEdge(11));
		hds.getEdge(11).linkNextEdge(hds.getEdge(8));
		
		hds.getEdge(12).linkNextEdge(hds.getEdge(13));
		hds.getEdge(13).linkNextEdge(hds.getEdge(14));
		hds.getEdge(14).linkNextEdge(hds.getEdge(15));
		hds.getEdge(15).linkNextEdge(hds.getEdge(12));
		
		hds.getEdge(16).linkNextEdge(hds.getEdge(17));
		hds.getEdge(17).linkNextEdge(hds.getEdge(18));
		hds.getEdge(18).linkNextEdge(hds.getEdge(19));
		hds.getEdge(19).linkNextEdge(hds.getEdge(16));
		
		hds.getEdge(20).linkNextEdge(hds.getEdge(21));
		hds.getEdge(21).linkNextEdge(hds.getEdge(22));
		hds.getEdge(22).linkNextEdge(hds.getEdge(23));
		hds.getEdge(23).linkNextEdge(hds.getEdge(20));
		
		HalfEdgeUtils.fillAllHoles(hds);
		HalfEdgeUtils.isValidSurface(hds, true);
		System.out.println("genus: " + HalfEdgeUtils.getGenus(hds));
		List<CoEdge> auxEdges = Triangulator.triangulateSingleSource(hds);
		
		int index = 0;
		for (CoVertex v : hds.getVertices()) {
			v.setSolverIndex(index++);
			v.setTheta(2 * PI);
		}
		for (CoEdge e : hds.getPositiveEdges()) {
			e.setSolverIndex(index);
			e.getOppositeEdge().setSolverIndex(index++);
			if (auxEdges.contains(e)) {
				e.setTheta(PI);
				e.getOppositeEdge().setTheta(PI);
			} else {
				e.setTheta(PI/2);
				e.getOppositeEdge().setTheta(PI/2);
			}
		}
		
//		AdapterSet aSet = new ConformalAdapterSet();
//		
//		ZeroU zeroU = new ZeroU();
		
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, Math.abs(rnd.nextDouble()));
		}
		MyDomainValue u = new MyDomainValue(x);
		
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
		setEps(eps);
		setError(error);
	}
	
	
}
