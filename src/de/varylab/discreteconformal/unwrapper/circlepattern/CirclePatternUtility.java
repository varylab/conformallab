package de.varylab.discreteconformal.unwrapper.circlepattern;

import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.Method;
import de.varylab.discreteconformal.unwrapper.numerics.CPEuclideanApplication;
import de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility;

public class CirclePatternUtility {

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<F, Double> calculateCirclePatternRhos(HDS hds, Map<E, Double> thetaMap, Map<F, Double> phiMap) {
		double thetaStarSum = 0.0;
		for (E e : hds.getPositiveEdges()) {
			double thStar = PI-thetaMap.get(e);
			thetaStarSum += 2*thStar;
		}
		System.out.println("theta star sum: " + thetaStarSum/(2*PI));
		for (E e : hds.getPositiveEdges()) {
			double theta = thetaMap.get(e);
			if (theta < 0) {
				System.err.println("negative theta at " + e + ": " + theta);
			}
		}
		
		int dim = hds.numFaces();
		
		Tao.Initialize();
		Vec rho = new Vec(dim);
		rho.zeroEntries();
		rho.assemble();
		
		CPEuclideanApplication<V, E, F, HDS> app = new CPEuclideanApplication<V, E, F, HDS>(hds, thetaMap, phiMap);
		app.setInitialSolutionVec(rho);
		Mat H = Mat.createSeqAIJ(dim, dim, PETSc.PETSC_DEFAULT, QuasiisothermicUtility.getPETScNonZeros(hds, app.getFunctional()));
		H.assemble();
		app.setHessianMat(H, H);
		
		Tao tao = new Tao(Method.NTR);
		tao.setApplication(app);
		tao.setMaximumIterates(200);
		tao.setTolerances(1E-15, 0, 0, 0);
		tao.setGradientTolerances(1E-15, 0, 0);
		tao.solve();
		System.out.println(tao.getSolutionStatus());
		
		Map<F, Double> rhos = new HashMap<F, Double>();
		for (F f : hds.getFaces()) {
			rhos.put(f, rho.getValue(f.getIndex()));
		}
		return rhos;
	}

}
