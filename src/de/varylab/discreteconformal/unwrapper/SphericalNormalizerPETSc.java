package de.varylab.discreteconformal.unwrapper;

import static de.jreality.math.Pn.EUCLIDEAN;

import java.lang.annotation.Annotation;
import java.util.List;
import java.util.logging.Logger;

import de.jreality.math.P3;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.GetSolutionStatusResult;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.functional.MobiusCenteringFunctional;
import de.varylab.discreteconformal.util.SparseUtility;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class SphericalNormalizerPETSc {
	
	private static Logger
		log = Logger.getLogger(SphericalNormalizerPETSc.class.getName());
	private static double
		tolerance = 1E-9;
	private static int
		maxIterations = 400;
	
	static {
		NativePathUtility.set("native");
		NativePathUtility.set("../DiscreteConformalLab/native");
		Tao.Initialize();		
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>,
		DATAGET extends Annotation, 
		DATASET extends Annotation
	> void normalize(HDS hds, AdapterSet a, Class<DATAGET> get, Class<DATASET> set) {
		normalize(hds, hds.getVertices(), a, get, set);
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>,
		DATAGET extends Annotation, 
		DATASET extends Annotation
	> void normalize(HDS hds, List<V> include, AdapterSet a, Class<DATAGET> get, Class<DATASET> set) {
		MobiusCenteringFunctional<V, E, F, DATAGET> f = new MobiusCenteringFunctional<>(include, get, a);
		TaoApplication app = f.getTaoApplication(hds);
		Vec x = new Vec(f.getDimension(hds)); x.set(0.0);
		app.setInitialSolutionVec(x);
		Mat H = SparseUtility.getHessianTemplate(f, hds);
		app.setHessianMat(H, H);
		
		log.info("staring mobius normization");
		Tao optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setGradientTolerances(tolerance, 0, 0);
		optimizer.setMaximumIterates(maxIterations);
		optimizer.solve();
		
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		UnwrapUtility.logSolutionStatus(optimizer, log);
		if (status.reason.ordinal() > 4) {
			log.warning("Mobius normalization did not converge: " + status);
		}
		double[] cm = getCenterOfMass(include, a, get);
		log.info("|CoM| before normalization: " + Pn.norm(cm, EUCLIDEAN));
		int_normalize(x, hds, a, get, set);
		cm = getCenterOfMass(include, a, get);
		log.info("|CoM| after normalization: " + Pn.norm(cm, EUCLIDEAN));
	}

	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>,
		DATAGET extends Annotation, 
		DATASET extends Annotation
	> void int_normalize(Vec c, HDS hds, AdapterSet a, Class<DATAGET> get, Class<DATASET> set){
		double[] xp = Pn.homogenize(null, c.getArray());
		double[] TInv = P3.makeTranslationMatrix(null, xp, Pn.HYPERBOLIC);
		double[] T = Rn.inverse(null, TInv);
		for (V v : hds.getVertices()) {
			double[] p = a.getD(get, v);
			Rn.matrixTimesVector(p, T, p);
			a.set(set, v, p);
		}
	}
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		DATAGET extends Annotation
	> double[] getCenterOfMass(List<V> vList, AdapterSet a, Class<DATAGET> get) {
		double[] cm = new double[4];
		double[] tmp = new double[4];
		for (V v : vList) {
			double[] p = a.getD(get, v);
			Pn.dehomogenize(tmp, p);
			Rn.add(cm, tmp, cm);
		}
		return Pn.dehomogenize(cm, cm);
	}
	
}
