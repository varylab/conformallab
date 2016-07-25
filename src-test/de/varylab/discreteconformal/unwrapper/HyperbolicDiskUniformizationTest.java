package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryEdges;
import static de.varylab.discreteconformal.unwrapper.BoundaryMode.QuantizedAngles;
import static de.varylab.discreteconformal.unwrapper.QuantizationMode.Straight;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.toRadians;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.TypedAdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.algorithm.subdivision.LoopLinear;
import de.jtem.halfedgetools.jreality.ConverterHeds2JR;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomVertexInfo;

public class HyperbolicDiskUniformizationTest extends Assert {

	private CoHDS
		hds = null;
	private AdapterSet
		a = null;
	
	static {
		NativePathUtility.set("native");
		Tao.Initialize();
	}
	
	public static CoHDS createTriangle() {
		CoHDS hds = new CoHDS();
		CoHDS hds2 = new CoHDS();
		HalfEdgeUtils.addNGon(hds, 3);
		LoopLinear loop = new LoopLinear();
		TypedAdapterSet<double[]> a = new ConformalAdapterSet().querySet(double[].class);
		a.set(Position.class, hds.getVertex(0), new double[]{cos(2*PI/3), sin(2*PI/3)});
		a.set(Position.class, hds.getVertex(1), new double[]{cos(4*PI/3), sin(4*PI/3)});
		a.set(Position.class, hds.getVertex(2), new double[]{cos(6*PI/3), sin(6*PI/3)});
		loop.subdivide(hds, hds2, a);
		hds.clear();
		loop.subdivide(hds2, hds, a);
		return hds;
	}
	
	@Before
	public void setUp() throws Exception {
		hds = createTriangle();
		a = new ConformalAdapterSet();
	}

	@Test
	public void test() throws Exception {
		// set boundary conditions
		hds.getVertex(0).info = new CustomVertexInfo();
		hds.getVertex(1).info = new CustomVertexInfo();
		hds.getVertex(2).info = new CustomVertexInfo();
		hds.getVertex(0).info.theta = toRadians(50.0);
		hds.getVertex(1).info.theta = toRadians(50.0);
		hds.getVertex(2).info.theta = toRadians(80.0);
		EuclideanUnwrapperPETSc unwrapper = new EuclideanUnwrapperPETSc();
		unwrapper.setBoundaryQuantMode(Straight);
		unwrapper.setBoundaryMode(QuantizedAngles);
		unwrapper.setGradientTolerance(1E-12);
		unwrapper.setMaxIterations(50);
		unwrapper.unwrap(hds, 0, a);
		for (CoEdge e : boundaryEdges(hds)) {
			CoVertex v0 = e.getTargetVertex();
			CoVertex v1 = e.getStartVertex();
			CoVertex v2 = e.getNextEdge().getTargetVertex();
			double[] t0 = Pn.dehomogenize(null, v0.T);
			double[] t1 = Pn.dehomogenize(null, v1.T);
			double[] t2 = Pn.dehomogenize(null, v2.T);
			double[] u = Rn.subtract(null, t1, t0);
			double[] v = Rn.subtract(null, t2, t0);
			double alpha = Pn.angleBetween(u, v, Pn.EUCLIDEAN);
			switch (v0.getIndex()) {
				case 0: case 1: assertEquals(toRadians(50.0), alpha, 1E-6); break;
				case 2: assertEquals(toRadians(80.0), alpha, 1E-6); break;
				default: assertEquals(toRadians(180.0), alpha, 1E-6); break;
			}
		}
	}

	public static void main(String[] args) {
		ConverterHeds2JR c = new ConverterHeds2JR();
		CoHDS hds = createTriangle();
		JRViewer.display(c.heds2ifs(hds, new ConformalAdapterSet()));
	}
	
}
