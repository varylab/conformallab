package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryVertices;
import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.outgoingEdges;
import static de.varylab.discreteconformal.util.UnwrapUtility.prepareInvariantDataEuclidean;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import no.uib.cipr.matrix.DenseVector;
import de.jreality.geometry.IndexedFaceSetUtility;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.data.Attribute;
import de.jreality.scene.data.DataList;
import de.jreality.scene.data.DoubleArrayArray;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomVertexInfo;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.heds.adapter.MappedEdgeLengthAdapter;
import de.varylab.discreteconformal.math.ComplexUtility;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.mtjoptimization.NotConvergentException;


public class CircleDomainUnwrapper implements Unwrapper{

	private Logger
		log = Logger.getLogger(CircleDomainUnwrapper.class.getName());
	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	private CoVertex
		unwrapRoot = null,
		layoutRoot = null;
	

	static {
		Tao.Initialize();		
	}
	
	
	public static void unwrap(IndexedFaceSet ifs) throws Exception {
		try {
		ifs.setVertexAttributes(Attribute.TEXTURE_COORDINATES, null);
		IndexedFaceSetUtility.makeConsistentOrientation(ifs);
		CoHDS hds = new CoHDS();
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		ConverterJR2Heds cTo = new ConverterJR2Heds();
		cTo.ifs2heds(ifs, hds, a);
		Set<CoVertex> oldVertices = new HashSet<CoVertex>(hds.getVertices());

		// store initial indices
		Map<CoVertex, Integer> indexMap = new HashMap<CoVertex, Integer>();
		for (CoVertex v : hds.getVertices()) {
			indexMap.put(v, v.getIndex());
		}
		
		// unwrap
		CircleDomainUnwrapper unwrapper = new CircleDomainUnwrapper();
		unwrapper.unwrap(hds, 0, a);
		
		oldVertices.removeAll(hds.getVertices());
		int oldV0Index = oldVertices.iterator().next().getIndex(); 
		
		double[][] texArr = new double[ifs.getNumPoints()][];
		for (CoVertex v : hds.getVertices()) {
			Integer oldI = indexMap.get(v);
			if (oldI == null) {
				oldI = oldV0Index;
			}
			texArr[oldI] = v.T;
		}
		
		DataList newTexCoords = new DoubleArrayArray.Array(texArr);
		ifs.setVertexAttributes(Attribute.TEXTURE_COORDINATES, newTexCoords);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		CoVertex v0 = unwrapRoot;
		// find boundary vertex with high valence
		if (v0 == null || !HalfEdgeUtils.isBoundaryVertex(v0)) {
			int maxValence = 0;
			for (CoVertex v : boundaryVertices(surface)) {
				int val = incomingEdges(v).size();
				if (val > maxValence) {
					v0 = v;
					maxValence = val;
				}
			}
		}
		
		// find reference vertices
		CoVertex v1 = null;
		CoVertex v2 = null;
		CoVertex v3 = null;
		double meanLength = 0.0;
		int numIncident = 0;
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			meanLength += aSet.get(Length.class, e, Double.class);
			numIncident++;
			if (HalfEdgeUtils.isBoundaryEdge(e)) {
				if (v1 == null) {
					v1 = e.getStartVertex();
				} else {
					v2 = e.getStartVertex();
				}
			} else {
				v3 = e.getStartVertex();
			}
		}
		// get cut border
		Set<CoEdge> cutEdges = new HashSet<CoEdge>();
		for (CoFace f : HalfEdgeUtils.facesIncidentWithVertex(v0)) {
			cutEdges.addAll(HalfEdgeUtils.boundaryEdges(f));
		}
		Set<CoEdge> inEdges = new HashSet<CoEdge>(incomingEdges(v0));
		Set<CoEdge> outEdges = new HashSet<CoEdge>(outgoingEdges(v0));
		cutEdges.removeAll(inEdges);
		cutEdges.removeAll(outEdges);
		
		// change edge lengths conformally such that the edges incident with vertex 0 are of equal length
		meanLength /= numIncident;
		Map<CoVertex, Double> uMap = new HashMap<CoVertex, Double>();
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			double l = aSet.get(Length.class, e, Double.class);
			double u = meanLength / l;
			uMap.put(e.getStartVertex(), u);
		}
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			CoVertex v = e.getStartVertex();
			for (CoEdge ee : HalfEdgeUtils.incomingEdges(v)) {
				double l = aSet.get(Length.class, ee, Double.class);
				CoVertex vs = ee.getStartVertex();
				CoVertex vt = ee.getTargetVertex();
				double us = uMap.containsKey(vs) ? uMap.get(vs) : 1.0;
				double ut = uMap.containsKey(vt) ? uMap.get(vt) : 1.0;
				double nl = l * us * ut;
				lMap.put(ee, nl);
				lMap.put(ee.getOppositeEdge(), nl);
			}
		}
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			CoVertex s = e.getStartVertex();
			CustomVertexInfo info = new CustomVertexInfo();
			info.boundaryMode = BoundaryMode.Isometric;
			s.info = info;
		}
		MappedEdgeLengthAdapter eAdapter = new MappedEdgeLengthAdapter(lMap, 100);
		aSet.add(eAdapter);
		
		// punch out the first vertex
		TopologyAlgorithms.removeVertex(v0);
		Map<CoEdge, CoEdge> nextOnCutBoundaryMap = new HashMap<CoEdge, CoEdge>();
		for (CoEdge e : cutEdges) {
			CoEdge next = e.getNextEdge();
			if (cutEdges.contains(next)) {
				nextOnCutBoundaryMap.put(e, next);
			}
		}
		
		// map to the upper half space
		CEuclideanOptimizable opt = new CEuclideanOptimizable(surface);
		int n = prepareInvariantDataEuclidean(opt.getFunctional(), surface, BoundaryMode.Conformal, QuantizationMode.Straight, aSet);
		CEuclideanApplication app = new CEuclideanApplication(surface);
		Vec u = new Vec(n);
		u.assemble();
		for (CoEdge e : surface.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				u.setValue(e.getSolverIndex(), e.getLambda(), InsertMode.INSERT_VALUES);
			}
		}
		app.setInitialSolutionVec(u);
		Mat H = app.getHessianTemplate();
		app.setHessianMat(H, H);
		
		Tao optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(gradTolerance, gradTolerance, gradTolerance); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(maxIterations);
		System.out.println("Using grad tolerance " + gradTolerance);
		optimizer.solve();
		if (optimizer.getSolutionStatus().reason != ConvergenceFlags.CONVERGED_ATOL) {
			throw new RuntimeException("Optimizer did not converge: \n" + optimizer.getSolutionStatus());
		}
		System.out.println(optimizer.getSolutionStatus());
		double[] uValues = u.getArray();
		DenseVector uVec = new DenseVector(uValues);
		layoutRoot = EuclideanLayout.doLayout(surface, opt.getFunctional(), uVec);
		
		
		// pre-normalize before projection
		double[] v1T = Pn.dehomogenize(v1.T, v1.T);
		double[] v2T = Pn.dehomogenize(v2.T, v2.T);
		double[] v3T = Pn.dehomogenize(v3.T, v3.T);
		double dist = Rn.euclideanDistance(v1T, v2T);
		double angle = Math.atan2(v1T[1] - v2T[1], v1T[0] - v2T[0]);
		MatrixBuilder mb = MatrixBuilder.euclidean();
		mb.scale(2 / dist);
		mb.rotate(-angle, 0, 0, 1);
		mb.translateFromTo(v1T, new double[]{0,0,0,1});
		Matrix M = mb.getMatrix();
		double[] check = M.multiplyVector(v3T);
		MatrixBuilder nb = MatrixBuilder.euclidean();
		nb.translate(check[0] < 0 ? 1 : -1, 0, 0);
		Matrix N = nb.getMatrix();
		M = Matrix.times(N, M);
		nb = MatrixBuilder.euclidean().scale(1, check[1] > 0 ? -1 : 1, 1);
		N = nb.getMatrix();
		M = Matrix.times(N, M);
		double[] baryCenter = new double[2];
		double[] tmp = new double[2];
		for (CoVertex v : surface.getVertices()) {
			M.transformVector(v.T);
			Pn.dehomogenize(v.T, v.T);
			tmp[0] = v.T[0];
			tmp[1] = v.T[1];
			Rn.add(baryCenter, baryCenter, tmp);
		}
		Rn.times(baryCenter, 1.0 / surface.numVertices(), baryCenter);
		nb = MatrixBuilder.euclidean();
		nb.translate(-baryCenter[0], 0, 0);
		N = nb.getMatrix();
		
		// project stereographically
		double[] spherePoint = new double[3];
		for (CoVertex v : surface.getVertices()) {
			N.transformVector(v.T);
			Pn.dehomogenize(v.T, v.T);
			Complex hp = new Complex(v.T[0], v.T[1]);
			hp = hp.times(10);
			double[] sp = ComplexUtility.inverseStereographic(hp, spherePoint);
			double tmp1 = sp[1];
			sp[1] = sp[2];
			sp[2] = tmp1;
			v.T[0] = sp[0];
			v.T[1] = sp[1];
			v.T[2] = sp[2];
			v.T[3] = 1.0;
		}
		
		// reinsert vertex 0
		CoVertex v0New = surface.addNewVertex();
		v0New.P = v0.P;
		v0New.T = new double[] {0,1,0,1};
		CoEdge onBoundaryIn = null;
		CoEdge onBoundaryOut = null;
		CoEdge inNextOpp = null;
		CoEdge outPrevOpp = null;
		for (CoEdge e : cutEdges) {
			CoEdge in = surface.addNewEdge();
			CoEdge out = surface.addNewEdge();
			CoFace face = surface.addNewFace();
			if (!cutEdges.contains(e.getNextEdge()) && e.getNextEdge() != null) {
				onBoundaryOut = e.getNextEdge();
				outPrevOpp = in;
			}
			if (!cutEdges.contains(e.getPreviousEdge()) && e.getPreviousEdge() != null) {
				onBoundaryIn = e.getPreviousEdge();
				inNextOpp = out;
			}
			in.setLeftFace(face);
			out.setLeftFace(face);
			e.setLeftFace(face);
			in.setTargetVertex(v0New);
			out.setTargetVertex(e.getStartVertex());
			in.linkNextEdge(out);
			out.linkNextEdge(e);
			e.linkNextEdge(in);
		}
		for (CoEdge e : nextOnCutBoundaryMap.keySet()) {
			CoEdge next = nextOnCutBoundaryMap.get(e);
			e.getNextEdge().linkOppositeEdge(next.getPreviousEdge());
		}
		CoEdge inNext = surface.addNewEdge();
		CoEdge outPrev = surface.addNewEdge();
		inNext.setTargetVertex(v0New);
		outPrev.setTargetVertex(onBoundaryOut.getStartVertex());
		onBoundaryIn.linkNextEdge(inNext);
		onBoundaryOut.linkPreviousEdge(outPrev);
		inNext.linkOppositeEdge(inNextOpp);
		outPrev.linkOppositeEdge(outPrevOpp);
		inNext.linkNextEdge(outPrev);
		
		// normalize
		try {
			List<CoVertex> bVerts = new LinkedList<CoVertex>(boundaryVertices(surface));
			SphericalNormalizerPETc.normalize(surface, bVerts);
		} catch (NotConvergentException e) {
			log.info("Sphere normalization did not succeed: " + e.getMessage());
		}
		
		// project back
		for (CoVertex v : surface.getVertices()) {
			Pn.dehomogenize(v.T, v.T);
			Complex cp = ComplexUtility.stereographic(v.T);
			v.T[0] = cp.re;
			v.T[1] = cp.im;
			v.T[2] = 0;
			v.T[3] = 1;	
		}
	}

	
	/**
	 * Project stereographically onto the sphere
	 * @param graph
	 * @param scale
	 */
	public static void inverseStereographicProjection(CoHDS hds, double scale){
		for (CoVertex v : hds.getVertices()){
			Pn.dehomogenize(v.T, v.T);
			double x = v.T[0] / scale;
			double y = v.T[1] / scale;
			double nx = 2 * x;
			double ny = x*x + y*y - 1;
			double nz = 2 * y;
			double nw = ny + 2;
			v.T = new double[] {nx / nw, ny / nw, nz / nw, 1.0};
		}
	}
	
	
	public static double[] baryCenter(CoHDS hds){
		double[] result = {0,0,0,1};
		for (CoVertex v : hds.getVertices()){
			Rn.add(result, v.T, result);
		}
		return Rn.times(result, 1.0 / hds.numVertices(), result);
	}
	
	
	public static double meanRadius(CoHDS hds){
		double result = 0;
		double[] zero = {0,0,0,1};
		for (CoVertex v : hds.getVertices()){
			result += Pn.distanceBetween(zero, v.T, Pn.EUCLIDEAN);
		}
		return result / hds.numVertices();
	}
	
	public static void normalizeBeforeProjection(CoHDS hds, double scale){
		double[] offset = baryCenter(hds);
		Pn.dehomogenize(offset, offset);
		for (CoVertex v : hds.getVertices()){
			Rn.subtract(v.T, v.T, offset);
			v.T[3] = 1.0;
		}
		scale = meanRadius(hds) / scale;
		for (CoVertex v : hds.getVertices()){
			v.T[3] *= scale;
		}
	
	}


	@Override
	public void setGradientTolerance(double tol) {
		gradTolerance = tol;
	}

	@Override
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}

	@Override
	public void setCutRoot(CoVertex root) {
		unwrapRoot = root;
	}

	@Override
	public CuttingInfo<CoVertex, CoEdge, CoFace> getCutInfo() {
		return new CuttingInfo<CoVertex, CoEdge, CoFace>();
	}
	@Override
	public Map<CoEdge, Double> getlengthMap() {
		return null;
	}
	@Override
	public CoVertex getLayoutRoot() {
		return layoutRoot;
	}
	
}
