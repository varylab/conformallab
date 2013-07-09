package de.varylab.discreteconformal.unwrapper;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
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
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.GetSolutionStatusResult;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomEdgeInfo;
import de.varylab.discreteconformal.heds.CustomVertexInfo;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.Search.DefaultWeightAdapter;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class EuclideanUnwrapperPETSc implements Unwrapper {

	private QuantizationMode
		conesMode = QuantizationMode.AllAngles,
		boundaryQuantMode = QuantizationMode.AllAngles;
	private BoundaryMode
		boundaryMode = BoundaryMode.Isometric;
	private int 
		maxIterations = 150,
		numCones = 0;
	private double
		gradTolerance = 1E-8;
	private boolean
		cutAndLayout = true;

	public static double
		lastGNorm = 0;
	private Collection<CoVertex>
		cones = new HashSet<CoVertex>();
	
	public CoVertex
		layoutRoot = null;
	public CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	public Map<CoEdge, Double>
		lengthMap = null;
	private Vector
		uVec = null;	
	private CEuclideanApplication
		app = null;
	
	static {
		Tao.Initialize();		
	}
	
	public EuclideanUnwrapperPETSc() {
		super();
	}

	public EuclideanUnwrapperPETSc(boolean cutAndLayout) {
		super();
		this.cutAndLayout = cutAndLayout;
	}

	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		app = new CEuclideanApplication(surface);
		double[] uValues = calculateConformalFactors(surface, aSet, app);
		uVec = new DenseVector(uValues);

		if (cutAndLayout) {
			switch (genus) {
			case 0:
				cutInfo = ConesUtility.cutMesh(surface);
				break;
			case 1:
				CoVertex cutRoot = surface.getVertex(0);
				DefaultWeightAdapter<CoEdge> constantWeight = new DefaultWeightAdapter<CoEdge>();
				cutInfo = CuttingUtility.cutTorusToDisk(surface, cutRoot, constantWeight);
				break;
			default:
				throw new RuntimeException("Cannot work on higher genus");
			}
			lengthMap = EuclideanLayout.getLengthMap(surface, app.getFunctional(), uVec);
			layoutRoot = EuclideanLayout.doLayout(surface, app.getFunctional(), uVec);
		}
	}

	
	public Map<CoVertex, Double> calculateConformalFactors(CoHDS surface, AdapterSet aSet) throws UnwrapException {
		CEuclideanApplication app = new CEuclideanApplication(surface);
		double[] u = calculateConformalFactors(surface, aSet, app);
		DenseVector uVec = new DenseVector(u);
		MTJDomain uDomain = new MTJDomain(uVec);
		Map<CoVertex, Double> result = new LinkedHashMap<CoVertex, Double>();
		for (CoVertex v : surface.getVertices()) {
			Double uVal = app.getFunctional().getVertexU(v, uDomain);
			result.put(v, uVal);
		}
		return result;
	}
	
	
	private double[] calculateConformalFactors(CoHDS surface, AdapterSet aSet, CEuclideanApplication app) throws UnwrapException {
		UnwrapUtility.prepareInvariantDataEuclidean(app.getFunctional(), surface, boundaryMode, boundaryQuantMode, aSet);
		// cones
		cones = ConesUtility.setUpCones(surface, numCones); 
		// optimization
		Vec u;
		Tao optimizer;
		int n = app.getDomainDimension();
		u = new Vec(n);
		// set variable lambda start values
		for (CoEdge e : surface.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				u.setValue(e.getSolverIndex(), e.getLambda(), InsertMode.INSERT_VALUES);
			}
		}
		app.setInitialSolutionVec(u);
		Mat H = app.getHessianTemplate();
		app.setHessianMat(H, H);

		optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setGradientTolerances(gradTolerance, gradTolerance, gradTolerance);
		optimizer.setMaximumIterates(maxIterations);
		optimizer.solve();
		
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		lastGNorm = status.gnorm;
		if (status.reason.cvalue() < 0) {
			throw new UnwrapException("Optimization did not succeed: " + status);
		}

		if (!cones.isEmpty()) {
			if (conesMode != QuantizationMode.AllAngles) {
				// calculating cones
				cones = ConesUtility.quantizeCones(surface, cones, conesMode);
				
				// optimizing conformal structure
				CEuclideanApplication app2 = new CEuclideanApplication(surface);
				n = app2.getDomainDimension();
	
				u = new Vec(n);
				H = app.getHessianTemplate();
				
				app2.setInitialSolutionVec(u);
				app2.setHessianMat(H, H);
				
				optimizer.setApplication(app2);
				optimizer.setGradientTolerances(1E-5, 0, 0);
				optimizer.solve();
				
				status = optimizer.getSolutionStatus();
				System.out.println("Cone Quantization: " + status);
				if (status.reason.cvalue() < 0) {
					throw new UnwrapException("Cone quantization did not succeed: " + status);
				}
			}
		}
		double [] uValues = u.getArray();
		u.restoreArray();
		return uValues;
	}
	
	public Collection<CoVertex> getCones() {
		return cones;
	}
	public void setNumCones(int numCones) {
		this.numCones = numCones;
	}
	public void setConeMode(QuantizationMode quantizationMode) {
		this.conesMode = quantizationMode;
	}
	@Override
	public void setGradientTolerance(double tol) {
		gradTolerance = tol;
	}
	@Override
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}
	public void setBoundaryQuantMode(QuantizationMode boundaryQuantMode) {
		this.boundaryQuantMode = boundaryQuantMode;
	}
	public void setBoundaryMode(BoundaryMode boundaryMode) {
		this.boundaryMode = boundaryMode;
	}
	@Override
	public void setCutRoot(CoVertex root) {
	}
	
	public Vector getUResult() {
		return uVec;
	}
	public boolean isCutAndLayout() {
		return cutAndLayout;
	}
	public void setCutAndLayout(boolean cutAndLayout) {
		this.cutAndLayout = cutAndLayout;
	}

	@Override
	public CuttingInfo<CoVertex, CoEdge, CoFace> getCutInfo() {
		return cutInfo;
	}
	@Override
	public Map<CoEdge, Double> getlengthMap() {
		return lengthMap;
	}
	@Override
	public CoVertex getLayoutRoot() {
		return layoutRoot;
	}
	
	public ConformalFunctional<CoVertex, CoEdge, CoFace> getFunctional() {
		return app.getFunctional();
	}
	
	public static void unwrapcg(IndexedFaceSet ifs, int layoutRoot, Map<Integer, Double> boundaryAngles, List<Integer> circularVerts) {

		IndexedFaceSetUtility.makeConsistentOrientation(ifs);
		CoHDS hds = new CoHDS();
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		ConverterJR2Heds cTo = new ConverterJR2Heds();
		cTo.ifs2heds(ifs, hds, a);
		CircleDomainUnwrapper.flipAtEars(hds, a);
		
		// set boundary conditions
		for (CoVertex v : HalfEdgeUtils.boundaryVertices(hds)) {
			if (boundaryAngles.containsKey(v.getIndex())) {
				Double angle = boundaryAngles.get(v.getIndex());
				v.info = new CustomVertexInfo();
				v.info.useCustomTheta = true;
				v.info.theta = angle;
			}
		}
		
		// find circular edges
		List<CoEdge> circularEdges = new LinkedList<CoEdge>();
		
		if (circularVerts != null) {
			// find the boundary components associated to each vertex in the list of circular vertices
			for (int vIndex : circularVerts) {
				System.err.println("Processing circ vertex "+vIndex);
				CoVertex v1 = hds.getVertex(vIndex);
				List<CoEdge> edges = HalfEdgeUtils.incomingEdges(v1);
				CoEdge bedge = null;
				for (CoEdge e : edges)	{
					if (HalfEdgeUtils.isBoundaryEdge(e)) {
						bedge = e;
						break;
					}
				}
				if (bedge == null)	{
					throw new IllegalStateException("No boundary edge on vertex "+vIndex);
				}
				CoFace filler = HalfEdgeUtils.fillHole(bedge);
				circularEdges.addAll(Triangulator.triangulateFace(filler, hds));
			}
			for (CoEdge ce : circularEdges) {
//				CoVertex v2 = hds.getVertex(eIndex[1]);
//				CoEdge ce = HalfEdgeUtils.findEdgeBetweenVertices(v1, v2);
//				if (ce == null) throw new RuntimeException("could not find circular edge between vertex " + v1 + " and " + v2);
				CoEdge ceOpp = ce.getOppositeEdge(); 
				ce.info = new CustomEdgeInfo();
				ce.info.circularHoleEdge = true;
				ceOpp.info = new CustomEdgeInfo();
				ceOpp.info.circularHoleEdge = true;
			}
		}
		
		EuclideanUnwrapperPETSc unwrap = new EuclideanUnwrapperPETSc();
		unwrap.setBoundaryQuantMode(QuantizationMode.Straight);
		unwrap.setBoundaryMode(BoundaryMode.Conformal);
		unwrap.setGradientTolerance(1E-6);
		unwrap.setMaxIterations(100);
		try {
			unwrap.unwrap(hds, 0, a);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		double[][] texArr = new double[hds.numVertices()][];
		for (int i = 0; i < texArr.length; i++) {
			CoVertex v = hds.getVertex(i);
			texArr[i] = a.getD(TexturePosition4d.class, v);
		}

		
		CoVertex layoutRootVertex = hds.getVertex(layoutRoot);
		CoEdge bdIn = layoutRootVertex.getIncomingEdge();
		// find a boundary edge if there is any
		for (CoEdge e : HalfEdgeUtils.incomingEdges(layoutRootVertex)) {
			if (e.getLeftFace() == null) {
				bdIn = e;
				break;
			}
		}
		CoVertex nearVertex = bdIn.getStartVertex();
		int nearIndex = nearVertex.getIndex();
		
		double[] t = Pn.dehomogenize(null, texArr[layoutRoot]);
		double[] t2 = Pn.dehomogenize(null, texArr[nearIndex]);
		double angle = Math.atan2(t2[1] - t[1], t2[0] - t[0]);
		MatrixBuilder bT = MatrixBuilder.euclidean();
		bT.rotate(-angle, 0, 0, 1);
		bT.translate(-t[0], -t[1], 0);
		Matrix T = bT.getMatrix();
		for (int i = 0; i < texArr.length; i++) {
			double[] tc = texArr[i];
			double e = tc[3];
			Pn.dehomogenize(tc, tc);
			T.transformVector(tc);
			Rn.times(tc, e, tc);
		}
		
		DataList newTexCoords = new DoubleArrayArray.Array(texArr);
		ifs.setVertexAttributes(Attribute.TEXTURE_COORDINATES, null);
		ifs.setVertexAttributes(Attribute.TEXTURE_COORDINATES, newTexCoords);
	}
	
	
	public static void unwrap(IndexedFaceSet ifs, int layoutRoot, Map<Integer, Double> boundaryAngles, List<int[]> circularEdges) {
		IndexedFaceSetUtility.makeConsistentOrientation(ifs);
		CoHDS hds = new CoHDS();
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		ConverterJR2Heds cTo = new ConverterJR2Heds();
		cTo.ifs2heds(ifs, hds, a);
		CircleDomainUnwrapper.flipAtEars(hds, a);
		
		// set boundary conditions
		for (CoVertex v : HalfEdgeUtils.boundaryVertices(hds)) {
			if (boundaryAngles.containsKey(v.getIndex())) {
				Double angle = boundaryAngles.get(v.getIndex());
				v.info = new CustomVertexInfo();
				v.info.useCustomTheta = true;
				v.info.theta = angle;
			}
		}
		
		// find circular edges
		if (circularEdges != null) {
			for (int[] eIndex : circularEdges) {
				CoVertex v1 = hds.getVertex(eIndex[0]);
				CoVertex v2 = hds.getVertex(eIndex[1]);
				CoEdge ce = HalfEdgeUtils.findEdgeBetweenVertices(v1, v2);
				if (ce == null) throw new RuntimeException("could not find circular edge between vertex " + v1 + " and " + v2);
				CoEdge ceOpp = ce.getOppositeEdge(); 
				ce.info = new CustomEdgeInfo();
				ce.info.circularHoleEdge = true;
				ceOpp.info = new CustomEdgeInfo();
				ceOpp.info.circularHoleEdge = true;
			}
		}
		
		EuclideanUnwrapperPETSc unwrap = new EuclideanUnwrapperPETSc();
		unwrap.setBoundaryQuantMode(QuantizationMode.Straight);
		unwrap.setBoundaryMode(BoundaryMode.Conformal);
		unwrap.setGradientTolerance(1E-12);
		unwrap.setMaxIterations(50);
		try {
			unwrap.unwrap(hds, 0, a);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		double[][] texArr = new double[hds.numVertices()][];
		for (int i = 0; i < texArr.length; i++) {
			CoVertex v = hds.getVertex(i);
			texArr[i] = a.getD(TexturePosition4d.class, v);
		}

		
		CoVertex layoutRootVertex = hds.getVertex(layoutRoot);
		CoEdge bdIn = layoutRootVertex.getIncomingEdge();
		// find a boundary edge if there is any
		for (CoEdge e : HalfEdgeUtils.incomingEdges(layoutRootVertex)) {
			if (e.getLeftFace() == null) {
				bdIn = e;
				break;
			}
		}
		CoVertex nearVertex = bdIn.getStartVertex();
		int nearIndex = nearVertex.getIndex();
		
		double[] t = Pn.dehomogenize(null, texArr[layoutRoot]);
		double[] t2 = Pn.dehomogenize(null, texArr[nearIndex]);
		double angle = Math.atan2(t2[1] - t[1], t2[0] - t[0]);
		MatrixBuilder bT = MatrixBuilder.euclidean();
		bT.rotate(-angle, 0, 0, 1);
		bT.translate(-t[0], -t[1], 0);
		Matrix T = bT.getMatrix();
		for (int i = 0; i < texArr.length; i++) {
			double[] tc = texArr[i];
			double e = tc[3];
			Pn.dehomogenize(tc, tc);
			T.transformVector(tc);
			Rn.times(tc, e, tc);
		}
		
		DataList newTexCoords = new DoubleArrayArray.Array(texArr);
		ifs.setVertexAttributes(Attribute.TEXTURE_COORDINATES, null);
		ifs.setVertexAttributes(Attribute.TEXTURE_COORDINATES, newTexCoords);
	}
	
}
