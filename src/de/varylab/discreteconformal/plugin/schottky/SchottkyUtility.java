package de.varylab.discreteconformal.plugin.schottky;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position2d;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.math.ComplexUtility;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.HyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.unwrapper.Unwrapper;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.NodeIndexComparator;
import de.varylab.discreteconformal.util.SurgeryUtility;

public class SchottkyUtility {

	public static List<SchottkyCircle> getAllCircles(List<SchottkyGenerator> pairs) {
		List<SchottkyCircle> r = new LinkedList<SchottkyCircle>();
		for (SchottkyGenerator p : pairs) {
			SchottkyCircle sCircle = p.getCycle();
			SchottkyCircle tCircle = p.mapCircle(sCircle);
			r.add(sCircle);
			r.add(tCircle);
		}
		return r;
	}

	public static void cutIdentificationHoles(CoHDS hds, List<ArrayList<CoVertex>> circles) {
		for (ArrayList<CoVertex> circle : circles) {
			int len = circle.size();
			for (int i = 0; i < len; i++) {
				CoVertex pre = circle.get((i + len - 1) % len);
				CoVertex act = circle.get((i + len + 0) % len);
				CoVertex pos = circle.get((i + len + 1) % len);
				List<CoEdge> star = incomingEdges(act);
				for (CoEdge e : star) {
					if (!e.isValid()) continue;
					if (e.getStartVertex() == pre) continue;
					if (e.getStartVertex() == pos) continue;
					if (circle.contains(e.getStartVertex())) {
						TopologyAlgorithms.removeEdge(e);
					}
				}
			}
		}
	}

	public static CoVertex generateSurface(
		CoHDS hds, 
		List<SchottkyGenerator> pairs, 
		Complex rootPos, 
		Map<CoEdge, Double> lMap, 
		List<Set<CoEdge>> cyclesReturn, 
		Map<CoVertex, double[]> mapCycleMap,
		long randomSeed,
		int circleResolution,
		int numExtraPoints, 
		int numSpreadIterations
	) {
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		
		Map<CoVertex, CoVertex> sMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, CoVertex> sInvMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, SchottkyGenerator> vertexPairMap = new HashMap<CoVertex, SchottkyGenerator>(); 
		Map<CoVertex, SchottkyGenerator> vertexPairInvMap = new HashMap<CoVertex, SchottkyGenerator>(); 
		Map<CoEdge, CoEdge> edgeMap = new HashMap<CoEdge, CoEdge>();
		Map<CoEdge, SchottkyGenerator> edgePairMap = new HashMap<CoEdge, SchottkyGenerator>();
		Map<CoEdge, SchottkyGenerator> edgePairInvMap = new HashMap<CoEdge, SchottkyGenerator>();
		
		// equalized extra vertices
		Random rnd = new Random(randomSeed);
		for (int i = 0; i < numExtraPoints; i++) {
			CoVertex v = hds.addNewVertex();
			double[] vPos = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			Rn.normalize(vPos, vPos);
			a.set(Position.class, v, vPos);
		}
		if (numSpreadIterations > 0) {
			Set<CoVertex> emptySet = Collections.emptySet();
			SphereUtility.equalizeSphereVertices(hds, emptySet, numSpreadIterations, 1E-6);
		}
		for (CoVertex v : hds.getVertices()) {
			double[] vPos = a.getD(Position.class, v);
			Complex z = ComplexUtility.stereographic(vPos);
			a.set(Position.class, v, new double[] {z.re, z.im, 0});
		}
		
//		// the root vertex 
		CoVertex root = hds.addNewVertex();
		a.set(Position.class, root, new double[]{rootPos.re, rootPos.im, 0});
		
		// add the vertices on the source and target circles
		Set<CoVertex> fixedVertices = new HashSet<CoVertex>();
		List<ArrayList<CoVertex>> sourceVertexCircles = new ArrayList<ArrayList<CoVertex>>();
		List<ArrayList<CoVertex>> targetVertexCircles = new ArrayList<ArrayList<CoVertex>>();
		for (SchottkyGenerator p : pairs) {
			ArrayList<CoVertex> sourceCircle = new ArrayList<CoVertex>();
			ArrayList<CoVertex> targetCircle = new ArrayList<CoVertex>();
			sourceVertexCircles.add(sourceCircle);
			targetVertexCircles.add(targetCircle);
			for (int i = 0; i < circleResolution; i++) {
				SchottkyCircle c = p.getCycle();
				double phi = 2*i*PI / circleResolution;
				double x = c.getRadius() * cos(phi) + c.getCenter().re;
				double y = c.getRadius() * sin(phi) + c.getCenter().im;
				Complex z = new Complex(x, y);
				Complex sz = p.getMoebius().applyTo(z);
				double[] zPos = new double[] {z.re, z.im};
				double[] szPos = new double[] {sz.re, sz.im};
				CoVertex v = hds.addNewVertex();
				CoVertex sv = hds.addNewVertex();
				sourceCircle.add(v);
				targetCircle.add(sv);
				a.set(Position.class, v, zPos);
				a.set(Position.class, sv, szPos);
				sMap.put(v, sv);
				sInvMap.put(sv, v);
				vertexPairMap.put(v, p);
				vertexPairInvMap.put(sv, p);
				fixedVertices.add(v);
				fixedVertices.add(sv);
			}
		}

		// exclude extra vertices from the circles
		for (CoVertex v : new HashSet<CoVertex>(hds.getVertices())) {
			if (sMap.containsKey(v) || sMap.containsValue(v)) {
				continue; // we are a circle vertex
			}
			for (SchottkyCircle c : getAllCircles(pairs)) {
				double[] vp = a.getD(Position2d.class, v); 
				Complex vz = new Complex(vp[0], vp[1]);
				if (c.isInside(vz, c.getRadius() * 1E-4)) {
					if (!v.isValid()) {
						continue;
					}
					hds.removeVertex(v);
				}
			}
		}
		
		// scale and project 
		Map<CoVertex, Complex> zMap = new HashMap<CoVertex, Complex>();
		for (CoVertex v : hds.getVertices()) {
			double[] p = a.getD(Position2d.class, v);
			Complex z = new Complex(p[0], p[1]);
			zMap.put(v, z); // store original positions
			double[] pProjected = ComplexUtility.inverseStereographic(z);
			a.set(Position.class, v, pProjected);
		}
		
		// create triangulation and remove hole caps
		ConvexHull.convexHull(hds, a, 1E-8, true);
		List<ArrayList<CoVertex>> vertexCircles = new LinkedList<ArrayList<CoVertex>>();
		vertexCircles.addAll(sourceVertexCircles);
		vertexCircles.addAll(targetVertexCircles);
		cutIdentificationHoles(hds, vertexCircles);
		
		// back to C
		for (CoVertex v : hds.getVertices()) {
			Complex z = zMap.get(v);
			a.set(Position.class, v, new double[] {z.re, z.im});
		}
				
		// build edge identification opposites maps
		for (List<CoVertex> circle : sourceVertexCircles) {
			Set<CoVertex> vSet = new HashSet<CoVertex>(circle);
			for (CoVertex vt : circle) {
				CoVertex svt = sMap.get(vt);
				for (CoEdge e : incomingEdges(vt)) {
					CoVertex vs = e.getStartVertex();
					if (vSet.contains(vs)) {
						CoVertex svs = sMap.get(vs);
						for (CoEdge se : incomingEdges(svs)) {
							if (se.getStartVertex() == svt) {
								SchottkyGenerator p = vertexPairMap.get(vt);
								edgeMap.put(e, se);
								edgeMap.put(e.getOppositeEdge(), se.getOppositeEdge());
								edgeMap.put(se, e);
								edgeMap.put(se.getOppositeEdge(), e.getOppositeEdge());
								edgePairMap.put(e, p);
								edgePairMap.put(e.getOppositeEdge(), p);
								edgePairInvMap.put(se, p);
								edgePairInvMap.put(se.getOppositeEdge(), p);
							}
						}
					}
				}
			}
		}
		
		// we know how many edges there are on the circles
		assert edgeMap.size() == 2 * pairs.size() * circleResolution * 2 : "lengths of cycles";
		
		// calculate length cross ratios
		Map<CoEdge, Double> crMap = calculateLengthCrossRatiosFromEdgePairs(hds, edgeMap, edgePairMap, edgePairInvMap, zMap, a);
		
		// back to Chat again
		for (CoVertex v : hds.getVertices()) {
			double[] p = a.getD(Position2d.class, v);
			Complex z = new Complex(p[0], p[1]);
			double[] pProjected = ComplexUtility.inverseStereographic(z);
			a.set(Position.class, v, pProjected);
		}
		
		// store to be lost coordinates, also deleted vertices
		for (CoVertex v : sMap.keySet()) {
			CoVertex sv = sMap.get(v);
			double[] posv = a.getD(Position3d.class, v);
			double[] possv = a.getD(Position3d.class, sv);
			mapCycleMap.put(v, possv);
			mapCycleMap.put(sv, posv);
		}

		// identify circles
		try {
			Set<CoEdge> mappedEdgesSet = new TreeSet<CoEdge>(new NodeIndexComparator<CoEdge>());
			mappedEdgesSet.addAll(edgeMap.keySet());
			for (CoEdge e : mappedEdgesSet) {
				if (!e.isValid() || e.getLeftFace() != null) {
					continue;
				}
				CoEdge se = edgeMap.get(e);
				try {
					System.out.println("Identify " + e + " " + se);
					Set<CoEdge> cycle = SurgeryUtility.glueAlongBoundaries(e, se);
					cyclesReturn.add(cycle);
				} catch (RuntimeException re) {
					re.printStackTrace();
				}
			}
		} catch (Exception e) {
			System.err.println("Invalid Schottky configuration.");
		}
		
		// check cross-ratios
		checkCrossRatioAssignments(hds, edgeMap, crMap);
		
		// define lengths from cross-ratios
		calculateLengthsFromCrossRatios(hds, crMap, lMap);

		System.out.println("Generated surface of genus " + HalfEdgeUtils.getGenus(hds));
		return root.isValid() ? root : null;
	}

	protected static void checkCrossRatioAssignments(
		CoHDS hds,
		Map<CoEdge, CoEdge> edgeMap, 
		Map<CoEdge, Double> crMap
	) {
		// check if identification was correct
		for (CoEdge e : edgeMap.keySet()) {
			if (!e.isValid()) continue;
			CoEdge se = edgeMap.get(e);
			assert e == se.getOppositeEdge() : "opposite";
		}
		// assert mapped and opposite cross ratios are the same
		for (CoEdge e : edgeMap.keySet()) {
			CoEdge se = edgeMap.get(e);
			double cr = crMap.get(e);
			double scr = crMap.get(se);
			assert abs(cr - scr) < 1E-8 : "Mapped cross ratio: " + cr + " - " + scr;
		}
		for (CoEdge e : hds.getEdges()) {
			CoEdge opp = e.getOppositeEdge();
			double cr = crMap.get(e);
			double scr = crMap.get(opp);
			assert abs(cr - scr) < 1E-8 : "Opposite cross ratios: " + cr + " - " + scr;
		}
		// check cross-ratio product around each vertex
		for (CoVertex v : hds.getVertices()) {
			double pr = 1.0;
			for (CoEdge e : incomingEdges(v)) {
				pr *= crMap.get(e);
			}
			assert abs(pr - 1.0) < 1E-8 : "Cross-ratio product around vertex not 1.0";
		}
	}

	protected static void calculateLengthsFromCrossRatios(
		CoHDS hds,
		Map<CoEdge, Double> crMap, 
		Map<CoEdge, Double> lMap
	) {
		Map<CoEdge, Double> aMap = new HashMap<CoEdge, Double>();
		for (CoVertex v : hds.getVertices()) {
			double ai = 100.0;
			for (CoEdge e : incomingEdges(v)) {
				double q = crMap.get(e);
				ai *= q;
				aMap.put(e, ai);
			}
		}
		for (CoEdge e : hds.getEdges()) {
			double ai = aMap.get(e);
			double aj = aMap.get(e.getPreviousEdge());
			double l = sqrt(1 / (ai * aj));
			lMap.put(e, l);
		}
		
		// check opposite lengths
		for (CoEdge e : hds.getEdges()) {
			double l = lMap.get(e);
			double lo = lMap.get(e.getOppositeEdge());
			assert abs(l - lo) < 1E-8 : "lengths";
		}
	}

	protected static Map<CoEdge, Double> calculateLengthCrossRatiosFromEdgePairs(
		CoHDS hds, 
		Map<CoEdge, CoEdge> edgeMap,
		Map<CoEdge, SchottkyGenerator> edgePairMap,
		Map<CoEdge, SchottkyGenerator> edgePairInvMap,
		Map<CoVertex, Complex> zMap, 
		AdapterSet a
	) {
		Map<CoEdge, Double> crMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getEdges()) {
			CoVertex vi = e.getStartVertex();
			CoVertex vj = e.getTargetVertex();
			CoVertex vl = e.getNextEdge().getTargetVertex();
			CoVertex vk = e.getOppositeEdge().getNextEdge().getTargetVertex();
			Complex zi = zMap.get(vi);
			Complex zj = zMap.get(vj);
			Complex zl = zMap.get(vl);
			Complex zk = zMap.get(vk);
			
			Moebius pullTransform = null;
			// check if we have to pull a vertex position
			if (edgePairMap.containsKey(e)) {
				assert e.getLeftFace() == null || e.getRightFace() == null;
				// we are at a source circle edge
				SchottkyGenerator pair = edgePairMap.get(e);
				pullTransform = pair.getMoebius().invert();
			} 
			if (edgePairInvMap.containsKey(e)) {
				assert e.getLeftFace() == null || e.getRightFace() == null;
				// we are at a target circle edge
				SchottkyGenerator pair = edgePairInvMap.get(e);
				pullTransform = pair.getMoebius();
			}
			// pull one vertex from a different domain
			if (pullTransform != null) {
				CoEdge se = edgeMap.get(e);
				if (se.getLeftFace() == null) {
					CoVertex mappedL = se.getOppositeEdge().getNextEdge().getTargetVertex();
					zl = ComplexUtility.stereographic(a.getD(Position3d.class, mappedL));
					zl = pullTransform.applyTo(zl);
				} else {
					CoVertex mappedL = se.getNextEdge().getTargetVertex();
					zk = ComplexUtility.stereographic(a.getD(Position3d.class, mappedL));
					zk = pullTransform.applyTo(zk);
				}
			}
			double l_ik = zi.dist(zk); 
			double l_kj = zk.dist(zj); 
			double l_jl = zj.dist(zl); 
			double l_li = zl.dist(zi); 
			double cr = (l_kj * l_li) / (l_ik * l_jl);
			// store the cross ratio
			crMap.put(e, cr);
		}
		return crMap;
	}

	public static CuttingInfo<CoVertex, CoEdge, CoFace> unwrapSchottkySurface(
		CoHDS hds, 
		List<Set<CoEdge>> cycles,
		Map<CoVertex, double[]> mapCycleMap, 
		CoVertex rootVertex,
		AdapterSet aSet,
		boolean onSphere,
		boolean cutFromRoot
	) throws Exception {
		int genus = HalfEdgeUtils.getGenus(hds);
		Unwrapper unwrapper = null;
		if (genus > 1) {
			unwrapper = new HyperbolicUnwrapperPETSc(false);
		} else {
			unwrapper = new EuclideanUnwrapperPETSc(false);
		}
		unwrapper.setGradientTolerance(1E-8);
		unwrapper.setMaxIterations(2000);
		unwrapper.unwrap(hds, genus, aSet);
	
		// cut along schottky circles
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
		for (Set<CoEdge> cycle : cycles) {
			CuttingUtility.cutAlongPath(cycle, cutInfo);
			for (CoEdge ce : cycle) {
				CoVertex bv = ce.getStartVertex();
				CoVertex cv = cutInfo.vertexCopyMap.get(bv);
				double[] pos = mapCycleMap.get(bv);
				aSet.set(Position.class, cv, pos);
			}
		}
		
		if (!onSphere) {
			for (CoVertex v : hds.getVertices()) {
				double[] pos = aSet.getD(Position3d.class, v);
				Complex zPos = ComplexUtility.stereographic(pos);
				double[] newPos = {zPos.re, zPos.im};
				aSet.set(Position.class, v, newPos);
			}
		}
	
		// cut between circles
		if (cutFromRoot) {
			CuttingUtility.cutToSimplyConnected(hds, rootVertex, cutInfo);
		} else {
			CuttingUtility.cutToSimplyConnected(hds, null, cutInfo);
		}
		
		if (rootVertex == null) {
			rootVertex = hds.getVertex(0);
		}
		if (genus > 1) {
			HyperbolicUnwrapperPETSc hypUnwrapper = (HyperbolicUnwrapperPETSc)unwrapper;
			ConformalFunctional<CoVertex, CoEdge, CoFace> fun = hypUnwrapper.getFunctional();
			Vector u = hypUnwrapper.getUResult();
			rootVertex = HyperbolicLayout.doLayout(hds, null, fun,  u);
		} else {
			EuclideanUnwrapperPETSc eucUnwrapper = (EuclideanUnwrapperPETSc)unwrapper;
			ConformalFunctional<CoVertex, CoEdge, CoFace> fun = eucUnwrapper.getFunctional();
			Vector u = eucUnwrapper.getUResult();
			rootVertex = EuclideanLayout.doLayout(hds, fun, u);
		}
		return cutInfo;
	}

}
