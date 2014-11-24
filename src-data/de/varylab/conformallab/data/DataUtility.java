package de.varylab.conformallab.data;

import java.lang.annotation.Annotation;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.logging.Logger;

import de.jreality.geometry.IndexedFaceSetFactory;
import de.jreality.math.Matrix;
import de.jreality.math.Pn;
import de.jreality.scene.IndexedFaceSet;
import de.jtem.discretegroup.core.DiscreteGroup;
import de.jtem.discretegroup.core.DiscreteGroupElement;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.riemann.surface.BranchPoint;
import de.varylab.conformallab.data.types.Circle;
import de.varylab.conformallab.data.types.Complex;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMap;
import de.varylab.conformallab.data.types.DiscreteMetric;
import de.varylab.conformallab.data.types.EdgeIdentification;
import de.varylab.conformallab.data.types.EmbeddedTriangle;
import de.varylab.conformallab.data.types.EmbeddedVertex;
import de.varylab.conformallab.data.types.EmbeddingSelection;
import de.varylab.conformallab.data.types.EmbeddingSelection.EdgeSelection;
import de.varylab.conformallab.data.types.EmbeddingSelection.FaceSelection;
import de.varylab.conformallab.data.types.EmbeddingSelection.VertexSelection;
import de.varylab.conformallab.data.types.FundamentalEdge;
import de.varylab.conformallab.data.types.FundamentalVertex;
import de.varylab.conformallab.data.types.HalfedgeEdge;
import de.varylab.conformallab.data.types.HalfedgeEmbedding;
import de.varylab.conformallab.data.types.HalfedgeFace;
import de.varylab.conformallab.data.types.HalfedgeMap;
import de.varylab.conformallab.data.types.HalfedgeSelection;
import de.varylab.conformallab.data.types.HalfedgeVertex;
import de.varylab.conformallab.data.types.HyperEllipticAlgebraicCurve;
import de.varylab.conformallab.data.types.IsometryPSL2R;
import de.varylab.conformallab.data.types.MetricEdge;
import de.varylab.conformallab.data.types.MetricTriangle;
import de.varylab.conformallab.data.types.ObjectFactory;
import de.varylab.conformallab.data.types.SchottkyData;
import de.varylab.conformallab.data.types.UniformizationData;
import de.varylab.conformallab.data.types.UniformizingGroup;
import de.varylab.conformallab.data.types.VertexIdentification;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.MappedEdgeLengthAdapter;
import de.varylab.discreteconformal.math.P2Big;
import de.varylab.discreteconformal.math.PnBig;
import de.varylab.discreteconformal.math.RnBig;
import de.varylab.discreteconformal.plugin.hyperelliptic.Curve;
import de.varylab.discreteconformal.plugin.schottky.SchottkyGenerator;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class DataUtility {

	private static Logger
		log = Logger.getLogger(DataUtility.class.getName());
	
	public static HalfedgeEmbedding toHalfedgeEmbedding(DiscreteEmbedding de) {
		CoHDS hds = new CoHDS();
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
		AdapterSet a = new ConformalAdapterSet();
		toHalfedge(de, a, Position.class, hds, cutInfo);
		return toHalfedgeEmbedding(de.getName(), hds, a, Position4d.class, cutInfo);
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> HalfedgeEmbedding toHalfedgeEmbedding(
		String name, 
		HDS hds, 
		AdapterSet a, 
		Class<? extends Annotation> type,
		CuttingInfo<V, E, F> identifications
	) {
		ObjectFactory of = new ObjectFactory();
		HalfedgeEmbedding he = of.createHalfedgeEmbedding();
		he.setName(name);
		for (V v : hds.getVertices()) {
			HalfedgeVertex hv = of.createHalfedgeVertex();
			double[] p = a.getD(type, v); 
			hv.setX(p[0]);
			hv.setY(p[1]);
			hv.setZ(p[2]);
			hv.setW(p[3]);
			hv.setIndex(v.getIndex());
			he.getVertices().add(hv);
		}
		for (E e : hds.getEdges()) {
			HalfedgeEdge hed = of.createHalfedgeEdge();
			hed.setIndex(e.getIndex());
			hed.setLeft(e.getLeftFace() == null ? -1 : e.getLeftFace().getIndex());
			hed.setOpposite(e.getOppositeEdge().getIndex());
			hed.setNext(e.getNextEdge().getIndex());
			hed.setTarget(e.getTargetVertex().getIndex());
			he.getEdges().add(hed);
		}
		for (F f : hds.getFaces()) {
			HalfedgeFace hf = of.createHalfedgeFace();
			hf.setIndex(f.getIndex());
			he.getFaces().add(hf);
		}
		if (identifications != null && !identifications.vertexCopyMap.isEmpty()) {
			Set<V> ready = new HashSet<V>();
			for (V v : HalfEdgeUtils.boundaryVertices(hds)) {
				if (ready.contains(v)) continue;
				Set<V> idSet = identifications.getCopies(v);
				VertexIdentification vi = new VertexIdentification();
				for (V iv : idSet) {
					vi.getVertices().add(iv.getIndex());
					ready.add(iv);
				}
				he.getIdentifications().add(vi);
			}
			Set<E> readyEdges = new HashSet<>();
			for (E e : HalfEdgeUtils.boundaryEdges(hds)) {
				if (readyEdges.contains(e)) continue;
				E eCopy = null;
				if (identifications.edgeCutMap.containsKey(e)) {
					eCopy = identifications.edgeCutMap.get(e);
				} else {
					e = e.getOppositeEdge();
					eCopy = identifications.edgeCutMap.get(e);
				}
				if (eCopy == null) {
					log.warning("illegal edge copy map, could not find copy of edge " + e);
					continue;
				}
				EdgeIdentification ei = of.createEdgeIdentification();
				ei.setEdge1(e.getIndex());
				ei.setEdge2(eCopy.getIndex());
				ei.setEdge3(e.getOppositeEdge().getIndex());
				ei.setEdge4(eCopy.getOppositeEdge().getIndex());
				he.getEdgeIdentifications().add(ei);
				readyEdges.add(e);
				readyEdges.add(eCopy);
				readyEdges.add(e.getOppositeEdge());
				readyEdges.add(eCopy.getOppositeEdge());
			}
		}
		return he;
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void toHalfedge(
		HalfedgeEmbedding he, 
		AdapterSet a,
		Class<? extends Annotation> type,
		HDS hdsOUT,
		CuttingInfo<V, E, F> cutInfoOUT
	) {
		Map<Integer, HalfedgeVertex> vMap = new HashMap<>();
		Map<Integer, HalfedgeEdge> eMap = new HashMap<>();
		Map<Integer, HalfedgeFace> fMap = new HashMap<>();
		for (HalfedgeVertex v : he.getVertices()) {
			vMap.put(v.getIndex(), v);
		}
		for (HalfedgeEdge e : he.getEdges()) {
			eMap.put(e.getIndex(), e);
		}
		for (HalfedgeFace f : he.getFaces()) {
			fMap.put(f.getIndex(), f);
		}
		hdsOUT.clear();
		hdsOUT.addNewVertices(he.getVertices().size());
		hdsOUT.addNewEdges(he.getEdges().size());
		hdsOUT.addNewFaces(he.getFaces().size());
		for (Integer i : vMap.keySet()) {
			HalfedgeVertex v = vMap.get(i);
			V vv = hdsOUT.getVertex(i);
			double[] p = {v.getX(), v.getY(), v.getZ(), v.getW()};
			a.set(type, vv, p);
		}
		for (Integer i : eMap.keySet()) {
			HalfedgeEdge e = eMap.get(i);
			E ee = hdsOUT.getEdge(i);
			V target = hdsOUT.getVertex(e.getTarget());
			E opposite = hdsOUT.getEdge(e.getOpposite());
			E next = hdsOUT.getEdge(e.getNext());
			F left = e.getLeft() < 0 ? null : hdsOUT.getFace(e.getLeft());
			ee.setTargetVertex(target);
			ee.linkOppositeEdge(opposite);
			ee.linkNextEdge(next);
			ee.setLeftFace(left);
		}
		if (he.getIdentifications().size() != 0 && cutInfoOUT != null) {
			createCuttingInfo(hdsOUT, he, cutInfoOUT);
		}
	}
	
	
	public static EmbeddingSelection toEmbeddingSelection(Selection s) {
		ObjectFactory of = new ObjectFactory();
		EmbeddingSelection ems = of.createEmbeddingSelection();
		VertexSelection vs = new VertexSelection();
		EdgeSelection es = new EdgeSelection();
		FaceSelection fs = new FaceSelection();
		ems.setVertexSelection(vs);
		ems.setEdgeSelection(es);
		ems.setFaceSelection(fs);
		for (Vertex<?, ?, ?> v : s.getVertices()) {
			de.varylab.conformallab.data.types.EmbeddingSelection.VertexSelection.Vertex vv = new de.varylab.conformallab.data.types.EmbeddingSelection.VertexSelection.Vertex();
			vv.setIndex(v.getIndex());
			vv.setChannel(s.getChannel(v));
			vs.getVertices().add(vv);
		}
		for (Edge<?, ?, ?> e : s.getEdges()) {
			de.varylab.conformallab.data.types.EmbeddingSelection.EdgeSelection.Edge ee = new de.varylab.conformallab.data.types.EmbeddingSelection.EdgeSelection.Edge();
			if (e.getLeftFace() != null) {
				ee.setFace(e.getLeftFace().getIndex());
			} else {
				ee.setFace(-1);
			}
			ee.setVertex1(e.getStartVertex().getIndex());
			ee.setVertex2(e.getTargetVertex().getIndex());
			ee.setChannel(s.getChannel(e));
			es.getEdges().add(ee);
		}
		for (Face<?, ?, ?> f : s.getFaces()) {
			de.varylab.conformallab.data.types.EmbeddingSelection.FaceSelection.Face ff = new de.varylab.conformallab.data.types.EmbeddingSelection.FaceSelection.Face();
			ff.setIndex(f.getIndex());
			ff.setChannel(s.getChannel(f));
			fs.getFaces().add(ff);
		}
		return ems;
	}
	
	public static HalfedgeSelection toHalfedgeSelection(Selection s) {
		ObjectFactory of = new ObjectFactory();
		HalfedgeSelection ems = of.createHalfedgeSelection();
		de.varylab.conformallab.data.types.HalfedgeSelection.VertexSelection vs = new de.varylab.conformallab.data.types.HalfedgeSelection.VertexSelection();
		de.varylab.conformallab.data.types.HalfedgeSelection.EdgeSelection es = new de.varylab.conformallab.data.types.HalfedgeSelection.EdgeSelection();
		de.varylab.conformallab.data.types.HalfedgeSelection.FaceSelection fs = new de.varylab.conformallab.data.types.HalfedgeSelection.FaceSelection();
		ems.setVertexSelection(vs);
		ems.setEdgeSelection(es);
		ems.setFaceSelection(fs);
		for (Vertex<?, ?, ?> v : s.getVertices()) {
			de.varylab.conformallab.data.types.HalfedgeSelection.VertexSelection.Vertex vv = new de.varylab.conformallab.data.types.HalfedgeSelection.VertexSelection.Vertex();
			vv.setIndex(v.getIndex());
			vv.setChannel(s.getChannel(v));
			vs.getVertices().add(vv);
		}
		for (Edge<?, ?, ?> e : s.getEdges()) {
			de.varylab.conformallab.data.types.HalfedgeSelection.EdgeSelection.Edge ee = new de.varylab.conformallab.data.types.HalfedgeSelection.EdgeSelection.Edge();
			ee.setIndex(e.getIndex());
			ee.setChannel(s.getChannel(e));
			es.getEdges().add(ee);
		}
		for (Face<?, ?, ?> f : s.getFaces()) {
			de.varylab.conformallab.data.types.HalfedgeSelection.FaceSelection.Face ff = new de.varylab.conformallab.data.types.HalfedgeSelection.FaceSelection.Face();
			ff.setIndex(f.getIndex());
			ff.setChannel(s.getChannel(f));
			fs.getFaces().add(ff);
		}
		return ems;
	}
	
	
	public static Selection toSelection(EmbeddingSelection es, CoHDS hds) {
		Selection s = new Selection();
		for (de.varylab.conformallab.data.types.EmbeddingSelection.VertexSelection.Vertex v : es.getVertexSelection().getVertices()) {
			CoVertex vv = hds.getVertex(v.getIndex());
			s.add(vv, v.getChannel());
		}
		for (de.varylab.conformallab.data.types.EmbeddingSelection.EdgeSelection.Edge e : es.getEdgeSelection().getEdges()) {
			List<CoEdge> candidates = null;
			if (e.getFace() >= 0) {
				CoFace f = hds.getFace(e.getFace());
				candidates = HalfEdgeUtils.boundaryEdges(f);
			} else {
				candidates = HalfEdgeUtils.boundaryEdges(hds);
			}
			boolean found = false;
			for (CoEdge be : candidates) {
				if (
					be.getTargetVertex().getIndex() == e.getVertex1() && be.getStartVertex().getIndex() == e.getVertex2() ||
					be.getTargetVertex().getIndex() == e.getVertex2() && be.getStartVertex().getIndex() == e.getVertex1() 
				) {
					s.add(be, e.getChannel());
					found = true;
					break;
				}
			}
			if (!found) {
				log.warning("edge selection not found: " + e);
			}
		}
		for (de.varylab.conformallab.data.types.EmbeddingSelection.FaceSelection.Face f : es.getFaceSelection().getFaces()) {
			CoFace ff = hds.getFace(f.getIndex());
			s.add(ff, f.getChannel());
		}
		return s;
	}
	
	public static Selection toSelection(HalfedgeSelection es, CoHDS hds) {
		Selection s = new Selection();
		for (de.varylab.conformallab.data.types.HalfedgeSelection.VertexSelection.Vertex v : es.getVertexSelection().getVertices()) {
			CoVertex vv = hds.getVertex(v.getIndex());
			s.add(vv, v.getChannel());
		}
		for (de.varylab.conformallab.data.types.HalfedgeSelection.EdgeSelection.Edge e : es.getEdgeSelection().getEdges()) {
			CoEdge ee = hds.getEdge(e.getIndex());
			s.add(ee, e.getChannel());
		}
		for (de.varylab.conformallab.data.types.HalfedgeSelection.FaceSelection.Face f : es.getFaceSelection().getFaces()) {
			CoFace ff = hds.getFace(f.getIndex());
			s.add(ff, f.getChannel());
		}
		return s;
	}
	
	
	public static UniformizationData toUniformizationData(String name, FundamentalPolygon p) {
		ObjectFactory of = new ObjectFactory();
		UniformizingGroup fg = of.createUniformizingGroup();
		DiscreteGroup G = p.getDiscreteGroup();
		for (DiscreteGroupElement s : G.getGenerators()) {
			IsometryPSL2R m = toIsometryPSL2R(s);
			fg.getGenerators().add(m);
		}
		de.varylab.conformallab.data.types.FundamentalPolygon P = of.createFundamentalPolygon();
		for (de.varylab.discreteconformal.uniformization.FundamentalEdge e : p.getEdges()) {
			FundamentalEdge fe = new FundamentalEdge();
			fe.setIndex(e.index);
			fe.setNextEdge(e.nextEdge.index);
			fe.setPreviousEdge(e.prevEdge.index);
			fe.setIdentifiedEdge(e.partner.index);
			fe.setStartVertex(e.start.index);
			Complex sp = of.createComplex();
			BigDecimal[] spBig = PnBig.dehomogenize(null, e.startPosition, FundamentalPolygonUtility.context); 
			sp.setRe(spBig[0].doubleValue());
			sp.setIm(spBig[1].doubleValue());
			fe.setStartPosition(sp);
			P.getEdges().add(fe);
		}
		for (de.varylab.discreteconformal.uniformization.FundamentalVertex v : p.getVertices()) {
			FundamentalVertex fv = new FundamentalVertex();
			fv.setIndex(v.index);
			P.getVertices().add(fv);
		}
		UniformizationData result = of.createUniformizationData();
		result.setName(name);
		result.setUniformizingGroup(fg);
		result.setFundamentalPolygon(P);
		return result;
	}


	public static IsometryPSL2R toIsometryPSL2R(DiscreteGroupElement s) {
		IsometryPSL2R m = new IsometryPSL2R();
		Matrix M = s.getMatrix();
		m.setM11(M.getEntry(0, 0));
		m.setM12(M.getEntry(0, 1));
		m.setM13(M.getEntry(0, 3));
		
		m.setM21(M.getEntry(1, 0));
		m.setM22(M.getEntry(1, 1));
		m.setM23(M.getEntry(1, 3));
		
		m.setM31(M.getEntry(3, 0));
		m.setM32(M.getEntry(3, 1));
		m.setM33(M.getEntry(3, 3));
		return m;
	}
	
	
	public static FundamentalPolygon toFundamentalPolygon(UniformizationData data) {
		int numVertices = data.getFundamentalPolygon().getVertices().size();
		int numEdges = data.getFundamentalPolygon().getEdges().size();
		de.varylab.discreteconformal.uniformization.FundamentalVertex[] vertices = new de.varylab.discreteconformal.uniformization.FundamentalVertex[numVertices];
		de.varylab.discreteconformal.uniformization.FundamentalEdge[] edges = new de.varylab.discreteconformal.uniformization.FundamentalEdge[numEdges];
		for (FundamentalVertex v : data.getFundamentalPolygon().getVertices()) {
			vertices[v.getIndex()] = new de.varylab.discreteconformal.uniformization.FundamentalVertex(v.getIndex());
		}
		for (FundamentalEdge e : data.getFundamentalPolygon().getEdges()) {
			edges[e.getIndex()] = new de.varylab.discreteconformal.uniformization.FundamentalEdge(e.getIndex());
		}
		// create combinatorial information and polygon
		for (FundamentalEdge e : data.getFundamentalPolygon().getEdges()) {
			de.varylab.discreteconformal.uniformization.FundamentalEdge fe = edges[e.getIndex()];
			fe.partner = edges[e.getIdentifiedEdge()];
			fe.nextEdge = edges[e.getNextEdge()];
			fe.prevEdge = edges[e.getPreviousEdge()];
			fe.sourceEdgeCount = 1;
			double[] sp = {e.getStartPosition().getRe(), e.getStartPosition().getIm(), 0.0, 1.0};
			fe.startPosition = RnBig.toBig(null, sp);
			fe.start = vertices[e.getStartVertex()];
			fe.prevEdge.end = fe.start;
		}
		FundamentalPolygon P = new FundamentalPolygon(Arrays.asList(edges));
		int genus = P.getGenus();
		int signature = 0;
		if (genus == 0) signature = Pn.ELLIPTIC;
		if (genus == 1) signature = Pn.EUCLIDEAN;
		if (genus > 1) signature = Pn.HYPERBOLIC;
		// create isometries
		for (de.varylab.discreteconformal.uniformization.FundamentalEdge e : P.getEdges()) {
			BigDecimal[] s1 = e.startPosition;
			BigDecimal[] t1 = e.nextEdge.startPosition;
			BigDecimal[] s2 = e.partner.startPosition;
			BigDecimal[] t2 = e.partner.nextEdge.startPosition;
			e.motionBig = P2Big.imbedMatrixP2InP3(null, P2Big.makeDirectIsometryFromFrames(null, 
				P2Big.projectP3ToP2(null, s1), 
				P2Big.projectP3ToP2(null, t1), 
				P2Big.projectP3ToP2(null, t2), 
				P2Big.projectP3ToP2(null, s2), 
				signature,
				FundamentalPolygonUtility.context
			));
			e.motion = new Matrix(RnBig.toDouble(null, e.motionBig));
			System.out.println(e.motion);
		}
		return P;
	}
	
	
	public static HyperEllipticAlgebraicCurve toHyperEllipticAlgebraicCurve(String name, Curve curve) {
		HyperEllipticAlgebraicCurve c = new HyperEllipticAlgebraicCurve();
		c.setName(name);
		for (BranchPoint bp : curve.getBranchPoints()) {
			Complex cbp = new Complex();
			cbp.setRe(bp.getRe());
			cbp.setIm(bp.getIm());
			c.getBranchPoints().add(cbp);
		}
		return c;
	}
	
	public static Curve toCurve(HyperEllipticAlgebraicCurve curve) {
		BranchPoint[] bpArray = new BranchPoint[curve.getBranchPoints().size()];
		int bpIndex = 0;
		for (Complex cbp : curve.getBranchPoints()) {
			BranchPoint bp = new BranchPoint(cbp.getRe(), cbp.getIm());
			bpArray[bpIndex++] = bp;
		}
		return new Curve(bpArray);
	}
	
	public static IndexedFaceSet toIndexedFaceSet(DiscreteEmbedding de) {
		int numVertices = de.getVertices().size();
		int numFaces = de.getTriangles().size();
	
		double[][] vData = new double[numVertices][];
		for (EmbeddedVertex ev : de.getVertices()) {
			vData[ev.getIndex()] = new double[]{ev.getX(), ev.getY(), ev.getZ(), ev.getW()};
		}
		int faceIndex = 0;
		int[][] fData = new int[numFaces][];
		for (EmbeddedTriangle et : de.getTriangles()) {
			fData[faceIndex++] = new int[]{et.getVertex1(), et.getVertex2(), et.getVertex3()};
		}
		
		IndexedFaceSetFactory f = new IndexedFaceSetFactory();
		f.setVertexCount(numVertices);
		f.setFaceCount(numFaces);
		f.setVertexCoordinates(vData);
		f.setFaceIndices(fData);
		f.setGenerateEdgesFromFaces(true);
		f.setGenerateFaceNormals(true);
		f.setGenerateVertexNormals(true);
		f.update();
		IndexedFaceSet ifs = f.getIndexedFaceSet();
		return ifs;
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void toHalfedge(
		DiscreteEmbedding de, 
		AdapterSet a, 
		Class<? extends Annotation> type,
		HDS hdsOUT,
		CuttingInfo<V, E, F> cutInfoOUT
	) {
		hdsOUT.clear();
		IndexedFaceSet ifs = toIndexedFaceSet(de);
		ConverterJR2Heds c = new ConverterJR2Heds();
		c.ifs2heds(ifs, hdsOUT, a);
		for (V v : hdsOUT.getVertices()) {
			EmbeddedVertex vv = de.getVertices().get(v.getIndex());
			double[] p = {vv.getX(), vv.getY(), vv.getZ(), vv.getW()};
			a.set(type, v, p);
		}
		if (de.getIdentifications().size() != 0) {
			createCuttingInfo(hdsOUT, de, cutInfoOUT);
		}
	}
	
	
	public static int calculateGenus(DiscreteEmbedding de) {
		int v = de.getVertices().size();
		int f = de.getTriangles().size();
		int e = f * 3 / 2;
		// compensate for vertex identifications
		for (VertexIdentification vi : de.getIdentifications()) {
			v -= vi.getVertices().size() - 1;
		}
		int X = v - e + f;
		return (2 - X) / 2;
	}
	
	public static int calculateGenus(HalfedgeEmbedding de) {
		int v = de.getVertices().size();
		int f = de.getFaces().size();
		int e = de.getEdges().size() / 2;
		// compensate for vertex identifications
		for (VertexIdentification vi : de.getIdentifications()) {
			v -= vi.getVertices().size() - 1;
		}
		e -= de.getEdgeIdentifications().size();
		int X = v - e + f;
		return (2 - X) / 2;
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void createCuttingInfo(HDS hds, DiscreteEmbedding de, CuttingInfo<V, E, F> cutInfoOUT) {
		List<E> boundary = HalfEdgeUtils.boundaryEdges(hds);
		if (boundary.size() == 0) return;
		cutInfoOUT.vertexCopyMap.clear();
		for (VertexIdentification vi : de.getIdentifications()) {
			V copyRoot = hds.getVertex(vi.getVertices().get(0));
			for (int i : vi.getVertices()) {
				V copy = hds.getVertex(i);
				if (copy == copyRoot) continue;
				cutInfoOUT.vertexCopyMap.put(copyRoot, copy);
				copyRoot = copy;
			}
		}
		cutInfoOUT.edgeCutMap.clear();
		Set<E> readyEdges = new HashSet<>();
		Map<V, VertexIdentification> idMap = new HashMap<V, VertexIdentification>();
		for (VertexIdentification vi : de.getIdentifications()) {
			for (int iv : vi.getVertices()) {
				V v = hds.getVertex(iv);
				idMap.put(v, vi);
			}
		}
		for (E e : HalfEdgeUtils.boundaryEdges(hds)) {
			if (readyEdges.contains(e)) continue;
			V vs = e.getStartVertex();
			V vt = e.getTargetVertex();
			VertexIdentification si = idMap.get(vs);
			VertexIdentification ti = idMap.get(vt);
			Set<Integer> sIdSet = new HashSet<Integer>(si.getVertices());
			Set<Integer> tIdSet = new HashSet<Integer>(ti.getVertices());
			E copyEdge = null;
			for (int i : sIdSet) {
				if (vs.getIndex() == i) continue;
				V sCopy = hds.getVertex(i);
				for (E ee : HalfEdgeUtils.incomingEdges(sCopy)) {
					if (tIdSet.contains(ee.getStartVertex().getIndex())) {
						copyEdge = ee;
						break;
					}
				}
				if (copyEdge != null) break;
			}
			assert copyEdge != null;
			cutInfoOUT.edgeCutMap.put(e, copyEdge);
			cutInfoOUT.edgeCutMap.put(copyEdge, e);
			cutInfoOUT.edgeCutMap.put(e.getOppositeEdge(), copyEdge.getOppositeEdge());
			cutInfoOUT.edgeCutMap.put(copyEdge.getOppositeEdge(), e.getOppositeEdge());
			readyEdges.add(e);
			readyEdges.add(copyEdge);
		}
		cutInfoOUT.cutRoot = HalfEdgeUtils.boundaryVertices(hds).iterator().next();
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void createCuttingInfo(HDS hds, HalfedgeEmbedding he, CuttingInfo<V, E, F> cutInfoOUT) {
		List<E> boundary = HalfEdgeUtils.boundaryEdges(hds);
		if (boundary.size() == 0) return;
		cutInfoOUT.vertexCopyMap.clear();
		for (VertexIdentification vi : he.getIdentifications()) {
			V copyRoot = hds.getVertex(vi.getVertices().get(0));
			V lastCopy = copyRoot;
			for (int i : vi.getVertices()) {
				V copy = hds.getVertex(i);
				if (copy == copyRoot) continue;
				cutInfoOUT.vertexCopyMap.put(lastCopy, copy);
				lastCopy = copy;
			}
		}
		cutInfoOUT.edgeCutMap.clear();
		Set<E> readyEdges = new HashSet<>();
		Map<V, VertexIdentification> idMap = new HashMap<>();
		for (VertexIdentification vi : he.getIdentifications()) {
			for (int iv : vi.getVertices()) {
				V v = hds.getVertex(iv);
				idMap.put(v, vi);
			}
		}
		Map<E, EdgeIdentification> edgeIdMap = new HashMap<>();
		for (EdgeIdentification ei : he.getEdgeIdentifications()) {
			E e1 = hds.getEdge(ei.getEdge1());
			E e2 = hds.getEdge(ei.getEdge2());
			E e3 = hds.getEdge(ei.getEdge3());
			E e4 = hds.getEdge(ei.getEdge4());
			EdgeIdentification old1 = edgeIdMap.put(e1, ei);
			EdgeIdentification old2 = edgeIdMap.put(e2, ei);
			EdgeIdentification old3 = edgeIdMap.put(e3, ei);
			EdgeIdentification old4 = edgeIdMap.put(e4, ei);
			if (old1 != null || old2 != null || old3 != null || old4 != null) {
				log.warning("edge index occurring in more than one edge identification");
			}
		}
		for (E e : HalfEdgeUtils.boundaryEdges(hds)) {
			if (readyEdges.contains(e)) continue;
			EdgeIdentification ei = edgeIdMap.get(e);
			if (ei == null) continue;
			E e1 = hds.getEdge(ei.getEdge1());
			E e2 = hds.getEdge(ei.getEdge2());
			E e3 = hds.getEdge(ei.getEdge3());
			E e4 = hds.getEdge(ei.getEdge4());
			E be2 = null;
			be2 = e1.getLeftFace() == null && e1 != e ? e1 : be2;
			be2 = e2.getLeftFace() == null && e2 != e ? e2 : be2;
			be2 = e3.getLeftFace() == null && e3 != e ? e3 : be2;
			be2 = e4.getLeftFace() == null && e4 != e ? e4 : be2;
			if (be2 == null) {
				log.warning("could not find identified boundary edge");
				continue;
			}
			cutInfoOUT.edgeCutMap.put(e, be2);
			cutInfoOUT.edgeCutMap.put(be2, e);
			cutInfoOUT.edgeCutMap.put(e.getOppositeEdge(), be2.getOppositeEdge());
			cutInfoOUT.edgeCutMap.put(be2.getOppositeEdge(), e.getOppositeEdge());
			readyEdges.add(e1);
			readyEdges.add(e2);
			readyEdges.add(e3);
			readyEdges.add(e4);
		}
		int maxCopies = -1;
		for (V v : HalfEdgeUtils.boundaryVertices(hds)) {
			Set<V> copies = cutInfoOUT.getCopies(v);
			if (copies.size() > maxCopies) {
				maxCopies = copies.size();
				cutInfoOUT.cutRoot = v;
			}
		}
		
	}

	public static DiscreteMetric toDiscreteMetric(String name, CoHDS surface, AdapterSet a) {
		DiscreteMetric dm = new DiscreteMetric();
		dm.setName(name);
		int edgeIndex = 0;
		Map<Integer, Integer> edgeIndexMap = new HashMap<Integer, Integer>();
		for (CoEdge e : surface.getPositiveEdges()) {
			MetricEdge me = new MetricEdge();
			me.setIndex(edgeIndex);
			me.setLength(a.get(Length.class, e, Double.class));
			dm.getEdges().add(me);
			edgeIndexMap.put(e.getIndex(), edgeIndex);
			edgeIndex++;
		}
		for (CoFace f : surface.getFaces()) {
			MetricTriangle mt = new MetricTriangle();
			CoEdge e1 = f.getBoundaryEdge();
			CoEdge e2 = e1.getNextEdge();
			CoEdge e3 = e2.getNextEdge();
			e1 = e1.isPositive() ? e1 : e1.getOppositeEdge(); 
			e2 = e2.isPositive() ? e2 : e2.getOppositeEdge(); 
			e3 = e3.isPositive() ? e3 : e3.getOppositeEdge();
			mt.setEdge1(edgeIndexMap.get(e1.getIndex()));
			mt.setEdge2(edgeIndexMap.get(e2.getIndex()));
			mt.setEdge3(edgeIndexMap.get(e3.getIndex()));
			dm.getTriangles().add(mt);
		}
		return dm;
	}
	
	
	public static HalfedgeMap toHalfedgeMap(DiscreteMap dm) {
		HalfedgeEmbedding domain = toHalfedgeEmbedding(dm.getDomain());
		HalfedgeEmbedding image = toHalfedgeEmbedding(dm.getImage());
		ObjectFactory of = new ObjectFactory();
		HalfedgeMap map = of.createHalfedgeMap();
		map.setDomain(domain);
		map.setImage(image);
		map.setName(dm.getName());
		return map;
	}
	
	public static HalfedgeMap toHalfedgeMap(
		String name, 
		CoHDS surface, 
		AdapterSet a, 
		Class<? extends Annotation> domainType, 
		Class<? extends Annotation> imageType, 
		CuttingInfo<CoVertex, CoEdge, CoFace> identifications
	) {
		String domainName = name + "_domain";
		String imageName = name + "_image";
		HalfedgeEmbedding domain = toHalfedgeEmbedding(domainName, surface, a, domainType, identifications);
		HalfedgeEmbedding image = toHalfedgeEmbedding(imageName, surface, a, imageType, identifications);
		HalfedgeMap map = new HalfedgeMap();
		map.setDomain(domain);
		map.setImage(image);
		map.setName(name);
		return map;
	}
	
	public static DiscreteMap toDiscreteMap(
		String name, 
		CoHDS surface, 
		AdapterSet a, 
		Class<? extends Annotation> domainType, 
		Class<? extends Annotation> imageType, 
		CuttingInfo<CoVertex, CoEdge, CoFace> identifications
	) {
		String domainName = name + "_domain";
		String imageName = name + "_image";
		DiscreteEmbedding domain = toDiscreteEmbedding(domainName, surface, a, domainType, identifications);
		DiscreteEmbedding image = toDiscreteEmbedding(imageName, surface, a, imageType, identifications);
		DiscreteMap map = new DiscreteMap();
		map.setDomain(domain);
		map.setImage(image);
		map.setName(name);
		return map;
	}

	public static DiscreteEmbedding toDiscreteEmbedding(
		String name, 
		CoHDS surface, 
		AdapterSet a, 
		Class<? extends Annotation> type, 
		CuttingInfo<CoVertex, CoEdge, CoFace> identifications
	) {
		DiscreteEmbedding de = new DiscreteEmbedding();
		de.setName(name);
		for (CoVertex v : surface.getVertices()) {
			EmbeddedVertex ev = new EmbeddedVertex();
			double[] pos = a.getD(type, v);
			ev.setX(pos[0]);
			ev.setY(pos[1]);
			ev.setZ(pos[2]);
			ev.setW(pos[3]);
			ev.setIndex(v.getIndex());
			de.getVertices().add(ev);
		}
		for (CoFace f : surface.getFaces()) {
			List<CoEdge> b = HalfEdgeUtils.boundaryEdges(f);
			assert b.size() == 3 : "only triangulations are supported";
			CoVertex v0 = b.get(0).getStartVertex();
			CoVertex v1 = b.get(1).getStartVertex();
			CoVertex v2 = b.get(2).getStartVertex();
			EmbeddedTriangle et = new EmbeddedTriangle();
			et.setVertex1(v0.getIndex());
			et.setVertex2(v1.getIndex());
			et.setVertex3(v2.getIndex());
			de.getTriangles().add(et);
		}
		if (identifications != null && !identifications.vertexCopyMap.isEmpty()) {
			Set<CoVertex> ready = new HashSet<CoVertex>();
			for (CoVertex v : HalfEdgeUtils.boundaryVertices(surface)) {
				if (ready.contains(v)) continue;
				Set<CoVertex> idSet = identifications.getCopies(v);
				VertexIdentification vi = new VertexIdentification();
				for (CoVertex iv : idSet) {
					vi.getVertices().add(iv.getIndex());
					ready.add(iv);
				}
				de.getIdentifications().add(vi);
			}
		}
		return de;
	}

	public static CoHDS toHalfedgeAndLengths(DiscreteMetric dm, MappedEdgeLengthAdapter lengthAdapterOutput) {
		CoHDS hds = new CoHDS();
		Map<MetricEdge, CoEdge> edgeMap = new HashMap<MetricEdge, CoEdge>();
		for (MetricEdge me : dm.getEdges()) {
			CoEdge e1 = hds.addNewEdge();
			CoEdge e2 = hds.addNewEdge();
			e1.linkOppositeEdge(e2);
			edgeMap.put(me, e1);
		}
		Set<CoEdge> usedEdges = new HashSet<CoEdge>();
		for (MetricTriangle mt : dm.getTriangles()) {
			CoFace f = hds.addNewFace();
			MetricEdge e1 = dm.getEdges().get(mt.getEdge1());
			MetricEdge e2 = dm.getEdges().get(mt.getEdge2());
			MetricEdge e3 = dm.getEdges().get(mt.getEdge3());
			CoEdge ce1 = edgeMap.get(e1);
			CoEdge ce2 = edgeMap.get(e2);
			CoEdge ce3 = edgeMap.get(e3);
			if (usedEdges.contains(ce1)) {
				ce1 = ce1.getOppositeEdge();
			}
			if (usedEdges.contains(ce2)) {
				ce2 = ce2.getOppositeEdge();
			}
			if (usedEdges.contains(ce3)) {
				ce3 = ce3.getOppositeEdge();
			}
			usedEdges.add(ce1);
			usedEdges.add(ce2);
			usedEdges.add(ce3);
			ce1.linkNextEdge(ce2);
			ce2.linkNextEdge(ce3);
			ce3.linkNextEdge(ce1);
			ce1.setLeftFace(f);
			ce2.setLeftFace(f);
			ce3.setLeftFace(f);
		}
		// link boundary
		Set<CoEdge> boundaryEdges = new HashSet<CoEdge>(hds.getEdges());
		boundaryEdges.removeAll(usedEdges);
		for (CoEdge be : boundaryEdges) {
			CoEdge e = be.getOppositeEdge();
			CoEdge ne = e.getPreviousEdge().getOppositeEdge();
			while (ne.getLeftFace() != null) {
				ne = ne.getPreviousEdge().getOppositeEdge();
			}
			be.linkNextEdge(ne);
		}
		Queue<CoEdge> linkQueue = new LinkedList<CoEdge>(hds.getEdges());
		while (!linkQueue.isEmpty()) {
			CoEdge e = linkQueue.poll();
			if (e.getTargetVertex() != null) {
				continue;
			}
			CoVertex v = hds.addNewVertex();
			CoEdge ne = e;
			do {
				ne.setTargetVertex(v);
				ne = ne.getNextEdge().getOppositeEdge();
			} while (ne != e);
		}
		
		for (MetricEdge e : dm.getEdges()) {
			CoEdge e1 = edgeMap.get(e);
			CoEdge e2 = e1.getOppositeEdge();
			lengthAdapterOutput.setEdgeValue(e1, e.getLength(), null);
			lengthAdapterOutput.setEdgeValue(e2, e.getLength(), null);
		}
		return hds;
	}

	public static SchottkyData toSchottkyData(String name, List<SchottkyGenerator> generators) {
		ObjectFactory of = new ObjectFactory();
		SchottkyData schottkyData = of.createSchottkyData();
		schottkyData.setName(name);
		for (SchottkyGenerator s : generators) {
			de.varylab.conformallab.data.types.SchottkyGenerator g = of.createSchottkyGenerator();
			de.varylab.conformallab.data.types.Complex Ad = of.createComplex();
			de.varylab.conformallab.data.types.Complex Bd = of.createComplex();
			de.varylab.conformallab.data.types.Complex Mud = of.createComplex();
			de.varylab.conformallab.data.types.Complex center = of.createComplex();
			Circle Cd = of.createCircle();
			Ad.setRe(s.getA().re);
			Ad.setIm(s.getA().im);
			Bd.setRe(s.getB().re);
			Bd.setIm(s.getB().im);
			Mud.setRe(s.getMu().re);
			Mud.setIm(s.getMu().im);
			center.setRe(s.getCycle().getCenter().re);
			center.setIm(s.getCycle().getCenter().im);
			Cd.setCenter(center);
			Cd.setRadius(s.getCycle().getRadius());
			g.setA(Ad);
			g.setB(Bd);
			g.setMu(Mud);
			g.setCircle(Cd);
			schottkyData.getGenerators().add(g);
		}
		return schottkyData;
	}

	public static List<SchottkyGenerator> toSchottkyGeneratorsList(SchottkyData data) {
		List<SchottkyGenerator> result = new LinkedList<SchottkyGenerator>();
		for (de.varylab.conformallab.data.types.SchottkyGenerator s : data.getGenerators()) {
			SchottkyGenerator g = new SchottkyGenerator();
			g.getA().setRe(s.getA().getRe());
			g.getA().setIm(s.getA().getIm());
			g.getB().setRe(s.getB().getRe());
			g.getB().setIm(s.getB().getIm());
			g.getMu().setRe(s.getMu().getRe());
			g.getMu().setIm(s.getMu().getIm());
			g.getCycle().getCenter().setRe(s.getCircle().getCenter().getRe());
			g.getCycle().getCenter().setIm(s.getCircle().getCenter().getIm());
			g.getCycle().setRadius(s.getCircle().getRadius());
			result.add(g);
		}
		return result;
	}
	
}