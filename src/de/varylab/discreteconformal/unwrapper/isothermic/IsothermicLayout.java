package de.varylab.discreteconformal.unwrapper.isothermic;

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;

public class IsothermicLayout {

	/**
	 * Calculates the layout of a mesh from given edge angles
	 * @param hds
	 * @param angleMap
	 * @param a
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void  doTexLayout(HDS hds, Map<E, Double> angleMap, AdapterSet aSet) {
		Map<E, Double> lengthMap = new HashMap<E, Double>();
		Set<V> visited = new HashSet<V>(hds.numVertices());
		Queue<V> Qv = new LinkedList<V>();
		Queue<E> Qe = new LinkedList<E>();
		Queue<Double> Qa = new LinkedList<Double>();
		
		// start at boundary edge if there is one
		E e0 = hds.getEdge(0);
		E e1 = e0.getOppositeEdge();
		Double e0a = angleMap.get(e0);
		V v1 = e0.getStartVertex();
		V v2 = e0.getTargetVertex();
		// queued data
		Qv.offer(v1);
		Qv.offer(v2);
		Qe.offer(e1);
		Qe.offer(e0);
		Qa.offer(e0a + PI);
		Qa.offer(e0a);

		// vertices
		Double l = 2.0;
		lengthMap.put(e0, l);
		lengthMap.put(e1, l);
		aSet.set(TexturePosition.class, v1, new double[] {-cos(e0a), -sin(e0a), 0, 1});
		aSet.set(TexturePosition.class, v2, new double[] {cos(e0a), sin(e0a), 0, 1});		
		visited.add(v1);
		visited.add(v2);
		
		while (!Qv.isEmpty() && hds.numVertices() > visited.size()) {
			V v = Qv.poll();
			E inE = Qe.poll();
			Double a = Qa.poll();
			E outE = inE.getOppositeEdge();
			double[] tp = aSet.getD(TexturePosition4d.class, v);
			
			E e = inE.getNextEdge();
			Double globalAngle = a + PI;
			while (e != outE) {
				V nearVertex = e.getTargetVertex();
				
				E next = e.getNextEdge();
				Double alpha = IsothermicUtility.getOppositeBeta(next, angleMap);
				if (e.getLeftFace() == null) { // a boundary edge
					//TODO do layout in the other direction too
					break;
//					alpha = 2*PI - calculateAngleSum(v, angleMap);
				}
				
				globalAngle -= alpha;
				
				if (!lengthMap.containsKey(e)) {
					double prevLength = lengthMap.get(e.getPreviousEdge());
					l = IsothermicUtility.getEdgeLength(e, prevLength, angleMap);
					lengthMap.put(e, l);
					lengthMap.put(e.getOppositeEdge(), l);
				}
				
				if (!visited.contains(nearVertex)) {
					visited.add(nearVertex);
					Qv.offer(nearVertex);
					Qe.offer(e);
					Qa.offer(globalAngle);
					
					double[] dif = {cos(globalAngle), sin(globalAngle), 0.0, 1.0};
					Rn.times(dif, l, dif);
					double[] t = Rn.add(null, tp, dif);
					t[3] = 1.0;
					aSet.set(TexturePosition.class, nearVertex, t);
				}
				e = e.getOppositeEdge().getNextEdge();
			}
		}
		
		assert (visited.size() == hds.numVertices());
	}
	
	
}
