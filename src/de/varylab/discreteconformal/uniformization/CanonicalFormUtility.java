package de.varylab.discreteconformal.uniformization;

import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;


public class CanonicalFormUtility {

	/**
	 * Input: Minimal fundamental polygon with one vertex and 4g sides
	 * @param P
	 * @return A new fundamental polygon where opposite sides are glued
	 */
	public static FundamentalPolygon canonicalizeOpposite(FundamentalPolygon P) {
		FundamentalPolygon R = FundamentalPolygonUtility.copyPolygon(P);
		
		List<FundamentalEdge> cluster = findLargestCluster(R);
		System.out.println("Largest Cluster #####: " + cluster);
		FundamentalEdge e = findExternalEdge(cluster, R);
		while (e != null) {
			System.out.println("moving " + e + " or " + e.partner + " to cluster");
			moveIntoCluster(e, cluster);
			e = findExternalEdge(cluster, R);
		}
		cluster = findLargestCluster(R);
		System.out.println("Largest Cluster #####: " + cluster);
		List<FundamentalEdge> coCluster = findCoCluster(cluster);
		System.out.println("Co-Cluster #####: " + coCluster);
		sortClusters(cluster, coCluster);
		cluster = findLargestCluster(R);
		System.out.println("After SORT: Largest Cluster #####: " + cluster);
		coCluster = findCoCluster(cluster);
		System.out.println("After SORT: Co-Cluster #####: " + coCluster);
		R.normalizeEdgeList();
		return R;
	}
	
	
	public static FundamentalEdge findExternalEdge(List<FundamentalEdge> cluster, FundamentalPolygon P) {
		for (FundamentalEdge e : P.getEdges()) {
			if (!cluster.contains(e) && !cluster.contains(e.partner)) {
				return e;
			}
		}
		return null;
	}
	
	
	public static void sortClusters(List<FundamentalEdge> C, List<FundamentalEdge> coC) {
		for (int i = 0; i < C.size() - 1; i++) {
			FundamentalEdge e = C.get(i);
			FundamentalEdge next = C.get(i + 1);			
			FundamentalEdge coE = e.partner;
			FundamentalEdge coNext = next.partner;
			int coIndex = coC.indexOf(coE);
			int coNextIndex = coC.indexOf(coNext);
			if (coIndex != coNextIndex - 1) {
				if (canMoveBehind(coNext, coE)) {
					moveBehind(coNext, coE);
					coC.remove(coNext);
					int insertionIndex = coC.indexOf(coE) + 1;
					coC.add(insertionIndex, coNext);
				} else {
					FundamentalEdge markerEdge = null;
					if (coIndex == coC.size() - 1) { 
						markerEdge = C.get(0);
					} else {
						markerEdge = coC.get(coIndex + 1);
					}
					if (canMoveInFront(coNext, markerEdge)) {
						moveInFront(coNext, markerEdge);
						coC.remove(coNext);
						int insertionIndex = coC.indexOf(coE) + 1;
						coC.add(insertionIndex, coNext);
					} else {
						throw new RuntimeException("unable to sort cluster");
					}
				}
			}
		}
	}
//	
//	
//	public static List<FundamentalEdge> findLongestSortedSubCluster(List<FundamentalEdge> C) {
//		
//		
//	}
//	
//	public static void moveIntoSortedSubCluster(FundamentalEdge e, List<FundamentalEdge> C) {
//		
//	}
	

	/**
	 * moves e or e.partner into the given cluster
	 * @param e
	 * @param cluster
	 */
	public static void moveIntoCluster(FundamentalEdge e, List<FundamentalEdge> cluster) {
		FundamentalEdge P = e;
		FundamentalEdge Q = null;
		int minCost = Integer.MAX_VALUE;
		for (FundamentalEdge active : cluster) {
			// TODO check if we can do better by using moveBehind
			if (canMoveInFront(P, active)) {
				int cost = moveInFrontCost(P, active); 
				if (minCost > cost) {
					Q = active;
					minCost = cost;
				}
			}
		}
		assert Q != null;
		if (Q == null) {
			throw new RuntimeException("Cannot move genarator into cluster.");
		}
		moveInFront(P, Q);
		int index = cluster.indexOf(Q);
		cluster.add(index, P);
	}
	
	
	public static boolean canMoveInFront(FundamentalEdge P, FundamentalEdge Q) {
		FundamentalEdge e = P;
		while (e != Q.prevEdge) {
			if (e == P.partner) {
				return false;
			}
			e = e.prevEdge;
		}
		return true;
	}
	public static boolean canMoveBehind(FundamentalEdge P, FundamentalEdge Q) {
		FundamentalEdge e = P.partner;
		while (e != Q.partner.nextEdge) {
			if (e == P) {
				return false;
			}
			e = e.nextEdge;
		}
		return true;
	}

	public static int moveInFrontCost(FundamentalEdge P, FundamentalEdge Q) {
		FundamentalEdge e = P;
		int cost = 1;
		while (e != Q.prevEdge) {
			cost++;
			e = e.prevEdge;
		}
		return cost;
	}
	public static int moveBehindCost(FundamentalEdge P, FundamentalEdge Q) {
		FundamentalEdge e = P.partner;
		int cost = 1;
		while (e != Q.partner.nextEdge) {
			cost++;
			e = e.nextEdge;
		}
		return cost;
	}
	
	public static void moveInFront(FundamentalEdge P, FundamentalEdge Q) {
		FundamentalPolygonUtility.bringTogetherViaB(Q.prevEdge, P);
	}
	public static boolean moveBehind(FundamentalEdge P, FundamentalEdge Q) {
		FundamentalPolygonUtility.bringTogetherViaA(P.partner, Q.partner.nextEdge);
		return true;
	}

	
	private static List<FundamentalEdge> findLargestCluster(FundamentalPolygon P) {
		Set<FundamentalEdge> cluster = new LinkedHashSet<FundamentalEdge>();
		Set<FundamentalEdge> actCluster = new LinkedHashSet<FundamentalEdge>();
		for (FundamentalEdge se : P.getEdges()) {
			actCluster.clear();
			FundamentalEdge actEdge = se;
			while (!actCluster.contains(actEdge.partner)) {
				actCluster.add(actEdge);
				actEdge = actEdge.partner.nextEdge;
			}
			if (cluster.size() < actCluster.size()) {
				cluster.clear();
				cluster.addAll(actCluster);
			}
		}
		return new LinkedList<FundamentalEdge>(cluster);
	}
	
	
	private static List<FundamentalEdge> findCoCluster(List<FundamentalEdge> cluster) {
		List<FundamentalEdge> result = new LinkedList<FundamentalEdge>();
		FundamentalEdge last = cluster.get(cluster.size() - 1);
		FundamentalEdge first = cluster.get(0);
		FundamentalEdge active = last.partner.nextEdge;
		while (active != first) {
			result.add(active);
			active = active.partner.nextEdge;
		}
		return result;
	}
		
	
}
