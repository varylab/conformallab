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
		List<FundamentalEdge> sorted = findLongestSortedSubCluster(C, coC);
		System.out.println("Sorted SubCluster ######### " + sorted);
		
		List<FundamentalEdge> toSort = new LinkedList<FundamentalEdge>(C);
		toSort.removeAll(sorted);
		
		int numIterations = toSort.size();
		
		for (int i = 0; i < numIterations; i++) {
			sorted = findLongestSortedSubCluster(C, coC);
			toSort = new LinkedList<FundamentalEdge>(C);
			toSort.removeAll(sorted);
			FundamentalEdge e = toSort.get(0);
			System.out.println("Inserting " + e + " or " + e.partner + " into sorted cluster");
			insertIntoSortedCluster(e, sorted, C, coC);
			System.out.println("IN SORT: Largest Cluster #####: " + C);
			System.out.println("IN SORT: Co-Cluster #####: " + coC);
		}
	}
	
	
	public static void insertIntoSortedCluster(FundamentalEdge e, List<FundamentalEdge> sorted, List<FundamentalEdge> C, List<FundamentalEdge> coC) {
		List<FundamentalEdge> sortedPartners = getPartnerEdges(sorted);
		int sIndex = -1;
		int sCoIndex = -1;
		int tIndex = C.size();
		int tCoIndex = coC.size();
		
		
		int ePartnerIndex = coC.indexOf(e.partner);
		for (int i = ePartnerIndex + 1; i < coC.size(); i++) {
			FundamentalEdge check = coC.get(i);
			if (sortedPartners.contains(check)) {
				tIndex = C.indexOf(check.partner);
				break;
			}
		}
		for (int i = ePartnerIndex - 1; i >= 0; i--) {
			FundamentalEdge check = coC.get(i);
			if (sortedPartners.contains(check)) {
				sIndex = C.indexOf(check.partner);
				break;
			}
		}
		
		int eIndex = C.indexOf(e);
		for (int i = eIndex + 1; i < C.size(); i++) {
			FundamentalEdge check = C.get(i);
			if (sorted.contains(check)) {
				tCoIndex = coC.indexOf(check.partner);
				break;
			}
		}
		for (int i = eIndex - 1; i >= 0; i--) {
			FundamentalEdge check = C.get(i);
			if (sorted.contains(check)) {
				sCoIndex = coC.indexOf(check.partner);
				break;
			}
		}
		System.out.println(sIndex + ", " + tIndex + " - " + sCoIndex + ", " + tCoIndex);
		
		
		FundamentalEdge frontEdge = null;
		int cost = Integer.MAX_VALUE;
		FundamentalEdge frontCoEdge = null;
		int coCost = Integer.MAX_VALUE;
		for (int i = sIndex; i < tIndex; i++) {
			FundamentalEdge actEdge = i != -1 ? C.get(i) : coC.get(coC.size() - 1);
			if (canMoveBehind(e, actEdge)) {
				int c = moveBehindCost(e, actEdge);
				if (cost > c) {
					cost = c;
					frontEdge = actEdge;
					break;
				}
			}
		}
		for (int i = sCoIndex; i < tCoIndex; i++) {
			FundamentalEdge actEdge = i != -1 ? coC.get(i) : C.get(C.size() - 1);
			if (canMoveBehind(e.partner, actEdge)) {
				int c = moveBehindCost(e.partner, actEdge);
				if (coCost > c) {
					coCost = c;
					frontCoEdge = actEdge;
					break;
				}
			}
		}
		System.out.println("could move edge " + e + " behind " + frontEdge + " cost " + cost);
		System.out.println("could move edge " + e.partner + " behind " + frontCoEdge + " cost " + coCost);
		assert (cost != Integer.MAX_VALUE || coCost != Integer.MAX_VALUE);
		if (cost < coCost) {
			moveBehind(e, frontEdge);
			C.remove(e);
			int i = C.indexOf(frontEdge);
			C.add(i + 1, e);
		} else {
			moveBehind(e.partner, frontCoEdge);
			coC.remove(e.partner);
			int i = coC.indexOf(frontCoEdge);
			coC.add(i + 1, e.partner);
		}
		// TODO check if moveToFront is cheaper
	}
	
	
	
	public static List<FundamentalEdge> findLongestSortedSubCluster(List<FundamentalEdge> C, List<FundamentalEdge> coC) {
		int l = 0;
		List<FundamentalEdge> R = null;
		for (FundamentalEdge e : C) {
			List<FundamentalEdge> sC = getSortedSubClusterAt(e, C, coC);
			if (l < sC.size()) {
				R = sC;
				l = R.size();
			}
		}
		return R;
	}
	
	
	public static List<FundamentalEdge> getSortedSubClusterAt(FundamentalEdge start, List<FundamentalEdge> C, List<FundamentalEdge> coC) {
		List<FundamentalEdge> R = new LinkedList<FundamentalEdge>();
		R.add(start);
		int coIndex = coC.indexOf(start.partner);
		int i = C.indexOf(start) + 1;
		List<FundamentalEdge> recSubCluster = new LinkedList<FundamentalEdge>();
		for (;i < C.size(); i++) {
			FundamentalEdge next = C.get(i);
			int coNextIndex = coC.indexOf(next.partner); 
			if (coIndex < coNextIndex) {
				List<FundamentalEdge> recCluster = getSortedSubClusterAt(next, C, coC);
				if (recCluster.size() > recSubCluster.size()) {
					recSubCluster = recCluster;
				}
			}
		}
		R.addAll(recSubCluster);
		return R;
	}
	
	
	public static void moveIntoSortedSubCluster(FundamentalEdge e, List<FundamentalEdge> C) {
		
	}
	

	/**
	 * moves e or e.partner into the given cluster
	 * @param e
	 * @param cluster
	 */
	public static void moveIntoCluster(FundamentalEdge e, List<FundamentalEdge> cluster) {
		FundamentalEdge pre = cluster.get(0).prevEdge.partner;
		FundamentalEdge post = cluster.get(cluster.size() - 1).partner.nextEdge;
		List<FundamentalEdge> inFrontEdges = new LinkedList<FundamentalEdge>();
		List<FundamentalEdge> behindEdges = new LinkedList<FundamentalEdge>();
		inFrontEdges.addAll(cluster);
		inFrontEdges.add(post);
		behindEdges.add(pre);
		behindEdges.addAll(cluster);
		FundamentalEdge P = e;
		FundamentalEdge Q = null;
		int minCost = Integer.MAX_VALUE;
		boolean moveInFront = true;
		boolean movePartner = false;
		for (FundamentalEdge active : inFrontEdges) {
			if (canMoveInFront(P, active)) {
				int cost = moveInFrontCost(P, active); 
				if (minCost > cost) {
					Q = active;
					minCost = cost;
					moveInFront = true;
					movePartner = false;
				}
			}
			if (canMoveInFront(P.partner, active)) {
				int cost = moveInFrontCost(P.partner, active); 
				if (minCost > cost) {
					Q = active;
					minCost = cost;
					moveInFront = true;
					movePartner = true;
				}
			}
		}
		for (FundamentalEdge active : behindEdges) {
			if (canMoveBehind(P, active)) {
				int cost = moveBehindCost(P, active);
				if (minCost > cost) {
					Q = active;
					minCost = cost;
					moveInFront = false;
					movePartner = false;
				}
			}
			if (canMoveBehind(P.partner, active)) {
				int cost = moveBehindCost(P.partner, active);
				if (minCost > cost) {
					Q = active;
					minCost = cost;
					moveInFront = false;
					movePartner = true;
				}
			}			
		}
		assert Q != null;
		if (Q == null) {
			throw new RuntimeException("Cannot move generator into cluster.");
		}
		FundamentalEdge Pmove = movePartner ? P.partner : P;
		if (moveInFront) {
			moveInFront(Pmove, Q);
			int index = cluster.indexOf(Q);
			cluster.add(index, Pmove);
		} else {
			moveBehind(Pmove, Q);
			int index = cluster.indexOf(Q);
			cluster.add(index + 1, Pmove);
		}
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
		System.out.println("move in front " + P + " <-> " + Q);
		FundamentalPolygonUtility.bringTogetherViaB(Q.prevEdge, P);
	}
	public static boolean moveBehind(FundamentalEdge P, FundamentalEdge Q) {
		System.out.println("move behind " + P + " <-> " + Q);
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
		
	
	public static List<FundamentalEdge> getPartnerEdges(List<FundamentalEdge> edges) {
		List<FundamentalEdge> R = new LinkedList<FundamentalEdge>();
		for (FundamentalEdge e : edges) {
			R.add(e.partner);
		}
		return R;
	}
	
}
