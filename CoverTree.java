package com.sdsolve.tools;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public class CoverTree<E> {
	
	int maxLevel;
	int minLevel;
	double base;
	int maxMinLevel;
	Node<E> rootNode;
	boolean hasBounds;
	double[] min;
	double[] max;
	int[] numNodes;
	
	/**
	 * Create a cover tree at level zero which automatically expands above and below.
	 */
	public CoverTree () {
		this.maxMinLevel = Integer.MIN_VALUE;
		this.numNodes = new int[1000];
		this.base = 1.2;
	}
	
	/**
	 * Create a cover tree which stops increasing the minimumLevel as soon as the given number of nodes is reached.
	 */
	public CoverTree (double base, int maxMinLevel) {
		this.base = base;
		this.maxMinLevel = maxMinLevel;
		if (maxMinLevel>0) {
			this.maxLevel = maxMinLevel;
			this.minLevel = maxMinLevel;
		}
		this.numNodes = new int[1000];
	}
	
	/**
	 * Points outside of the bounding box, will not be included. This allows for easy truncation.
	 * @param min
	 * @param max
	 */
	public void setBounds(double[] min, double[] max) {
		hasBounds = true;
		this.min = min;
		this.max = max;
	}
	
	/**
	 * Returns the maximum level of this tree.
	 * @return
	 */
	public int maxLevel() {
		return maxLevel;
	}
	
	/**
	 * Returns the minimum level of this tree.
	 * @return
	 */
	public int minLevel() {
		return minLevel;
	}
	
	void incNodes(int level) {
		numNodes[level+500]++;
	}
	
	void decNodes(int level) {
		numNodes[level+500]--;
	}
	
	/**
	 * Returns the size of the cover tree up to the given level (inclusive)
	 * @param level
	 * @return
	 */
	public int size(int level) {
		int sum = 0;
		for (int i=maxLevel; i>=level; i--)
			sum += numNodes[i+500];
		return sum;
	}
	
	/**
	 * Returns the size of the cover tree
	 * @return
	 */
	public int size() {
		return size(minLevel);
	}
	
	void insertAtRoot(E element, double[] point) {
		//inserts the point above the root by successively increaasing the cover of the root node until it
		//contains the new point, the old root is added as child of the new root
		//inserts the point above the root node by increasing the maxlevel and replacing the old root with a node that contains the new point as root
		Node<E> oldRoot = rootNode;
		double dist = oldRoot.distance(point);
		while (dist > Math.pow(base,maxLevel)) {
			Node<E> newRoot = new Node<E>(null,rootNode.element,rootNode.point);
			rootNode.setParent(newRoot);
			newRoot.addChild(rootNode);
			rootNode = newRoot;
			decNodes(maxLevel);
			incNodes(++maxLevel);
			//System.out.println("insert "+rootNode.point+" as root at level "+(maxLevel));
		}
		Node<E> newNode = new Node<E>(rootNode,element,point);
		rootNode.addChild(newNode);
		incNodes(maxLevel-1);
		//System.out.println("insert "+point+" below root "+rootNode.point+" at level "+(maxLevel-1));
	}
	
	/**
	 * Insert a point into the tree. If the tree size is greater than k the lowest cover will be removed as long as it does not decrease tree size below k.
	 * @param point
	 */
	public boolean insert(E element, double[] point, int k) {
		boolean inserted = insert(element,point);
		//only do this if there are more than two levels
		if (maxLevel-minLevel>2) {
			//remove lowest cover if the cover before has a sufficient number of nodes
			if (size(minLevel+1)>=k) {
				removeLowestCover();
				//do not accept new nodes at the minimum level
				maxMinLevel = minLevel+1;
			}
			//remove redundant nodes from the minimum level
			if (size(minLevel)>=2*k) {
				removeNodes(k);
			}
		}
		return inserted;
	}
	
	/**
	 * Insert a point into the tree.
	 * @param point
	 */
	public boolean insert(E element, double[] point) {
		if (hasBounds) {
			//points outside of the bounding box will not be added to the tree
			for (int d=0; d<point.length; d++) {
				if (point[d]>max[d])
					return false;
				if (point[d]<min[d])
					return false;
			}
		}
		//if this is the first node make it the root node
		if (rootNode == null) {
			rootNode = new Node<E>(null,element,point);
			//List<Node> newCover = new LinkedList<Node>();
			//newCover.add(rootNode);
			//System.out.println("insert "+point+" as root at level "+maxLevel);
			incNodes(maxLevel);
			return true;
		}
		//do not add if the new node is identical to the root node
		rootNode.distance = rootNode.distance(point);
		if (rootNode.distance == 0.)
			return false;
		//if the node lies outside the cover of the root node and its decendants then insert the node above the root node
		if (rootNode.distance > Math.pow(base,maxLevel+1)) {
			insertAtRoot(element,point);
			return true;
		}
		//usually insertion begins here
		List<Node<E>> coverset = new LinkedList<Node<E>>();
		//the initial coverset contains only the root node
		coverset.add(rootNode);
		int level = maxLevel;
		Node<E> parent = null; //the root node does not have a parent
		int parentLevel = maxLevel;
		while (true) {
			boolean parentFound = true;
			List<Node<E>> candidates = new LinkedList<Node<E>>();
			for (Node<E> n1 : coverset) {
				for (Node<E> n2 : n1.getChildren()) {
					if (n1.point!=n2.point) {
						//do not compute distance twice
						n2.distance = n2.distance(point) ;
						//do not add if node is already contained in the tree
						if (n2.distance == 0.)
							return false;
					}
					else 
						n2.distance = n1.distance;
					if (n2.distance < Math.pow(base,level)) {
						candidates.add(n2);
						parentFound = false;
					}
				}
			}
			//if the children of the coverset are further away the 2^level then an element of the
			//coverset is the parent of the new node
			if (parentFound)
				break;
			//select one node of the coverset as the parent of the node
			for (Node<E> n : coverset) {
				//dist = n.distance(point);
				if (n.distance < Math.pow(base,level)) {
					parent = n;
					parentLevel = level;
					break;
				}
			}
			//set all nodes as the new coverset
			level--;
			coverset = candidates;
		}
		//if the point is a sibling of the root node, then the cover of the root node is increased
		if (parent == null) {
			insertAtRoot(element,point);
			return true;
		}
		if (parentLevel-1 < minLevel) {
			//if the maximum size is reached and this would only increase the depth of the tree then stop
			if (parentLevel-1 < maxMinLevel)
				return false;
			minLevel = parentLevel-1;
		}
		//otherwise add child to the tree
		//System.out.print("regular ");
		Node<E> newNode = new Node<E>(parent,element,point);
		parent.addChild(newNode);
		//record distance to parent node and add to the sorted set of nodes where distance is used for sorting (needed for removal)
		incNodes(parentLevel-1);
		//System.out.println("insert "+point+" below "+newNode.parent.point+" at level "+level);
//		for (int i=maxLevel; i>= minLevel; i--)
//			System.out.print(size(i)+" ");
//		System.out.println();
		return true;
	}
	
	/**
	 * Removes the the cover at the lowest level of the tree.
	 */
	void removeLowestCover() {
		List<Node<E>> coverset = new LinkedList<Node<E>>();
//		System.out.println("rootnode: "+rootNode.hashCode());
		coverset.add(rootNode);
//		System.out.println(maxLevel+" "+coverset.size());
		int k = maxLevel;
		while(k-- > minLevel+1){
			List<Node<E>> nextCoverset = new LinkedList<Node<E>>();
			for (Node<E> n : coverset) 
				nextCoverset.addAll(n.getChildren());
			coverset = nextCoverset;
//			System.out.println(k+" "+coverset.size());
		}
		for (Node<E> n : coverset) 
			n.removeChildren();
		
		minLevel++;
	}
	
	
	/**
	 * Removes all but k points.
	 */
	List<Node<E>> removeNodes(int numCenters) {
		List<Node<E>> coverset = new LinkedList<Node<E>>();
//		System.out.println("rootnode: "+rootNode.hashCode());
		coverset.add(rootNode);
//		System.out.println(size(minLevel+1)+" "+size(minLevel));
		int k = maxLevel;
		while(k-- > minLevel+1){
			List<Node<E>> nextCoverset = new LinkedList<Node<E>>();
			for (Node<E> n : coverset) 
				nextCoverset.addAll(n.getChildren());
			coverset = nextCoverset;
//			System.out.println(k+" "+coverset.size());
		}
		int missing = numCenters-coverset.size();
		if (missing < 0)
			System.err.println("Error: negative missing="+missing+" in coverset");
		//sucessively pick the node with the largest distance to the coverset and add it to the coverset
		LinkedList<Node<E>> candidates = new LinkedList<Node<E>>();
		for (Node<E> n : coverset) 
			for (Node<E> n2 : n.getChildren())
				if (n.point!=n2.point)
					candidates.add(n2);
		//only add candidates when the coverset is yet smaller then the number of desired centers
		if (coverset.size()<numCenters) {
			//compute the distance of all candidates to their parents and uncles
			for (Node<E> n1 : candidates) {
				double minDist = Double.POSITIVE_INFINITY;
				for (Node<E> n2 : n1.getParent().getParent().getChildren()) {
					double dist = n1.distance(n2.point);
					if (dist < minDist)
						minDist = dist;
				}
				n1.distance = minDist;
				if (minDist==Double.POSITIVE_INFINITY)
					System.err.println("Error: Infinite distance in k centers computation.");
			}
			do {
				Collections.sort(candidates);
				Node<E> newNode = candidates.removeLast();
				coverset.add(newNode);
				//update the distance of all candidates in the neighborhood of the new node
				for (Node<E> n : newNode.getParent().getParent().getChildren()) {
					if (n!=newNode) {
						double dist = newNode.distance(n.point);
						if (dist < newNode.distance)
							newNode.distance = dist;
					}
				}
			} while (coverset.size()<numCenters);
		}
		//finally remove all nodes that have not been selected from the tree to avoid confusing the nearest neighbor computation
		for (Node<E> n : candidates) {
			n.getParent().removeChild(n);
			decNodes(minLevel);
		}
		return coverset;
	}
	
	public E getNearest(double[] point) {
		List<Node<E>> candidates = new LinkedList<Node<E>>();
		candidates.add(rootNode);
		rootNode.distance = rootNode.distance(point);
		double minDist = rootNode.distance;
		for (int i=maxLevel; i>minLevel; i--) {
			List<Node<E>> newCandidates = new LinkedList<Node<E>>();
			for (Node<E> n : candidates) {
				for (Node<E> n2 : n.getChildren()) {
					if (n.point!=n2.point) {
						n2.distance = n2.distance(point);
						if (n2.distance<minDist)
							minDist = n2.distance;
					}
					else
						n2.distance = n.distance;
					newCandidates.add(n2);
				}
			}
			candidates.clear();
			for (Node<E> n : newCandidates)
				if (n.distance < minDist + Math.pow(base,i))
					candidates.add(n);
		}
		//Node nearest = null;
		//minDist = rootNode.distance;
		for (Node<E> n : candidates) {
			if (n.distance == minDist) {
				return n.element;
				//minDist = n.distance;
				//nearest = n;
			}
		}
		return null;//nearest.point;
	}
	
	/**
	 * Get the cover of the given level. All points at this level are guaranteed to be 2^i apart from one another.
	 * @param level
	 * @return
	 */
	public List<E> getCover(int level) {
//		System.out.println(level);
		List<Node<E>> coverset = new LinkedList<Node<E>>();
//		System.out.println("rootnode: "+rootNode.hashCode());
		coverset.add(rootNode);
//		System.out.println(maxLevel+" "+coverset.size());
		int k = maxLevel;
//		System.out.println(maxLevel+" "+level+" "+minLevel);
		while(k-- > level){
			List<Node<E>> nextCoverset = new LinkedList<Node<E>>();
			for (Node<E> n : coverset) 
				nextCoverset.addAll(n.getChildren());
			coverset = nextCoverset;
//			System.out.println(k+" "+coverset.size());
		}
//		System.out.println(level+" has size "+coverset.size());
		List<E> cover = new LinkedList<E>();
//		System.out.println("get cover: "+coverset.size()+" "+size(level));
		for (Node<E> n: coverset) {
			cover.add(n.element);
//			System.out.println("get cover: "+n.point.hashCode()/1000000);
		}
		
		return cover;
	}
	
	/**
	 * Gets at least k centers which are maximally apart from each other. All remaining centers are removed from the tree. This function only works as designed
	 * when the function insert(point,k) has been used before to add points to the tree. Otherwise, it will return the cover one level above the bottom most level of the tree.
	 * @param number of centers
	 * @return
	 */
	public List<E> getKCenters(int numCenters) {
		List<Node<E>> coverset = removeNodes(numCenters);
		//create cover
		List<E> cover = new LinkedList<E>();
		for (Node<E> n: coverset) {
			cover.add(n.element);
//			System.out.println(n.point[0]+"\t"+n.point[1]);
		}
		return cover;
		
	}
	
	static double distance(double[] d1, double[] d2) {
		double sumSq = 0.;
		for (int i=0; i<d1.length; i++) {
			double d = d1[i]-d2[i];
			sumSq += d*d;
		}
		return sumSq;
	}
	
	@SuppressWarnings("hiding")
	class Node<E> implements Comparable<Node<E>> {
		
		Node<E> parent; 
		E element;
		List<Node<E>> children;
		double[] point;
		double distance;
		
		//use for a child
		Node (Node<E> parent, E element, double[] point) {
			this.parent = parent;
			this.children = new LinkedList<Node<E>>();
			this.element = element;
			this.point = point;
//			if (parent!=null)
//				System.out.println("insert "+Arrays.toString(point)+" ("+this.hashCode()/1000000+") as child of "+Arrays.toString(parent.point)+" ("+parent.hashCode()/1000000+") at level "+level);
//			else
//				System.out.println("insert "+Arrays.toString(point)+" ("+this.hashCode()/1000000+") as child of null at level "+level);
		}
		
		Node<E> getParent() {
			return parent;
		}
		
		void setParent(Node<E> node) {
			parent = node;
		}
		
		void addChild(Node<E> node) {
			children.add(node);
		}
		
		List<Node<E>> getChildren() {
			if (children.isEmpty()) {
				Node<E> n = new Node<E>(this,this.element,this.point);
				addChild(n);
			}
			return children;
		}
		
		void removeChild(Node<E> n) {
			children.remove(n);
		}
		
		void removeChildren() {
			children.clear();
		}
		
		double distance(double[] point) {
			double sumSq = 0.;
			for (int i=0; i<point.length; i++) {
				double d = this.point[i]-point[i];
				sumSq += d*d;
			}
			sumSq = Math.sqrt(sumSq);
			return sumSq;
		}

		@Override
		public int compareTo(Node<E> o) {
			if (distance < o.distance)
				return -1;
			if (distance > o.distance)
				return 1;
			return 0;
		}
		
		
	}

}
