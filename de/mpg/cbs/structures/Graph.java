package de.mpg.cbs.structures;

import java.io.*;
import java.util.*;

/**
 *
 *  This class handles general graphs, and perform simple operations
 *  (loop detection and breaking, etc)
 *	
 *
 *	@version    October 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class Graph {
	
	private Vector			nodes;
	private int				Nn;
	
	private int 		Nloops;
	private int 		Nparts;
	
	public static final	int		NOTFOUND = -1;
	
	public Graph() {
		nodes = new Vector();
		Nn = 0;
	}
	
	public Graph(GraphNode[] nod, int N) {
		nodes = new Vector();
		for (int n=0;n<N;n++) nodes.add(nod[n]);
		Nn = N;
	}
	
	public void finalize() {
		nodes.clear();
		nodes = null;
		Nn = 0;
	}
	
	/**
	 *  add a new node into the graph
	 */
	public final void addNode(GraphNode nod) {
		nodes.add(nod);
		Nn++;
	}//addNode
	
	/**
	 *  replace a node in the graph
	 */
	public final void setNodeAt(int indx, GraphNode nod) {
		nodes.set(indx, nod);
	}//addNode
	
	/**
	 *  add a new node into the graph
	 */
	public final void addNode(int label) {
		GraphNode nod = new GraphNode(label);
		nodes.add(nod);
		Nn++;
	}//addNode
	
	/**
	 *  remove a node from the graph
	 */
	public final void removeNodeAt(int indx) {
		nodes.remove(indx);
		Nn--;
	}// removeNode

	/**
	 *  remove a node from the graph
	 */
	public final void removeNode(int label) {
		for (int n=0;n<Nn;n++) {
			if (((GraphNode)nodes.get(n)).equals(label)) {
				nodes.remove(n);
				break;
			}
		}
		Nn--;
	}// removeNode

	/**
	 * return the first value and its coordinates
	 */
	public final GraphNode getNodeAt(int indx) {
		return (GraphNode)nodes.get(indx);
	}
	public final GraphNode getNode(int label) {
		for (int n=0;n<Nn;n++) {
			if (((GraphNode)nodes.get(n)).equals(label)) {
				return (GraphNode)nodes.get(n);
			}
		}
		return null;
	}
	public final int getNodeIndex(int label) {
		for (int n=0;n<Nn;n++) {
			if (((GraphNode)nodes.get(n)).equals(label)) {
				return n;
			}
		}
		return NOTFOUND;
	}
	public final boolean containsNode(int label) {
		for (int n=0;n<Nn;n++) {
			if (((GraphNode)nodes.get(n)).equals(label)) {
				return true;
			}
		}
		return false;
	}
	public final void addNodeLink(int indx1, int indx2) {
		GraphNode nod1 = getNodeAt(indx1);
		GraphNode nod2 = getNodeAt(indx2);
		nod1.addLinkTo(nod2.getLabel());
		nod2.addLinkTo(nod1.getLabel());
		setNodeAt(indx1,nod1);
		setNodeAt(indx2,nod2);
		
		return;
	}
	public final int getSize() { return Nn; }
	
	public final int getLoops() { return Nloops; }
	
	public final int getParts() { return Nparts; }
	
	public final String printGraph() {
		String message = new String("graph("+Nn+":\n");
		for (int n=0;n<Nn;n++) {
			message = message.concat("node "+n+": ");
			message = message.concat(getNodeAt(n).printNode()+"\n");
		}
		return message;
	}
	
	/**
	 *	count the number of sub-graphs in a graph
	 */
	public final int countSubGraphs() { 
		int                     n;
		boolean[]				visited = new boolean[Nn];
	
		for (n=0;n<Nn;n++) {
			visited[n] = false;
		}
    
		Nparts = 0;
		for (n=0;n<Nn;n++) {
			if (!visited[n]) {
				// new graph piece: visit it completely
				Nparts++;
				visitGraph(n, visited);
			}
		} 
		return Nparts;
	}
	
	/**
	 * recursive visit of every node in the graph linked to the current node
	 *  when loops are present, they are not taken into account
	 */
	private final void visitGraph(int indx, boolean[] visited) {
		int      n,next;
		
		// the node is now visited
		visited[indx] = true;

		// linked ?
		for (n=0;n<getNodeAt(indx).getLinkSize();n++) {
			next = getNodeIndex(getNodeAt(indx).getLinkAt(n));
			if (!visited[next]) {
				visitGraph(next, visited);
			}
		}
		return;
	}

	/**
	 *	count the number of separate sub-graph and the number of loops.
	 *	<p>
	 *  destroys all connections!
	 *  it is necessary to duplicate the graph first if you need to keep it.
	 */
	public final int countLoopsAndSubGraphs() { 
		int                     n,k;
		boolean[]				visited = new boolean[Nn];
	
		for (n=0;n<Nn;n++) {
			visited[n] = false;
		}

		Nloops = 0;
		Nparts = 0;

		for (n=0;n<Nn;n++) {
			if (!visited[n]) {
				// visit a complete graph each time
				Nparts++;
				loopedPath(n, NOTFOUND, visited);
			}
		}

		return Nloops;
	}

	/**
	 *	follow a path into the graph and breaks all connections
	 *	(allows to count the number of loops)
	 */
	void loopedPath(int curr, int prec, boolean[] visited) {
		int      n,k;
		GraphNode nod = getNodeAt(curr);
		GraphNode linked;
		int			next;
		//int		indx = getNodeIndex(curr);
		
		// the node is now visited
		visited[curr] = true;
    
		// check for the link
		for (n=0;n<nod.getLinkSize();n++) {
			next = getNodeIndex(nod.getLinkAt(n));
			// avoids precedent node
			if ( next != prec) {
				if (!visited[next]) {
					loopedPath(next, curr, visited);
				} else {
					// there is a loop!
					Nloops++;
					// the loop is destroyed by suppressing the link, and continuing..
					linked = getNodeAt(next);
					linked.removeLinkTo(curr);
					setNodeAt(next, linked);
					nod.removeLinkAt(n);
					setNodeAt(curr, nod);
				}
            }
        }
		return;
	}// loopedPath

}
