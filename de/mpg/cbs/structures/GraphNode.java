package de.mpg.cbs.structures;

/**
 *
 *  This class describes the nodes for the graph class
 *  for now, only handles simple graphs (bidirectional, no weights)
 *	
 *
 *	@version    October 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class GraphNode {
	
	//private float		value;
	private int			label;
	private int[]		connected;
	private int 		Nc;
	private	int			Nmax;
	
	private static final int MIN_NODES = 2;
	
	public GraphNode(int lb) {
		//value = 0.0f;
		label = lb;
		Nc = 0;
		Nmax = MIN_NODES;
		connected = new int[Nmax];
	}
	
	public GraphNode(int lb, int max) {
		//value = 0.0f;
		label = lb;
		Nc = 0;
		Nmax = max;
		connected = new int[max];
	}
	
	public GraphNode(int lb, int[] connect, int N) {
		//value = 0.0f;
		label = lb;
		connected = new int[N];
		for (int n=0;n<N;n++) connected[n] = connect[n];
		Nc = N;
		Nmax = N;
	}
		
	
	public void finalize() {
		connected = null;
		Nc = 0;
		Nmax = 0;
	}
	
	/**
	 *  set the value of the node
	 */
	//public final void setValue(float val) { value = val; }

	/**
	 *  add a link to the node
	 */
	public final void addLinkTo(int lb) {
		if (Nc>=Nmax) {
			int[] tmp = new int[Nc];
			for (int n=0;n<Nc;n++) tmp[n] = connected[n];
			Nmax = Math.max(2*Nmax,1);
			connected = new int[Nmax];
			for (int n=0;n<Nc;n++) connected[n] = tmp[n];
			tmp = null;
		}	
		connected[Nc] = lb;
		Nc++;
	}

	/**
	 *  remove a link from the node using its index
	 */
	public final void removeLinkAt(int indx) {
		for (int n=indx+1;n<Nc;n++) connected[n-1] = connected[n];
		Nc--;
	}

	/**
	 *  remove a link from the node using its label
	 */
	public final void removeLinkTo(int lb) {
		for (int n=0;n<Nc;n++) {
			if (connected[n]==lb) {
				removeLinkAt(n);
				break;
			}
		}
	}

	/**
	 * return the node label
	 */
	public final int getLabel() { return label; }
	public final int getLinkSize() { return Nc; }
	public final int getMaxSize() { return Nmax; }
	public final int[] getLinks() {
		return connected;
	}
	public final int getLinkAt(int n) {
		return connected[n];
	}
	public final void setLinkAt(int n, int lb) {
		connected[n] = lb;
	}
	public final boolean isLinkedTo(int lb) {
		for (int n=0;n<Nc;n++) {
			if (connected[n]==lb) return true;
		}
		return false;
	}
	public final boolean equals(GraphNode nod) {
		return (label==nod.getLabel());
	}
	public final boolean equals(int lb) {
		return (label==lb);
	}
	public final String printNode() {
		String message = new String();
		message = message.concat(""+label+" ->");
		for (int n=0;n<Nc;n++)
			message = message.concat(" "+connected[n]+",");
		return message;
	}
}
