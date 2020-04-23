package framework;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.configuration.BaseConfiguration;
import org.apache.jena.graph.Node;
import org.apache.jena.graph.Triple;
import org.apache.jena.riot.RDFParser;
import org.apache.jena.riot.lang.PipedRDFIterator;
import org.apache.jena.riot.lang.PipedRDFStream;
import org.apache.jena.riot.lang.PipedTriplesStream;
import org.apache.tinkerpop.gremlin.process.traversal.dsl.graph.GraphTraversalSource;
import org.apache.tinkerpop.gremlin.structure.Direction;
import org.apache.tinkerpop.gremlin.structure.Edge;
import org.apache.tinkerpop.gremlin.structure.Graph;
import org.apache.tinkerpop.gremlin.structure.T;
import org.apache.tinkerpop.gremlin.structure.Vertex;
import org.apache.tinkerpop.gremlin.structure.util.GraphFactory;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.math3.util.Pair;

import com.arangodb.ArangoCollection;
import com.arangodb.ArangoDB;
import com.arangodb.ArangoDatabase;
import com.arangodb.tinkerpop.gremlin.utils.ArangoDBConfigurationBuilder;


interface SynthToolkit {
	String getLabel();
	String randomLabel();
	boolean randomDirectedEdge(Vertex vx1, Vertex vx2, String label);
	boolean attemptAddEdge(Vertex vx1, Vertex vx2, String label);
}


public class ArangoBridge {
	private String name;
	private Graph graph;
	private GraphTraversalSource gts; 	//original traverser
	private GraphTraversalSource g; 	//clone traverser
	private int tid;
	private ArangoCollection edgeCol;
	private ArangoCollection nodeCol;
	
	private final static String HOST_IP = "127.0.0.1";
	private final static int HOST_PORT = 8529;
	private final static String DB_NAME = "";		// TODO: Fill in according to your ArangoDB configuration
	private final static String DB_USER = "root";
	private final static String DB_PASS = "";		// TODO: Fill in according to your ArangoDB configuration
	
	
	/**
	 * Constructor for creating/loading graphs.
	 * @param name - Name of graph to be created/loaded.
	 * @param parse - If a dataset should be parsed from RDF, this should be true.
	 */
	public ArangoBridge(String name, boolean parse) {
		ArangoDBConfigurationBuilder builder = new ArangoDBConfigurationBuilder();
		
		builder.graph(name)
			.withVertexCollection("node")
			.withEdgeCollection("edge");
	
		builder.arangoHosts(HOST_IP + ":" + HOST_PORT)
			.dataBase(DB_NAME)
			.arangoUser(DB_USER)
			.arangoPassword(DB_PASS);
		
		ArangoDB arango = new ArangoDB.Builder()
				.host(HOST_IP, HOST_PORT)
				.user(DB_USER)
				.password(DB_PASS)
				.build();
		
		ArangoDatabase adb = arango.db(DB_NAME);
		edgeCol = adb.collection(name + "_edge");
		nodeCol = adb.collection(name + "_node");
		
		BaseConfiguration conf = builder.build();
		this.name = name;
		graph = GraphFactory.open(conf);
		gts = new GraphTraversalSource(graph);
		cloneGTS();
		
		long node_cnt = (nodeCol.count().getCount());
		tid = (int) ((node_cnt*2)+1);
		
		if (parse) {
			if (tid == 1) { 
				System.out.println("PARSING!");
				parseRDF(name); 
			}
		}
	}
	

	public Graph getGraph() { return graph; }
	public int getTid() { return tid; }
	public String getName() { return name; }
	
	
	public ArrayList<Vertex> getNodes() {
		ArrayList<Vertex> nodes = new ArrayList<Vertex>();
		Iterator<Vertex> it = graph.vertices();
		while (it.hasNext()) {
			nodes.add(it.next());
		}
		return nodes;
	}
	
	
	public ArrayList<Edge> getEdges() {
		ArrayList<Edge> edges = new ArrayList<Edge>();
		Iterator<Edge> it = graph.edges();
		while (it.hasNext()) {
			edges.add(it.next());
		}
		return edges;
	}
	
	
	/**
	 * Get the cloned version of the GraphTraversalSource.
	 * @return cloned GraphTraversalSource.
	 */
	public GraphTraversalSource getG() {
		return g;
	}
	
	
	public void cloneGTS() {
		g = gts.clone();
	}
	
	
	/**
	 * Main function of the synthesizing graph process. The process 
	 * is divided into two sub-functions, so that the edge synthesizing 
	 * process can be run multiple times during runtime.
	 * @param nodes - Number of nodes.
	 * @param edges - Number of edges.
	 * @param nol - Standard deviation of normal distribution.
	 * @param mcs - Minimum community size.
	 */
	public void synthesize(int nodes, int edges, int nol, int mcs, ArrayList<HashMap<String,Integer>> lc) {	
		if (!verifySynth(nodes, edges, nol, mcs)) { return; }
		
		synthesizeNodes(nodes);
		synthesizeEdges(nodes, edges, nol, mcs, lc);
		
		System.out.println("Nodes: " + nodeCol.count().getCount() + "/" + nodes);
		System.out.println("Edges: " + edgeCol.count().getCount() + "/" + edges);
	}
	
	
	/**
	 * Help function to determine the amount of possible connections between a number of 
	 * nodes. Considers edge direction, but not labels.
	 * @param nodes - Number of nodes to consider.
	 * @return
	 */
	public long detCon(int nodes) {
		long n = 2; 
		long c = 1;
		while (n < nodes) {
			c += n;
			n++;
		}
		return 2*c;
	}
	
	
	/**
	 * Checks that parameters are "legal".
	 * A threshold of 0.75 is in place to avoid lockup when generating edges. 
	 * For values closer to 1, the chance of randomizing an edge that doesn't exist 
	 * becomes too low. The synthesizer runs slower for higher threshold values.
	 * Threshold values were tuned with a low number of nodes and very high number 
	 * of edges.
	 * @param nodes - Number of nodes.
	 * @param edges - Number of edges.
	 * @param nol - Standard deviation of normal distribution.
	 * @param mcs - Minimum community size.
	 * @return result of verification.
	 */
	public boolean verifySynth(int nodes, int edges, int nol, int mcs) {
		long con = detCon(nodes);
		
		if (nodes > edges) {
			System.out.println("Invalid graph. There are more nodes than edges: " + nodes + " > " + edges);
			return false;
		} else if (edges - (Math.floor((edges/nodes) - 1.0)) * nodes < Math.floor(nodes/mcs)) {
			System.out.println("Invalid graph. Communities can't be interconnected: " + (edges - (Math.floor((edges/nodes) - 1.0)) * nodes + " < " + Math.floor(nodes/mcs)));
			return false;
		} else if (0.75 * nol * con < edges) {
			System.out.println("Invalid graph. Not enough labels for that many edges: " + (int) (0.75 * nol * con) + " allowed relations vs. " + edges + " edges");
			return false;
		}
		return true;
	}
	
	
	/**
	 * Creates all nodes in the synthetic graph.
	 * @param nodes - Number of nodes.
	 */
	public void synthesizeNodes(int nodes) {
		int nodes_cnt = nodes;
		ArrayList<Vertex> all_nodes = new ArrayList<Vertex>();
		
		// Adding all nodes
		while (nodes_cnt > 0) {
			all_nodes.add(g.addV("node").property(T.id, Integer.toString(tid))
				.property("name", "node_" + tid).next());
			cloneGTS(); // ?
			tid += 2;
			nodes_cnt--;
			System.out.println("Node added: " + (nodes-nodes_cnt) + "/" + nodes);
		}
	}
	
	
	/**
	 * Synthesizes the graph edges in three different stages.
	 * Stage 1: form communities and add edges to ensure connectivity.
	 * Stage 2: interconnect communities to ensure connectivity.
	 * Stage 3: randomly apply rest of the edges.
	 * @param nodes - Number of nodes.
	 * @param edges - Number of edges.
	 * @param nol - Standard deviation of normal distribution.
	 * @param mcs - Minimum community size.
	 */
	public void synthesizeEdges(int nodes, int edges, int nol, int mcs, ArrayList<HashMap<String,Integer>> lc) {
		double edges_cnt = edges;
		int cmty_cnt = 0;
		int cmty_rest_cnt = 0; 
		int cmty_intr_cnt = 0; 
		int rnd_cnt = 0;
		int out_idx, in_idx;
		String label;
		Random rnd = new Random();
		ArrayList<Vertex> all_nodes = getNodes();
		ArrayList<String> all_labels = new ArrayList<String>();
		ArrayList<String> waste_labels = new ArrayList<String>();
		LinkedHashSet<ImmutableTriple<Vertex,Vertex,String>> edge_list = new LinkedHashSet<ImmutableTriple<Vertex,Vertex,String>>();
		HashMap<Pair<Vertex,String>,Integer> lc_map = new HashMap<Pair<Vertex,String>,Integer>();
		Random ld = new Random();
		
		for (String s : lc.get(0).keySet()) {
			all_labels.add(s);
		}
		
		// Initiating the method of keeping track of node & constraints
		for (Vertex v : all_nodes) {
			for (String s : lc.get(0).keySet()) {
				lc_map.put(new Pair<Vertex,String>(v,s), 0);
			}
		}
		
		
		/**
		 * Toolkit with help-functions for graph synthesizing.
		 * getLabel() - Retrieves a label from waste pile if there are any, 
		 * 				otherwise calls function to generate a new one.
		 * randomLabel() - Generate random label based on label distribution.
		 * randomDirectedEdge(...) - Randomizes edge direction and calls 
		 * 							 function to attempt to add edge.
		 * attemptAddEdge(...) - Attempts to add an edge between two nodes.
		 * 						 returns false if such edge already exist, true otherwise. 
		 */
		SynthToolkit st = new SynthToolkit() {
			@Override
			public String getLabel() {
				if (waste_labels.size() > 0) {
					return waste_labels.remove(rnd.nextInt(waste_labels.size()));
				} else {
					return randomLabel();
				}
			}
			
			
			@Override
			public String randomLabel() {
				return all_labels.get(ld.nextInt(nol));
			}

			
			@Override
			public boolean randomDirectedEdge(Vertex vx1, Vertex vx2, String label) {
				Random rnd = new Random();
				Vertex start, end;
				int coin = rnd.nextInt(2);
				if (coin == 1) {
					start = vx1;
					end = vx2;
				} else {
					end = vx1;
					start = vx2;
				}
				return attemptAddEdge(start, end, label);
			}
			
			
			@Override
			public boolean attemptAddEdge(Vertex vx1, Vertex vx2, String label) {
				Integer cval = lc_map.get(new Pair<Vertex,String>(vx1,label));
				Integer ub = lc.get(1).get(label);
				int curr;
				
				if (ub == null || cval < ub) {
					if (!edge_list.contains(new ImmutableTriple<Vertex, Vertex, String>(vx1, vx2, label))) {
						g.addE("edge").from(vx1).to(vx2).property("label", label).next();
						cloneGTS();
						edge_list.add(new ImmutableTriple<Vertex, Vertex, String>(vx1, vx2, label));
						curr = lc_map.get(new Pair<Vertex,String>(vx1,label));
						lc_map.put(new Pair<Vertex,String>(vx1,label), curr + 1);
						return true;
					} 
				}
				cval = lc_map.get(new Pair<Vertex,String>(vx2,label));
				
				if (ub == null || cval < ub) {
					if (!edge_list.contains(new ImmutableTriple<Vertex, Vertex, String>(vx2, vx1, label))) {
						g.addE("edge").from(vx2).to(vx1).property("label", label).next();
						cloneGTS();
						edge_list.add(new ImmutableTriple<Vertex, Vertex, String>(vx2, vx1, label));
						curr = lc_map.get(new Pair<Vertex,String>(vx2,label));
						lc_map.put(new Pair<Vertex,String>(vx2,label), curr + 1);
						return true;
					}
				}
				return false;
			}
		};
		
		
		Vertex node_1, node_2;
		int max_cmty_cnt = (int) Math.floor((double) nodes / (double) mcs);
		int cmty_limit = (int) Math.floor((edges_cnt-(nodes-max_cmty_cnt)-(max_cmty_cnt-1)) / (double) max_cmty_cnt); // edges available  for randomization
		int cmty_size = 0;
		ArrayList<Vertex> cmty = new ArrayList<Vertex>();
		ArrayList<Vertex> cmty_refs = new ArrayList<Vertex>();
		int idx = 0;
		int all_nodes_size = all_nodes.size();
		
		
		// Construct communities
		while (idx < all_nodes_size) {
			cmty_size = mcs + rnd.nextInt(6);
			
			if (idx + cmty_size + mcs - 1 >= all_nodes.size()) { 
				cmty_size += (all_nodes.size() - (idx + cmty_size));
			}
			
			// Fill community with nodes
			cmty_refs.add(all_nodes.get(idx));
			for (int i = 0; i < cmty_size; i++) {
				cmty.add(all_nodes.get(idx+i));
			}
			
			// Add edges in community to ensure connectivity
			for (int i = cmty_size-2; i >= 0; i--) {
				label = st.getLabel();
				while(true) {
					node_1 = cmty.get(rnd.nextInt(i+1));
					node_2 = cmty.get(i+1);
					if (st.randomDirectedEdge(node_1, node_2, label)) {	
						edges_cnt--;
						cmty_cnt++;
						System.out.println("CC-edge added: " + (int)(edges-edges_cnt) + "/" + edges);
						break;
					} else {
						waste_labels.add(label);
						label = st.randomLabel();
					}
				}
			}		
			
			long pc = detCon(cmty_size);							// possible connections in community
			int ccs = cmty_size-1;									// current community size
			int cre = 0;											// community random edges
			double rrp = (pc*nol - ccs)/(pc*(double)nol);			// remaining relations percentage
			double threshold = 0.25;								// limit for the rrp value, once passed then community is done
			
			// Add random edges to community
			while (rrp > threshold && cre < cmty_limit) {
				out_idx = rnd.nextInt(cmty.size());
				label = st.getLabel();
				while (true) {
					in_idx = rnd.nextInt(cmty.size());
					if (in_idx != out_idx) { 
						node_1 = cmty.get(out_idx);
						node_2 = cmty.get(in_idx);
						if (st.randomDirectedEdge(node_1, node_2, label)) {
							edges_cnt--;
							cmty_rest_cnt++;
							ccs++;
							cre++;
							rrp = (pc*nol - ccs)/(pc*(double)nol);
							System.out.println("CR-edge added: " + (int)(edges-edges_cnt) + "/" + edges);
							break;
						} else {
							waste_labels.add(label);
							label = st.randomLabel();
						}
					}
				}
			}
			idx += cmty.size();
			cmty.clear();
		}
		
		// Interconnect communities to achieve connectivity
		for (int i = 1; i < cmty_refs.size(); i++) {
			label = st.getLabel();
			while(true) {
				node_1 = cmty_refs.get(rnd.nextInt(i));
				node_2 = cmty_refs.get(i);
				if (st.randomDirectedEdge(node_1, node_2, label)) {
					edges_cnt--;
					cmty_intr_cnt++;
					System.out.println("Interconnecting edge added: " + (int)(edges-edges_cnt) + "/" + edges);
					break;
				} else {
					waste_labels.add(label);
					label = st.randomLabel();
				}
			}
		}
		
		// Add remaining edges randomly
		while (edges_cnt > 0) {
			out_idx = rnd.nextInt(all_nodes.size());
			label = st.getLabel();
			while (true) {
				in_idx = rnd.nextInt(all_nodes.size());
				if (in_idx != out_idx) { 
					node_1 = all_nodes.get(out_idx);
					node_2 = all_nodes.get(in_idx);
					if (st.randomDirectedEdge(node_1, node_2, label)) {
						edges_cnt--;
						rnd_cnt++;
						System.out.println("Random edge added: " + (int)(edges-edges_cnt) + "/" + edges);
						break;
					} else {
						waste_labels.add(label);
						label = st.randomLabel();
						System.out.println("Adding random edge failed!");
					}
				}
			}
		}
		
		System.out.println("Number of communities: " + cmty_refs.size());
		System.out.println("Minimum community size: " + mcs);
		System.out.println("Maximum number of communities: " + max_cmty_cnt);
		
		System.out.println("Community edges added: " + cmty_cnt);
		System.out.println("Community rest edges added: " + cmty_rest_cnt);
		System.out.println("Community interconnective edges added: " + cmty_intr_cnt);
		System.out.println("Random rest edges added: " + rnd_cnt);
		
		boolean valid = true;
		for (Vertex v : all_nodes) {
			for (String s : all_labels) {
				int lod = lc_map.get(new Pair<Vertex,String>(v,s));
				if (lc.get(1).containsKey(s)) {
					if (lod > lc.get(1).get(s)) {
						valid = false;
					}
				}
			}
		}
		
		if (valid) {
			System.out.println("Label constraints were not violated!");
		} else {
			System.out.println("VIOLATION OF LABEL CONSTRAINTS!");
		}
	}
	
	
	/**
	 * Erases existing edges of graph, but not their indexes. 
	 * For multiple synthetizations, this function takes longer 
	 * time due to iterating over more indexes each time.
	 */
	public void eraseEdges() {
		System.out.print("Clearing edge collection...");
		edgeCol.truncate();
		System.out.println("done!");
	}
	
	
	/**
	 * Generates an iterator over triples in a RDF-file.
	 * @param name - Name of RDF-file.
	 * @return iterator
	 */
	public PipedRDFIterator<Triple> loadRDF(String name) {
		final String file = "resources/" + name + ".rdf";
		PipedRDFIterator<Triple> iter = new PipedRDFIterator<>();
		final PipedRDFStream<Triple> inputStream = new PipedTriplesStream(iter);
		ExecutorService executor = Executors.newSingleThreadExecutor();
		
		Runnable parser = new Runnable() {
			@Override
			public void run() {
				RDFParser.source(file).parse(inputStream);;
			}
		};
		
		executor.submit(parser);
		return iter;
	}
	
	
	/**
	 * Parses a RDF-file to a graph structure in database.
	 * @param name - Name of RDF-file.
	 */
	public void parseRDF(String name) {
		Map<Node,Vertex> graph_nodes = new HashMap<Node,Vertex>();
		PipedRDFIterator<Triple> iter = loadRDF(name);
		
		while (iter.hasNext()) {
			Triple next = iter.next();
			
			Node ns = next.getMatchSubject();
			Node no = next.getMatchObject();
			Vertex v1 = graph_nodes.get(ns);
			Vertex v2 = graph_nodes.get(no);
			
			if (v1 == null) {
				v1 = g.addV("node").property(T.id, Integer.toString(tid))
						.property("name", ns.toString()).next();
				cloneGTS();
				tid += 2;
			} 
			
			if (v2 == null) {
				v2 = g.addV("node").property(T.id, Integer.toString(tid))
						.property("name", no.toString()).next();
				cloneGTS();
				tid += 2;
			}
			
			graph_nodes.put(ns, v1);
			graph_nodes.put(no, v2);
				
			g.addE("edge").from(v1).to(v2).property("label", next.getMatchPredicate().toString()).next();
			cloneGTS();
		}
	}
	
	
	/**
	 * Query in- and out-degree of node (vertex).
	 * @param vx - Node in scope.
	 * @return
	 */
	private Long nodeDegree(Vertex vx) {
		return nodeInDegree(vx) + nodeOutDegree(vx);
	}
	
	
	/**
	 * Query out-degree of node (vertex).
	 * @param vx - Node in scope.
	 * @return
	 */
	private long nodeOutDegree(Vertex vx) {
		return g.V(vx.id()).outE().count().next();
	}
	
	
	/**
	 * Query in-degree of node (vertex).
	 * @param vx - Node in scope.
	 * @return
	 */
	private long nodeInDegree(Vertex vx) {
		return g.V(vx.id()).inE().count().next();
	}
	
	
	/**
	 * Query labeled in- and out-degree of a node (vertex).
	 * @param vx - Node in scope.
	 * @param labels - Labels included in the query.
	 * @return
	 */
	private long labelNodeDegree(Vertex vx, String[] labels) {
		return labelInDegree(vx, labels) + labelOutDegree(vx, labels);
	}
	
	
	/**
	 * Query labeled in-degree of node (vertex).
	 * @param vx - Node in scope.
	 * @param labels - Labels included in the query.
	 * @return
	 */
	private long labelInDegree(Vertex vx, String[] labels) {
		long result = 0;
		Iterator<Edge> it = g.V(vx.id()).inE();
		while (it.hasNext()) {
			if (containsLabel(labels, it.next().value("label").toString())) {
				result++;
			}
		}
		return result;
	}
	
	
	/**
	 * Query labeled out-degree of node (vertex).
	 * @param vx - Node in scope.
	 * @param labels - Labels included in the query.
	 * @return
	 */
	private long labelOutDegree(Vertex vx, String[] labels) {
		long result = 0;
		Iterator<Edge> it = g.V(vx.id()).outE();
		while (it.hasNext()) {
			if (containsLabel(labels, it.next().value("label").toString())) {
				result++;
			}
		}
		return result;
	}
	
	
	/**
	 * Checks if an edge label exists in the list of 
	 * user-provided query labels.
	 * @param labels - Labels included in the query.
	 * @param label - Edge label.
	 * @return
	 */
	private boolean containsLabel(String[] labels, String label) {
		for (String str : labels) {
			if (str.equals(label)) {
				return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Help function for Q1. Retrieves maximum node degree from a graph.
	 * @return max node degree.
	 */
	public long queryMaxDegree() {
		long max_deg = 0;
		long node_deg = 0;
		ArrayList<Vertex> nodes = getNodes();
		for (Vertex vx : nodes) {
			node_deg = nodeDegree(vx);
			max_deg = Math.max(max_deg, node_deg);
		}
		return max_deg;
	}
	
	
	/**
	 * Help function for queries that analyzes degree distribution. 
	 * Given query q, in-going, out-going, label-specific in-going 
	 * or label-specific out-going node degree is determined.
	 * @param q - Query.
	 * @param vx - Node in scope.
	 * @param labels - Labels included in query q.
	 * @return Node degree.
	 */
	private long nodeDegreeType(Query q, Vertex vx, String[] labels) {
		long result = 0;
		
		if (q == Query.Q2 || q == Query.Q5) {
			if (labels == null) {
				result = nodeDegree(vx);
			} else {
				result = labelNodeDegree(vx, labels);
			}
		} else if (q == Query.Q3 || q == Query.Q6) {
			if (labels == null) {
				result = nodeInDegree(vx);
			} else {
				result = labelInDegree(vx, labels);
			}
		} else if (q == Query.Q4 || q == Query.Q7) {
			if (labels == null) {
				result = nodeOutDegree(vx);
			} else {
				result = labelOutDegree(vx, labels);
			}
		}
		
		return result;
	}
	
	
	/**
	 * Function for iteration over all nodes in a graph for 
	 * determining node degrees for all nodes. Node degrees 
	 * are stored as a histogram in a double[], where the 
	 * value interval for each bucket is determined given 
	 * the user-provided input b and maximum node degree.
	 * @param q - Query.
	 * @param b - Number of buckets in histogram.
	 * @param bw - Bucket width.
	 * @param labels - Labels included in query q.
	 * @return Histogram.
	 */
	private ArrayList<Long> queryNodeDegrees(Query q, int b, int bw, String[] labels, boolean seq) {
		ArrayList<Vertex> nodes = getNodes();
		HashMap<Vertex,Long> node_degrees = new HashMap<Vertex,Long>();
		
		int previous = 0;
		int progress;
		double cnt = 0;
		double nodes_size = nodes.size();
		for (Vertex vx : nodes) {
			long res = nodeDegreeType(q, vx, labels);
			node_degrees.put(vx, res);
			cnt++;
			progress = (int) (Math.floor(100.0 * (cnt / nodes_size)));
			if (progress % 5 == 0) {
				if (progress != previous) {
					previous = progress;
					System.out.print(progress + "% ");
				}
			}
		}
		System.out.println();

		if (seq) {
			return retSequence(node_degrees);
		}
		return retHistogram(node_degrees, b, bw);
	}
	
	
	/**
	 * Help function for the triangle count algorithm. 
	 * Returns a list of neighbors to the node vx.
	 * @param vx - Node.
	 * @return List of neighbors (adjacent nodes).
	 */
	private ArrayList<Vertex> directedNeighbors(Vertex vx, Direction d) {
		ArrayList<Vertex> neighbors = new ArrayList<Vertex>();
		Iterator<Vertex> it = vx.vertices(d);
		while(it.hasNext()) {
			neighbors.add(it.next());
		}
		return neighbors;
	}
	
	
	/**
	 * Triangle count algorithm (compact-forward) that 
	 * counts the triangle involvement for each node in the 
	 * graph. Works for undirected and directed triangles, 
	 * given the value on input parameter.
	 * @param d - Direction (BOTH,IN,OUT)
	 * @return
	 */
	private ArrayList<ArrayList<Vertex>> cfTriangleCount(Direction d) {
		ArrayList<Vertex> nodes = getNodes();
		ArrayList<ArrayList<Vertex>> triangle_sets = new ArrayList<ArrayList<Vertex>>();
		Iterator<Vertex> uit,vit;
		Vertex u_prime,v_prime;
		Direction d1 = Direction.BOTH;
		Direction d2 = Direction.BOTH;
		
		if (d == Direction.IN) {
			d1 = Direction.IN;
			d2 = Direction.OUT;
		} else if (d == Direction.OUT) {
			d1 = Direction.OUT;
			d2 = Direction.IN;
		} else {
			System.out.println("Triangle count - BOTH directions selected.");
		}
		
		Comparator<Vertex> comp = new Comparator<Vertex>() {
			@Override
			public int compare(Vertex v1, Vertex v2) {
				return nodeDegree(v2).compareTo(nodeDegree(v1));
			}
		};
		
		Comparator<Vertex> comp_n = new Comparator<Vertex>() {
			@Override
			public int compare(Vertex v1, Vertex v2) {
				Integer i1 = nodes.indexOf(v1);
				Integer i2 = nodes.indexOf(v2);
				return i1.compareTo(i2);
			}
		};
		
		Collections.sort(nodes, comp);
		
		for (Vertex v : nodes) {
			ArrayList<Vertex> neighbors = directedNeighbors(v, d1);
			Collections.sort(neighbors, comp_n);
			for (Vertex u : neighbors) {
				if (nodes.indexOf(u) > nodes.indexOf(v)) {
					ArrayList<Vertex> uNeighbors = directedNeighbors(u, d1);
					ArrayList<Vertex> vNeighbors = directedNeighbors(v, d2);
					Collections.sort(uNeighbors, comp_n);
					Collections.sort(vNeighbors, comp_n);
					uit = uNeighbors.iterator();
					vit = vNeighbors.iterator();
					if (!(uit.hasNext() && vit.hasNext())) {
						continue;
					}
					u_prime = uit.next();
					v_prime = vit.next();
					
					while (true) {
						if (nodes.indexOf(u_prime) < nodes.indexOf(v) &&
								nodes.indexOf(v_prime) < nodes.indexOf(v)) {
							if (nodes.indexOf(u_prime) < nodes.indexOf(v_prime)) {
								if (!uit.hasNext()) {
									break;
								}
								u_prime = uit.next();
							} else if (nodes.indexOf(u_prime) > nodes.indexOf(v_prime)) {
								if (!vit.hasNext()) {
									break;
								}
								v_prime = vit.next();
							} else {
								System.out.println("{v1,v2,v3}: " + u.id() + "," +  v.id() + "," + u_prime.id());
								ArrayList<Vertex> foundSet = new ArrayList<Vertex>();
								foundSet.add(u);
								foundSet.add(v);
								foundSet.add(u_prime);
								Collections.sort(foundSet, comp_n);
								triangle_sets.add(foundSet);
								
								if (vit.hasNext() && uit.hasNext()) {
									v_prime = vit.next();
									u_prime = uit.next();
								} else {
									break;
								}
							}
						} else {
							break;
						}
					}
				}
			}
		}
		
		return triangle_sets;
	}
	
	
	/**
	 * Calls the triangle count algorithm, processes the result 
	 * and returns as histogram or sequence, depending on the query.
	 * (Maybe) TODO: Implement logic for directed triangle count query.
	 * @param q - Query.
	 * @param b - Number of buckets in histogram.
	 * @param bw - Bucket width.
	 * @return
	 */
	private ArrayList<Long> triangleCountQuery(Query q, int b, int bw) {
		ArrayList<ArrayList<Vertex>> triangle_sets = cfTriangleCount(Direction.BOTH);
		HashMap<Vertex,Long> triangle_counts = new HashMap<Vertex,Long>();
		ArrayList<Vertex> nodes = getNodes();
		
		for (Vertex vx : nodes) {
			triangle_counts.put(vx, (long) 0);
		}
		
		for (ArrayList<Vertex> ts : triangle_sets) {
			for (Vertex vx : ts) {
				long previous = triangle_counts.get(vx);
				triangle_counts.put(vx, previous+1);
			}
		}
		
		//System.out.println(triangle_sets.toString());
		//System.out.println(triangle_counts.toString());
		
		if (q == Query.Q9) {
			return retSequence(triangle_counts);
		}
		return retHistogram(triangle_counts, b, bw);
	}
	
	
	/**
	 * Similar to the function above, only that the edge direction 
	 * when counting the triangles is considered.
	 * TODO: Look into optimization? (takes ~9 minutes for a paper-generated graph)
	 * @return
	 */
	private HashMap<Vertex,Long> directedTriangleCountQuery() {
		ArrayList<ArrayList<Vertex>> triangle_sets_in = cfTriangleCount(Direction.IN);
		ArrayList<ArrayList<Vertex>> triangle_sets_out = cfTriangleCount(Direction.OUT);
		HashSet<ArrayList<Vertex>> triangle_sets = new HashSet<ArrayList<Vertex>>();
		HashMap<Vertex,Long> triangle_counts = new HashMap<Vertex,Long>();
		ArrayList<Vertex> nodes = getNodes();
		
		for (Vertex vx : nodes) {
			triangle_counts.put(vx, (long) 0);
		}
		
		for (ArrayList<Vertex> ts : triangle_sets_in) {
			triangle_sets.add(ts);
		}

		for (ArrayList<Vertex> ts : triangle_sets_out) {
			triangle_sets.add(ts);
		}

		for (ArrayList<Vertex> ts : triangle_sets) {
			for (Vertex vx : ts) {
				long previous = triangle_counts.get(vx);
				triangle_counts.put(vx, previous+1);
			}
		}
		
		return triangle_counts;
	}
	
	
	/**
	 * Determines how many nodes are in specific intervals 
	 * of local clustering coefficients. Intervals are determined 
	 * by the value on k (which equals number of buckets for 
	 * histograms, although the query is not approached as a 
	 * histogram query). The resulting array is a list of the 
	 * amount of nodes within certain intervals, where the 
	 * possible range for a coefficient is [0,1].
	 * @param k - Number of intervals.
	 * @return
	 */
	private ArrayList<Long> localClusterCoeff(int k) {
		HashMap<Vertex,Long> triangle_counts = directedTriangleCountQuery();
		HashMap<Vertex,Double> clustering_coeffs = new HashMap<Vertex,Double>();
		
		for (Map.Entry<Vertex,Long> entry : triangle_counts.entrySet()) {
			long lambdaV = directedNeighbors(entry.getKey(), Direction.OUT).size() * directedNeighbors(entry.getKey(), Direction.IN).size();
			//long lambdaV = nodeOutDegree(entry.getKey()) * nodeInDegree(entry.getKey());	// TODO: Not right, should be out-neighbors * in-neighbors - directedNeighbors(entry.getKey(), Direction.?)
			double lcc;
			if (lambdaV == 0) {
				lcc = 0;
			} else {
				lcc = entry.getValue() / (double) lambdaV;
			}
			clustering_coeffs.put(entry.getKey(), lcc);
		}
		
		ArrayList<Long> lccd = new ArrayList<Long>();
		
		for (int i = 0; i < k; i++) {
			lccd.add((long) 0);
		}
		
		for (Map.Entry<Vertex,Double> entry : clustering_coeffs.entrySet()) {
			int idx = (int) Math.floor(entry.getValue()*k);
			if (idx == k) { 
				idx--; 
			}
			long prev = lccd.get(idx);
			lccd.set(idx, ++prev);
		}
		
		return lccd;
	}
	
	
	/**
	 * Function for creating a histogram out of the query result. 
	 * @param res - Map that maps each node to the resulting value.
	 * @param b - Number of buckets.
	 * @param bw - Bucket width.
	 * @return
	 */
	private ArrayList<Long> retHistogram(HashMap<Vertex,Long> res, int b, int bw) {
		ArrayList<Long> hg = new ArrayList<Long>();
		int idx; 
		
		for (int i = 0; i < b; i++) {
			hg.add((long) 0);
		}
		
		for (Map.Entry<Vertex, Long> entry : res.entrySet()) {
			double val = entry.getValue()/bw;
			idx = (int) Math.floor(val);
			if (idx > b-1) { idx = b-1; }
			long prev = hg.get(idx);
			hg.set(idx, ++prev);
		}
		return hg;
	}
	
	
	/**
	 * Function for creating a sequence out of the query result.
	 * @param res - Map that maps each node to the resulting value.
	 * @return
	 */
	private ArrayList<Long> retSequence(HashMap<Vertex,Long> res) {
		ArrayList<Long> seq = new ArrayList<Long>();
		
		for (Map.Entry<Vertex, Long> entry : res.entrySet()) {
			seq.add(entry.getValue());
		}
		
		Collections.sort(seq);
		return seq;
	}
	
	
	/**
	 * Determines which query to be sent to database.
	 * @param q - Query.
	 * @param b - Number of buckets in histogram.
	 * @param bw - Bucket width.
	 * @param labels - Labels included in query q.
	 * @return True query result.
	 */
	public ArrayList<Long> queryDb(Query q, int b, int bw, String[] labels) {
		ArrayList<Long> query_result = new ArrayList<Long>();
		System.out.print("Querying database... ");
		
		// Q1
		if (q == Query.Q1) {
			query_result.add(queryMaxDegree());
		// Q2-Q4
		} else if (q.ordinal() > Query.Q1.ordinal() && q.ordinal() < Query.Q5.ordinal()) {
			query_result = queryNodeDegrees(q, b, bw, labels, false);
		// Q5-Q7
		} else if (q.ordinal() > Query.Q4.ordinal() && q.ordinal() < Query.Q8.ordinal()) {
			query_result = queryNodeDegrees(q, b, bw, labels, true);
		// Q8-Q9
		} else if (q == Query.Q8 || q == Query.Q9) {
			query_result = triangleCountQuery(q, b, bw);
		// Q10
		} else if (q == Query.Q10) {
			query_result = localClusterCoeff(b);
		}
		return query_result;
	}
	
	
	/**
	 * Attempts to close both graph traversal source and 
	 * the graph itself.
	 */
	public void closeGraph() {
		try {
			gts.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		try {
			graph.close();
		} catch (Exception e) {
			e.printStackTrace();
		} 
	}
}
