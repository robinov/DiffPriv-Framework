package framework;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;
import java.util.Scanner;
import java.util.concurrent.TimeUnit;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.tinkerpop.gremlin.structure.Edge;


public class Main {
	private Logger logger;
	private Handler windowHandler;
	private HashMap<String,Query> queryCodes;
	
	/**
	 * Constructor for Main-class for holding Logger and Window.
	 */
	public Main() {
		logger = Logger.getLogger(Main.class.getName());
		windowHandler = new WindowHandler();
		logger.addHandler(windowHandler);
		queryCodes = setupCodes();
	}
	
	
	public HashMap<String,Query> setupCodes() {
		HashMap<String,Query> codes = new HashMap<String,Query>();
		
		codes.put("MD", Query.Q1);		// Maximum degree
		codes.put("DDH", Query.Q2);		// Degree distribution [Histogram]
		codes.put("IDDH", Query.Q3);	// In-degree distribution [Histogram]
		codes.put("ODDH", Query.Q4);	// Out-degree distribution [Histogram]
		codes.put("DDS", Query.Q5);		// Degree distribution [Sequence]
		codes.put("IDDS", Query.Q6);	// In-degree distribution [Sequence]
		codes.put("ODDS", Query.Q7);	// Out-degree distribution [Sequence]
		codes.put("TCH", Query.Q8);		// Triangle count [Histogram]
		codes.put("TCS", Query.Q9);		// Triangle count [Sequence]
		codes.put("LCCD", Query.Q10);	// Local clustering coefficient distribution [0:1]
		
		return codes;
	}
	
	
	/**
	 * Takes an ArrayList of Longs and converts it to a list of doubles.
	 * @param input - ArrayList.
	 * @return double[]
	 */
	public static double[] convertToDoubleList(ArrayList<Long> input) {
		double[] output = new double[input.size()];
		for (int i = 0; i < output.length; i++) {
			output[i] = input.get(i);
		}
		return output;
	}
	
	
	/**
	 * Takes an ArrayList of Strings and converts it to a list of Strings.
	 * @param input
	 * @return
	 */
	public static String[] convertToStringList(ArrayList<String> input) {
		String[] output = new String[input.size()];
		for (int i = 0; i < output.length; i++) {
			output[i] = input.get(i);
		}
		return output;
	}
	
	
	/**
	 * Sends a query to a database once. Prints result.
	 * @param name - Name of graph.
	 * @param q - Query.
	 * @param labels - Labels included in query q.
	 * @param b - Number of buckets.
	 * @param bw - Bucket width.
	 */
	public static void queryOnce(String name, Query q, String[] labels, int b, int bw) {
		ArangoBridge ab = new ArangoBridge(name, true);
		ArrayList<Long> true_results_list = ab.queryDb(q, b, bw, labels);
		for (int i = 1; i <= true_results_list.size(); i++) {
			System.out.println("b_" + i + ": " + true_results_list.get(i));
		}
	}
	
	
	public void logMessage(LogRecord lr) {
		logger.log(lr);
	}
	
	
	/**
	 * Help function, String -> double.
	 * @param s - String to convert.
	 * @return
	 */
	public static double strToDouble(String s) {
		return Double.parseDouble(s);
	}
	
	
	/**
	 * Help function, String[] -> double[].
	 * @param sl - List of String to convert.
	 * @return
	 */
	public static double[] strsToDoubles(String[] sl) {
		double[] res = new double[sl.length]; 
		for (int i = 0; i < sl.length; i++) {
			res[i] = strToDouble(sl[i]);
		}
		return res;
	}
	
	
	/**
	 * Synthesizes a graph based on nodes, edges, initial maximum degree 
	 * and a label distribution.
	 * @param nodes - Number of nodes.
	 * @param edges - Number of edges.
	 * @param mcs - Minimum community size.
	 * @param nsd - Standard deviation of normal distribution.
	 */
	public static void runSynth(String name, int nodes, int edges, int mcs, ArrayList<HashMap<String,Integer>> ld, String lc_filename) {
		int nol = ld.get(0).size();
		ArangoBridge ab = new ArangoBridge(name, false);
		ab.synthesize(nodes, edges, nol, mcs, ld);
		ab.closeGraph();
	}
	
	
	/**
	 * Reads from a label distribution file (csv-format) and parses each 
	 * triple into two maps of lower and upper bounds on number of specifically 
	 * labeled outgoing edges per node.
	 * @param filename - Name of label distribution file.
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<HashMap<String,Integer>> getLabelConstraints(String filename) throws IOException {
		ArrayList<String[]> lc = new ArrayList<String[]>();
		String lc_input = "resources/label_constraints/" + filename + ".csv";
		String line = "";
		BufferedReader br = new BufferedReader(new FileReader(lc_input));
		while ((line = br.readLine()) != null) {
			String[] triple = line.split(";");
			lc.add(triple);
		}
		br.close();
		
		HashMap<String,Integer> lower_bounds = new HashMap<String,Integer>();
		HashMap<String,Integer> upper_bounds = new HashMap<String,Integer>();
		
		for (String[] tr : lc) {
			int lb;
			if (NumberUtils.isParsable(tr[1])) {
				lb = Integer.parseInt(tr[1]);
				if (lb >= 0) {
					lower_bounds.put(tr[0], lb);
				}
			} else {
				System.out.println("ERROR! LOWER BOUND WAS NOT PARSABLE!");
				return null;
			}
			
			if (tr[2].toLowerCase().equals("inf")) {
				continue;
			} else if (NumberUtils.isParsable(tr[2])) {
				int ub = Integer.parseInt(tr[2]);
				if (ub >= lb) {
					upper_bounds.put(tr[0], ub);
				}
			} else {
				System.out.println("ERROR! UPPER BOUND WAS NOT PARSABLE!");
				return null;
			}
		}
		
		ArrayList<HashMap<String,Integer>> label_bounds = new ArrayList<HashMap<String,Integer>>();
		label_bounds.add(lower_bounds);
		label_bounds.add(upper_bounds);
		return label_bounds;
	}
	
	
	/**
	 * Main function containing the actual running sequence, 
	 * which is the intended use of the framework.
	 */
	public static void run() {
		Main m = new Main();
		m.logMessage(new LogRecord(Level.INFO, "Application process starting..."));
		
		Properties p = new Properties();
		InputStream input;
		
		// read(param)
		m.logMessage(new LogRecord(Level.INFO, "Loading configuration..."));
		try {
			input = new FileInputStream("resources/config.properties");
			p.load(input);
			m.logMessage(new LogRecord(Level.INFO, "Configuration loaded successfully!"));
		} catch (IOException e) {
			e.printStackTrace();
			m.logMessage(new LogRecord(Level.SEVERE, "Loading configuration failed!"));
		}
		
		boolean synthesize;
		synthesize = Boolean.parseBoolean(p.getProperty("SYNTHESIZE"));
		String graph_name;
		
		if (synthesize) {
			m.logMessage(new LogRecord(Level.INFO, "synthesize = true"));
			int nodes = Integer.parseInt(p.getProperty("NODES"));
			int edges = Integer.parseInt(p.getProperty("EDGES"));
			int mcs = Integer.parseInt(p.getProperty("MCS"));
			String lc_filename = p.getProperty("LC");
			graph_name = p.getProperty("GRAPH");
			try {
				ArrayList<HashMap<String,Integer>> lc = getLabelConstraints(lc_filename);
				if (lc == null) {
					m.logMessage(new LogRecord(Level.SEVERE, "Label distribution had invalid data!"));
					return;
				}
				long start = System.nanoTime();
				runSynth(graph_name, nodes, edges, mcs, lc, lc_filename);
				long end = System.nanoTime();
				long duration = end - start;
				System.out.println("Synthetization took " + (double) TimeUnit.NANOSECONDS.toSeconds(duration) + " seconds");
			} catch (IOException e) {
				e.printStackTrace();
			}
			return;
		} 
		m.logMessage(new LogRecord(Level.INFO, "synthesize = false"));
			
		int nl, b, bw, sens;
		String out_name, query_code;
		double[] eps;
		String[] al;
		Query q;
		
		nl = Integer.parseInt(p.getProperty("NL"));
		b = Integer.parseInt(p.getProperty("BUCKETS"));
		bw = Integer.parseInt(p.getProperty("BUCKET_WIDTH"));
		sens = Integer.parseInt(p.getProperty("SENS"));
		out_name = p.getProperty("OUTPUT");
		graph_name = p.getProperty("GRAPH");
		query_code = p.getProperty("QUERY").toUpperCase();
		eps = strsToDoubles(p.getProperty("EPSILON").split("\\s*,\\s*"));
		String al_str = p.getProperty("AL");

		q = m.queryCodes.get(query_code);
		
		if (al_str.startsWith("[") && al_str.endsWith("]")) {
			al = al_str.substring(1, al_str.length()-1).split("\\s*,\\s*");
			m.logMessage(new LogRecord(Level.INFO, "Parsed " + al.length + " items from label parameter!"));
		} else if (al_str.toLowerCase().equals("all")) {
			al = null;
		} else {
			m.logMessage(new LogRecord(Level.SEVERE, "Error. Make sure AL matches either [label1,label2,...] or \"all\"."));
			return;
		}
		
		String name = graph_name;
		m.logMessage(new LogRecord(Level.INFO, "Query selected: " + q.toString()));
		
		ArangoBridge ab = new ArangoBridge(name, !synthesize);
		QueryMechanism qm = new QueryMechanism();
		double[] noisy_result;
		double[] average_utility = new double[eps.length];
		ArrayList<double[]> out_query = new ArrayList<double[]>();
		ArrayList<double[]> out_analysis = new ArrayList<double[]>();
		ArrayList<double[]> all_utility_ratios = new ArrayList<double[]>();
		double[] utility_ratios_per_epsilon = new double[eps.length];	// Total Utility Per Epsilon
		boolean calc_dev = false;
		
		double[] true_result = convertToDoubleList(ab.queryDb(q, b, bw, al));
		out_query.add(true_result);
		m.logMessage(new LogRecord(Level.INFO, "Query complete!"));
		
		RealVector trv = new ArrayRealVector(true_result);
		RealVector nrv;
		for (int i = 0; i < eps.length; i++) {
			noisy_result = qm.generateNoise(true_result, sens, eps[i]);
			out_query.add(noisy_result);
			nrv = new ArrayRealVector(noisy_result);
			utility_ratios_per_epsilon[i] = nrv.getL1Distance(trv);
		}
		all_utility_ratios.add(utility_ratios_per_epsilon);
		
		int iterations = nl+1;
		m.logMessage(new LogRecord(Level.INFO, "Generating noise..."));
		
		if (nl > 1) { 
			calc_dev = true; 
			nl--;
		}
		
		// Calculate the utility for each noisy result
		while (nl > 0) { 
			utility_ratios_per_epsilon = new double[eps.length];
			for (int i = 0; i < eps.length; i++) {
				noisy_result = qm.generateNoise(true_result, sens, eps[i]);
				nrv = new ArrayRealVector(noisy_result);
				utility_ratios_per_epsilon[i] = nrv.getL1Distance(trv);
			}
			all_utility_ratios.add(utility_ratios_per_epsilon);
			nl--;
		}
		
		// Calculate total utility for each epsilon
		for (int i = 0; i < all_utility_ratios.size(); i++) {
			for (int j = 0; j < eps.length; j++) {
				average_utility[j] += all_utility_ratios.get(i)[j];
			}
		}
		
		// Average total utility
		for (int i = 0; i < average_utility.length; i++) {
			average_utility[i] /= iterations;
		}
		
		double[] sd_utility = new double[eps.length];
		for (int i = 0; i < sd_utility.length; i++) {
			sd_utility[i] = 0.0;
		}
		
		if (calc_dev) {
			// Standard deviation
			for (int i = 0; i < all_utility_ratios.size(); i++) {
				for (int j = 0; j < eps.length; j++) {
					sd_utility[j] += Math.pow(all_utility_ratios.get(i)[j] - average_utility[j], 2.0);
				}
			}
			for (int i = 0; i < sd_utility.length; i++) {
				sd_utility[i] = Math.sqrt(sd_utility[i] / iterations);
			}
		}
		
		out_analysis.add(eps);
		out_analysis.add(average_utility);
		out_analysis.add(sd_utility);
		m.logMessage(new LogRecord(Level.INFO, "Noise generation complete!"));
		
		try {
			qm.genSumOutput(out_query, out_analysis, out_name + "_SUM", eps, true_result.length);
			m.logMessage(new LogRecord(Level.INFO, "Output file: \"" + out_name + "_SUM.csv\" was generated!"));
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			qm.genUtilOutput(all_utility_ratios, out_name + "_UTIL", eps);
			m.logMessage(new LogRecord(Level.INFO, "Output file: \"" + out_name + "_UTIL.csv\" was generated!"));
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ab.closeGraph();
		m.logMessage(new LogRecord(Level.INFO, "Process is finished!"));
	}
	
	
	/**
	 * Exports a list of unique elements of all labels 
	 * in a graph to an output file.
	 * @param synthesized - true if the graph was synthesized.
	 */
	public static void extractLabelsToFile(String name) {
		
		ArangoBridge ab = new ArangoBridge(name,false);
		ArrayList<Edge> graph_edges = ab.getEdges();
		ArrayList<String> labels = new ArrayList<String>();
		String label;
		
		for (Edge edge : graph_edges) {
			label = edge.value("label").toString();
			if (!labels.contains(label)) {
				labels.add(label);
			}
		}
		
		try {
			FileWriter fw = new FileWriter("output/labels/labels_" + name + ".txt");
			for (String str : labels) {
				fw.write(str + "\n");
			}
			fw.close();
		} catch(Exception e) {
			System.out.println(e);
		}
	}
	
	
	/**
	 * Produces a list of unique elements of all labels in a graph.
	 * @param name - Name of graph.
	 * @return list of labels.
	 */
	public static String[] extractLabels(String name) {
		ArangoBridge ab = new ArangoBridge(name,false);
		ArrayList<Edge> graph_edges = ab.getEdges();
		ArrayList<String> labels_arr = new ArrayList<String>();
		String label;
		
		for (Edge edge : graph_edges) {
			label = edge.value("label").toString();
			if (!labels_arr.contains(label)) {
				labels_arr.add(label);
			}
		}
		
		String[] labels = convertToStringList(labels_arr);
		
		return labels;
	}
	
	
	/**
	 * Help function to assist as a "contains" function for int[]
	 * @param list - list of integers.
	 * @param item - value to check for.
	 * @return
	 */
	public static boolean includes(int[] list, int item) {
		for (int i : list) {
			if (i == item) {
				return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Call function for generating label constraints file.
	 * @param filename - Name of the label constraints file.
	 */
	public static void runGLC(String filename) {
		try {
			genLabelConstraints(filename);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Generates label constraints in the form of a csv file, where each 
	 * row is on the format [L;X;Y], where X is the minimum amount of 
	 * outgoing edges with label L each node in the graph must have, and 
	 * Y is the maximum amount of such edges.
	 * @param filename - Name of the label constraints file.
	 * @throws IOException
	 */
	public static void genLabelConstraints(String filename) throws IOException{
		FileWriter csv = new FileWriter("resources/label_constraints/" + filename + ".csv");
		
		int l = 1;
		String x,y;
		Scanner scn = new Scanner(System.in);
		while(true) {
			String label = "label" + l;
			System.out.print("Enter x and y for " + label + ": ");
			String iline = scn.nextLine();
			if (iline.isEmpty()) {
				break;
			}
			String[] str_vals = iline.split(" ");
			x = str_vals[0];
			y = str_vals[1];
			csv.append(label);
			csv.append(";");
			csv.append(x);
			csv.append(";");
			csv.append(y);
			csv.append("\n");
			System.out.println("Formed triplet: [" + label + "," + x + "," + y + "]");
			l++;
		}
		scn.close();
		csv.flush();
		csv.close();
	}

	
	public static void main(String[] args) {
		long start = System.nanoTime();
		run();
		long end = System.nanoTime();
		long duration = end - start;
		System.out.println("Execution took " + (double) TimeUnit.NANOSECONDS.toSeconds(duration) + " seconds");
		System.exit(1);
	}
}
