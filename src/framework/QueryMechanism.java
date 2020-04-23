package framework;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.math3.distribution.LaplaceDistribution;

public class QueryMechanism {
	
	
	/**
	 * Previously had data fields for query sensitivities. 
	 * Now this class is only a collection of functionality 
	 * appropriate to a query mechanism (noise & output).
	 */
	public QueryMechanism() {}

	
	/**
	 * Adds noise to the true result through LaplaceDistribution(mu, beta), 
	 * where mu = center of graph (0 for accuracy) and beta = lambda, 
	 * where lambda = worst case sensitivity / epsilon.
	 * @param q_res - True query results.
	 * @param sens - Sensitivity according to query and notion.
	 * @param eps - Epsilon constant.
	 * @return Noisy results.
	 */
	public double[] generateNoise(double[] tr, int sens, double eps) { 
		double [] noisy_result = tr.clone();
		
		if (sens == 0.0) {
			return noisy_result;
		}
		
		double lambda = sens/eps;
		LaplaceDistribution ld = new LaplaceDistribution(0, lambda);
		
		for (int i = 0; i < noisy_result.length; i++) {
			noisy_result[i] += ld.sample();
		}
		
		return noisy_result;
	}
	
	
	/**
	 * Generates output file containing average and standard deviation of noise, 
	 * as well as the true result and multiple instances of noisy result 
	 * corresponding to the user-provided epsilon values.
	 * @param q_out - Query-related output data.
	 * @param a_out - Privacy & utility trade-off output data.
	 * @param filename - Name of output file.
	 * @param eps - Epsilon values.
	 * @param b - Number of buckets.
	 * @throws IOException
	 */
	public void genSumOutput(ArrayList<double[]> q_out, ArrayList<double[]> a_out, String filename, double[] eps, int b) throws IOException{
		FileWriter csv = new FileWriter("output/summary/" + filename + ".csv");
		String[] a_headers = new String[] {"eps", "uti", "sd"};
		String[] q_headers = new String[b+1];
		q_headers[0] = "eps";
		for (int i = 0; i < b; i++) {
			q_headers[i+1] = "B" + (i+1);
		}
		
		// Privacy & utility trade-off
		for (int i = 0; i < a_headers.length; i++) {
			csv.append(a_headers[i]);
			if (i != a_headers.length-1) { csv.append(";"); }
		}
		csv.append("\n");
		
		for (int i = 0; i < a_out.get(0).length; i++) {
			for (int j = 0; j < a_headers.length; j++) {
				csv.append(String.format("%.4f", a_out.get(j)[i]));
				if (j != a_headers.length-1) { csv.append(";"); }
			}
			csv.append("\n");
		}
		csv.append("\n");
		
		// True result
		for (int i = 1; i < q_headers.length; i++) {
			csv.append(q_headers[i]);
			if (i != q_headers.length-1) { csv.append(";"); }
		}
		csv.append("\n");
		
		for (int i = 0; i < q_out.get(0).length; i++) {
			csv.append(String.format("%.2f", q_out.get(0)[i]));
			if (i != q_out.get(0).length) { csv.append(";"); }
		}
		csv.append("\n\n");
		
		// Query output
		for (int i = 0; i < q_headers.length; i++) {
			csv.append(q_headers[i]);
			if (i != q_headers.length-1) { csv.append(";"); }
		}
		csv.append("\n");
		
		for (int i = 1; i < q_out.size(); i++) {
			csv.append(String.format("%.4f", eps[i-1]));
			csv.append(";");
			for (int j = 0; j < b; j++) {
				csv.append(String.format("%.4f", q_out.get(i)[j]));
				if (j != q_headers.length-1) { csv.append(";"); }
			}
			csv.append("\n");
		}
		
		csv.flush();
		csv.close();
	}
	
	
	/**
	 * Generates an output file containing the utility measurements of the 
	 * corresponding epsilon values over a number of iterations.
	 * @param res - Utility output data.
	 * @param filename - Name of output file.
	 * @param eps - Epsilon values.
	 * @throws IOException
	 */
	public void genUtilOutput(ArrayList<double[]> res, String filename, double[] eps) throws IOException {
		FileWriter csv = new FileWriter("output/utility/" + filename + ".csv");
		
		String[] headers = new String[eps.length];
		for (int i = 0; i < eps.length; i++) {
			headers[i] = "eps = " + eps[i];
		}
		
		CSVPrinter printer = new CSVPrinter(csv, CSVFormat.DEFAULT.withHeader(headers).withDelimiter(';'));
		ArrayList<String> record = new ArrayList<String>();
		
		for (int i = 0; i < res.size(); i++) {
			for (int j = 0; j < eps.length; j++) {
				record.add(String.format("%.2f", res.get(i)[j]));
			}
			printer.printRecord(record);
			record.clear();
		}
		
		printer.close();
	}
}
