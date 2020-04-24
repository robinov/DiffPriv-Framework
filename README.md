# This framework was coded with Eclipse Java 2018-12 and the instructions below are based on that version of the IDE.
1. Download or git-clone the project
2. Convert to Maven project
3. Java Compiler -> Uncheck 'Use compliance from execution env..."
4. Java Compiler -> Uncheck 'Use '--release' option'
5. Set compliance level: 1.8
6. Java Build Path -> JRE System Library -> Edit... -> Set execution environment to 1.8
7. Set up your Arango database and fill in the configuration under the ArangoBridge class
8. Should be good to go!

# Framework description
This framework requires a local setup of ArangoDB in order to be used. 
Once the setup is ready, the framework can query existing graphs and generate 
noisy results for which utility is measured and saved for output. The framework 
generates two output files in query-mode, where one contains a summary of true 
and noisy results as well as average utility per epsilon value provided in the 
config.properties file (and standard deviation for each average value). The other 
output file contains the the utility measurements that was made for as many iterations 
as specified by the parameter in the configuration file (NL).

# How-to-use
Firstly a reminder that any type of execution requires that a local setup of ArangoDB 
is up. How the framework then is executed is based on the content of the output/config.properties 
file. If set to synthesize (SYNTHESIZE=true), there needs to be an existing label-constraints 
file which can be created manually, or with the runGLC() function in the Main class (replace 
run() in the main function). The result of the function ends up in the resources/label_constraints 
folder, which does not need to be moved at all. With a valid input on the LC-parameter (name 
of the label-constraints file), the framework can synthesize a graph. If a graph already 
exist with the specified name (GRAPH=title), that graph needs to be manually removed from 
ArangoDB before continuing, as there currently is no way to ensure the creation of a whole new 
graph in those cases. 

When not synthesizing, the framework is in query-mode. In query-mode, two output files will 
be generated according to OUTPUT=filename (filename_SUM & filename_UTIL), which are placed 
in output/summary/ and output/utility/ respectively. It is the output placed in output/utility/ 
that was used to generated initial boxplots for the paper. The remaining parameters in the 
config.properties file speak for themselves, except for NL and AL. NL is how many times noisy 
query results should be generated (and utility therefore measured). AL is specifying the set 
of query labels, which is either 'ALL' or [l1,l2,l3]. To find out the unique set of labels 
of a graph, the function extractLabelsToFile(graph_name) can be used. The set of unique 
labels are placed in a text-file in output/labels/.
