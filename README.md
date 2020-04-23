# DiffPriv-Framework
# This framework was coded with Eclipse Java 2018-12 and the instructions below are based on that version of the IDE.
1. Download the project
2. Convert to Maven project
3. Java Compiler -> Uncheck 'Use compliance from execution env..."
4. Java Compiler -> Uncheck 'Use '--release' option'
5. Set compliance level: 1.8
6. Java Build Path -> JRE System Library -> Edit... -> Set execution environment to 1.8
7. Set up your Arango database and fill in the configuration under the ArangoBridge class
8. Should be good to go!