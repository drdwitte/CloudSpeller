# CloudSpeller

Latest version of CloudSpeller uses maven to build the project. 

Installation instructions: 


mkdir projectDirectory
git clone https://github.com/drdwitte/CloudSpeller.git
cd CloudSpeller
mvn clean install

Now the target directory will contain 2 jar files:

original-CloudSpeller-1.0.jar is a skinny jar file containing the Speller code only, you will need to add Hadoop to the classpath
CloudSpeller-1.0.jar is a fat jar which has the hadoop libraries included in the jar file, this can run in stand-alone mode
