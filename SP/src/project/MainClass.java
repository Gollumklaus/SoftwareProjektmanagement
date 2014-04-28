package project;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

import javax.xml.stream.XMLStreamException;
import lpsolve.LpSolveException;

// 0.7mb H_pylori_iIT341.xml
// 1.7mb S_cerevisiae_iND750.xml
// 3.2mb Y_pestis_co92.xml

public class MainClass {

	public static void main(String args[]) throws FileNotFoundException {	
		boolean ignore = false;
		if(args.length>1){
				if("-n".equals(args[1])){
				PrintStream err = new PrintStream(args[2]);
				System.setErr(err);
				System.out.println("\"-n\" option for better readability is taken. error output is saved in: \""+ args[2]+"\"");
				ignore = true;
				}
				else{
					System.out.println("wrong function call. use: MainClass.java PATHOFSBML\n or to ignore error-output: MainClass.java PATHOFSBML -n PATHFORERROR");
				}
		}
		else{
			System.out.println("INFO: you can use option \"-n PATHOFERRORLOG\" for better readability");
		}
		
		// load SBML etc.
		SBMLLoad load = new SBMLLoad();

		if (args.length > 0) {
			try {
				System.out.println("Loading SBML file: \"" + args[0] + "\"" );
				load.loadSBML(args[0]);
			} catch (XMLStreamException | IOException e) {
				e.printStackTrace();
				System.out.println("Fail loading SBML File!");
				System.exit(1);
			}
		} 
		else {
			System.out.println("wrong function call. use: MainClass.java PATHOFSBML\n or to ignore: MainClass.java PATHOFSBML -n PATHFORERROR");
			System.exit(1);
		}

		// call FindFluxModules
		try {
			FindFluxModules ff = new FindFluxModules(load,ignore);
			ff.findFlux();
		} catch (FileNotFoundException | LpSolveException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		// more methods: e.g. pdf,graphviz,...

		System.out.println("\nfinish");

	}
}
