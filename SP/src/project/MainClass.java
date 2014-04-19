package project;


import java.io.FileNotFoundException;
import java.io.IOException;

import javax.xml.stream.XMLStreamException;
import lpsolve.LpSolveException;



public class MainClass {
	
	public static void main(String args[]){
		//load SBML etc.
			SBMLLoad load = new SBMLLoad();
		
		
			if(args.length>0){		
				try {
					System.out.println("Loading SBML file");
					load.loadSBML(args[0]);
				} catch (XMLStreamException | IOException e) {
					e.printStackTrace();
					System.out.println("Fail loading SBML File!");
				}
			}
			else{
				try {
					//0.7mb H_pylori_iIT341.xml 
					//1.7mb S_cerevisiae_iND750.xml
					//3.2mb Y_pestis_co92.xml
					System.out.println("Loading SBML file");
					load.loadSBML("/home/guru/Downloads/H_pylori_iIT341.xml");
				} catch (XMLStreamException | IOException e) {
					e.printStackTrace();
					System.out.println("Fail loading SBML File!");
				}
			}
			
			
		//call FindFluxModules
			try {
				FindFluxModules ff = new FindFluxModules(load);
				ff.findFlux();
			} catch (FileNotFoundException | LpSolveException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			
			
			
			
		//more methods: e.g. pdf,graphviz,...
	
		
		
		System.out.println("finish");
	
	}
}





