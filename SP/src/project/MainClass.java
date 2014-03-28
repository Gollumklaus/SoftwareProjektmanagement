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
					load.loadSBML(args[0]);
				} catch (XMLStreamException | IOException e) {
					e.printStackTrace();
					System.out.println("Fail loading SBML File!");
				}
			}
			else{
				try {
					load.loadSBML("/home/guru/Downloads/S_cerevisiae_iND750.xml");
					System.out.println("");
				} catch (XMLStreamException | IOException e) {
					e.printStackTrace();
					System.out.println("Fail loading SBML File!");
				}
			}
			
			
		//call FindFluxModules
			try {
				FindFluxModules.findFlux(load);
			} catch (FileNotFoundException | LpSolveException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			
			
			
			
		//more methods: e.g. pdf,graphviz,...
	
		
		
		System.out.println("finish");
	
	}
}





