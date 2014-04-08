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
					System.out.println("loading SBML file");
					load.loadSBML(args[0]);
				} catch (XMLStreamException | IOException e) {
					e.printStackTrace();
					System.out.println("Fail loading SBML File!");
				}
			}
			else{
				try {
					System.out.println("loading SBML file");
					load.loadSBML("/home/guru/Downloads/S_aureus_iSB619.xml");
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





