package project;

import java.io.File;
import java.io.IOException;

import javax.xml.stream.XMLStreamException;


import org.junit.runners.Parameterized;
import org.sbml.jsbml.*;



public class Test{
	
	

	public static void main(String args[]) throws XMLStreamException, IOException{
		
		SBMLDocument d = SBMLReader.read(new File("/home/guru/Downloads/PLUSPARAMETERM_barkeri_iAF692.xml"));
		Model m = d.getModel();
		int numR = m.getNumReactions();
		int numS = m.getNumSpecies();
		
		
		//create reference arrays for reactions and metabolites
		String[] rct = new String[m.getNumReactions()];
		String[] met = new String[m.getNumSpecies()];
		
		for(int i=0;i<m.getNumReactions();i++){
			Reaction r = m.getReaction(i);
			String id = r.getId();
			rct[i] = id;
		}
		
		for(int i=0;i<m.getNumSpecies();i++){
			Species s = m.getSpecies(i);
			String id = s.getId();
			met[i] = id;
		}
		
		//array with reactions and metabolites 
		int[][] rctMetArr = new int[m.getNumReactions()][m.getNumSpecies()];
		
			
		Reaction r = m.getReaction(0);
		SpeciesReference s = r.getReactant(0);
		Species in = s.getSpeciesInstance();
		String id = in.getId();
		System.out.println(s);
		System.out.println(met[0]);
		
		for(int k=0;k<numS;k++){
			if(met[k]== id ){
				System.out.print(id + " : ");
				System.out.println(met[k]);
			}
		}
		
		
		for(int i=0;i<numR;i++){
			for(int k=0;k<numS;k++){
		//			System.out.print(rctMetArr[i][k]);
			}
//			System.out.println();
		}

		
		/*
		for(int i=0;i<m.getNumReactions();i++){
			Reaction r = m.getReaction(i);
			
			
			Parameter para = m.getParameter(2);
			double ken = para.getValue();
			System.out.println("Parameter = " + ken);
					
			int rea = r.getNumReactants();
			int pro = r.getNumProducts();
			
			String id = r.getId();
			System.out.println(id + ":" );
						
			
			for(int k=0;k<rea;k++){
				double stoch = r.getReactant(k).getStoichiometry();
				SpeciesReference react = r.getReactant(k);
				System.out.println("\t" + react + " | " + stoch);
			}
		}
		*/
		
		
			
		
		
						//	new JSBMLvisualizer(SBMLReader.read(new File("/home/guru/Downloads/H_pylori_iIT341.xml")));
						//	SBMLReader.read(new File("/home/guru/Arbeitsfläche/blabla.xml"));
						//	new JSBMLvisualizer(SBMLReader.read(new File("/home/guru/Arbeitsfläche/blabla.xml")));
		
		
	System.out.println("0");	

	}	


	
}












 