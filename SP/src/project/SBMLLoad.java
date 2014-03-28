package project;

import java.io.File;
import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.LocalParameter;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.Reaction;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLReader;
import org.sbml.jsbml.Species;



public class SBMLLoad {


		//these variables get set by loadSBML
			//number of reactions
			private int numR;
			//number of metabolites (species)
			private int numS;
			//model
			private Model model;
			//id's of all reactions
			private  String[] rct;
			//id's of all metabolites
			private  String[] met;
			//sbml file
			private  File file;
			//objectiveFunction
			private double[] objectiveFunction;
			//upper bound and lower bound
			private double[] lowerBound;
			private double[] upperBound;
			private int biomassOptValuePos;
	//getter

			public int getNumR() {
						return numR;
					}
			public int getNumS() {
				return numS;
			}
			public Model getModel() {
				return model;
			}
			public String[] getRct(){ 
				return rct;
			}
			public String[] getMet() {
				return met;
			}
			public File getFile() {
				return file;
			}
			public double[] getObjectiveFunction(){
				return objectiveFunction;
			}
			public double[] getLowerBound(){
				return lowerBound;
			}
			public double[] getUpperBound(){
				return upperBound;
			}
			public int getBiomassOptValuePos() {
				return biomassOptValuePos;
			}
			
			
	//load model
	public void loadSBML(String path) throws XMLStreamException, IOException{
		file = new File(path);
		System.out.println("loading SBML file");
		SBMLDocument d = SBMLReader.read(file);
		model = d.getModel();
		
		numR = model.getNumReactions();
		numS = model.getNumSpecies();
		
		rct = new String[model.getNumReactions()];
		met = new String[model.getNumSpecies()];
		
		for(int i=0;i<model.getNumReactions();i++){
			Reaction r = model.getReaction(i);
			String id = r.getId();
			rct[i] = id;
		}
		
		for(int i=0;i<model.getNumSpecies();i++){
			Species s = model.getSpecies(i);
			String id = s.getId();
			met[i] = id;
		}
		
		
		//set objective function, set lower and upper bound
		objectiveFunction = new double[numR];
		lowerBound = new double[numR];
		upperBound = new double[numR];
int count =0;		
		for(int i=0;i<numR;i++){
			Reaction r = model.getReaction(i);
			KineticLaw k = r.getKineticLaw();
			LocalParameter p1 = k.getLocalParameter("OBJECTIVE_COEFFICIENT");
			LocalParameter p2 = k.getLocalParameter("LOWER_BOUND");
			LocalParameter p3 = k.getLocalParameter("UPPER_BOUND");
			
			double v = p1.getValue();
			
			lowerBound[i] = p2.getValue();
			upperBound[i] = p3.getValue();
			

			if(v==0.0){
				objectiveFunction[i] = 0;
			}
			else{
				objectiveFunction[i] = 1;
				biomassOptValuePos = i;
				System.out.println("para: "+v);
				System.out.println("rct: "+r.getId());
			}
			
		}
		 
		 
		//test
		if(model != null && numR !=0 && numS !=0 && (rct.length>0) && (met.length>0)){	}
		else{System.out.println("Error load data in loadSBML."); System.exit(1);}
			
	}

}
