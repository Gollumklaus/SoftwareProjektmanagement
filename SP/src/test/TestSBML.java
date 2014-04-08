package test;

import static org.junit.Assert.*;

import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import org.junit.Test;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.Reaction;

import project.FindFluxModules;
import project.SBMLLoad;
import scpsolver.util.SparseMatrix;

public class TestSBML {

	SBMLLoad load;
	
	@Test
	public void testLoadSBML() throws XMLStreamException, IOException {
		load = new SBMLLoad();
		load.loadSBML("SP/data/Ec_iJR904.xml");
	}
	
	@Test
	public void testBuildMatrix() throws XMLStreamException, IOException {
		//SBMLLoad load = new SBMLLoad();
		//load.loadSBML("SP/data/Ec_iJR904.xml");
		
		if(load == null) testLoadSBML();
		
		int numR = load.getNumR();
		int numS = load.getNumS();
		Model m = load.getModel();
		String[] met = load.getMet();
		double[] objectiveFunction = load.getObjectiveFunction();
		int biomassOptValuePos = load.getBiomassOptValuePos();
		double[] low = load.getLowerBound();
		double[] upp = load.getUpperBound();
		
		//build matrix
		FindFluxModules ff = new FindFluxModules(load);
		SparseMatrix rctMetArr = ff.matrixBuild();
	}
	
/*	
	
	@Test
	public void allMetabolitesInMatrix(){
	double wus=0;
	int error=0;
	
	for(int i=0;i<numR;i++){
		Reaction r = m.getReaction(i);
		
		try{
			int a = r.getNumReactants();
				wus += r.getNumReactants();
		}
		catch(NullPointerException tt){
			error++;
		}
		
		try{
			int b = r.getNumProducts();
			
			wus += r.getNumProducts();
		}
		catch(NullPointerException tt){
			error++;
		}
	
	System.out.println(wus);
	System.out.println(error);
	}
	
*/
	
	
}
