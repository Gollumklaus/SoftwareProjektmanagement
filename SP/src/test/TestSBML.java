package test;

import static org.junit.Assert.*;

import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import org.junit.Before;
import org.junit.Test;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.Reaction;

import project.FindFluxModules;
import project.SBMLLoad;
import scpsolver.util.SparseMatrix;

public class TestSBML {

	SBMLLoad load;
	
	@Before
	public void setUp() throws XMLStreamException, IOException {
		load = new SBMLLoad();
		load.loadSBML("data/Ec_iJR904.xml");
	}
	
	@Test
	public void isSetUp() {
		assertTrue(load != null);
	}

	
	//sbmlload
	 @Test
	 public void SBMLSetParameter(){
		 assertNotNull(load.getBiomassOptValuePos());
		 assertTrue(load.getNumR() != 0);
		 assertTrue(load.getNumS() != 0);
			assertNotNull(load.getModel());
			assertNotNull(load.getObjectiveFunction());
			assertTrue(load.getRct().length==load.getNumR());
			assertTrue(load.getMet().length==load.getNumS());
			assertTrue(load.getLowerBound().length==load.getNumR());
			assertTrue(load.getUpperBound().length==load.getNumR());
			int count = 0;
			for(int i=0;i<load.getNumR();i++){
				if(load.getObjectiveFunction()[i] == 1){
					count++;
				}
			}
			assertTrue(count == 1);
	 }
	 
	 //buildMatrix
	 @Test
	 public void correctRctMetMatrix(){	//compares parameter in matrix with real number of reactants/products
		
		FindFluxModules ff = new FindFluxModules(load);
		SparseMatrix sm = ff.matrixBuild();
		int count = 0;
		for(int i=0;i<load.getNumR();i++){
			for(int k=0;k<load.getNumS();k++){
				if(sm.get(i,k) != 0){
					count++;
				}
			}
		}
		Model m = load.getModel();
		int real = 0;
		for(int i=0;i<load.getNumR();i++){
			Reaction r = m.getReaction(i);
				real += r.getNumReactants();
				real+= r.getNumProducts();
		}
		assertTrue(count == real);
	}
		
}
