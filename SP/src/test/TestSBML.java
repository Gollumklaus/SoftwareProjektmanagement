package test;

import static org.junit.Assert.*;

import org.junit.Test;
import org.sbml.jsbml.Reaction;

public class TestSBML {

	@Test
	public void testLoadSBML() {
		fail("Not yet implemented");
	}
	
	@Test
	public void testLoadSBMLcorrect() {
		// does nothing
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
