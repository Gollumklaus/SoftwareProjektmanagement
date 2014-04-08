package test;

import static org.junit.Assert.*;

import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import lpsolve.LpSolveException;

import org.junit.Before;
import org.junit.Test;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.Reaction;

import project.FindFluxModules;
import project.SBMLLoad;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.util.SparseMatrix;

public class TestComputeModules {

	SBMLLoad load;
	
	@Before
	public void setUp() throws XMLStreamException, IOException {
		load = new SBMLLoad();
		load.loadSBML("/home/guru/Downloads/S_cerevisiae_iND750.xml");
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

//rctmetmatrix
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
	
//optimize
	public void FindFluxModules_optimizeParameterControl() throws LpSolveException{
		FindFluxModules ff = new FindFluxModules(load);
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		
		SparseMatrix rctMetArr = ff.matrixBuild();
		double[] solved = ff.optimize(solver,rctMetArr);
		assertNotNull(solved);
		double opt = ff.getBiomassOptValue();
		assertTrue(opt>0);
	}
	
//minMax
	@Test
	public void testFindingMinimalModules() throws LpSolveException{
		FindFluxModules ff = new FindFluxModules(load);
		SparseMatrix rctMetArr = ff.matrixBuild();
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		double[] solved = ff.optimize(solver,rctMetArr);
		boolean[] isConstant = new boolean[load.getNumR()];
		ff.computeMinModules(isConstant, rctMetArr, solver);
		
	}
	
	
	
	@Test
	public void testComputeModules() {
		FindFluxModules ff = new FindFluxModules(load);
		ff.getModules();
	}
	
	@Test
	public void feasible(){
		
	}
	
	@Test
	public void getBiomass(){
		int c =  load.getBiomassOptValuePos();
		System.out.println("position Biomass: " +c);	
		Reaction r = load.getModel().getReaction(c);
		String id = r.getId();
		System.out.println("idBiomass: " + id);
		System.out.println();
	}
	
	
	
	@Test //optimal biomass is biggest value in matrix?
	public void biomassOptValue() throws LpSolveException{
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		FindFluxModules ff = new FindFluxModules(load);
		double[] v = ff.optimize(solver,ff.matrixBuild());
		
		double max = 0;

		for(int x=0;x<load.getNumR();x++){
			if(v[x]>max){
				max = x;
				System.out.println("error");
			}
			if(v[x]>500){
				System.out.println("greatValues: "+v[x] + "\tposition: " + x);
			}
			
		}
		double biomass = v[load.getBiomassOptValuePos()];
		
		assertTrue(Math.abs(biomass-max)<1.0);
		System.out.println("max: "+max +"     biomassOpt: " + biomass );
	
	
	}
	
	@Test
	public void buildVectoren(){
		assertTrue(true);
	}
}
