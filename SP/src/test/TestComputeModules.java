package test;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.Arrays;

import javax.xml.stream.XMLStreamException;

import lpsolve.LpSolveException;

import org.junit.Before;
import org.junit.Test;
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
		load.loadSBML("/home/guru/Downloads/S_aureus_iSB619.xml");
	}
	
	//test minModules with pure true vector
	@Test
	public void trueVectorForMinModules() throws LpSolveException{
	
		boolean ignore = true;
		FindFluxModules ff = new FindFluxModules(load,ignore);
		SparseMatrix rctMetArr = ff.matrixBuild();
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		boolean[] var = new boolean[load.getNumR()];
		Arrays.fill(var,Boolean.TRUE);
		SparseMatrix spma = ff.computeMinModules(var, rctMetArr, solver);
		ff.dfs(var,spma);
		ff.printModules();
		
	}
	
	//test minModules with computed (**OUT_DATED**) variability vector of h_pylori 
	@Test
	public void testMinModules() throws LpSolveException{
		boolean[] var = {false,false,false,false,false,false,true,true,true,true,true,true,true,true,true,false,true,false,false,true,false,false,false,false,true,true,true,true,true,true,true,true,true,true,true,true,true,true,false,false,false,true,false,false,false,false,false,false,false,false,false,true,false,true,true,false,false,false,false,true,true,false,true,false,true,true,true,false,true,true,false,false,true,true,false,true,true,false,false,false,false,false,true,true,false,false,true,true,false,true,true,true,false,false,false,true,true,false,false,true,false,true,false,false,false,true,false,false,true,true,true,true,false,false,true,true,false,false,false,false,false,true,false,false,false,true,false,true,false,false,false,false,true,true,true,false,true,true,false,true,false,false,false,true,false,false,true,false,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,true,true,true,false,false,false,false,false,false,false,true,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,true,false,false,true,false,false,false,false,false,false,true,false,true,true,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,false,true,false,false,false,true,false,false,false,true,true,true,false,false,true,true,true,false,false,true,true,false,false,true,false,false,false,false,true,false,false,true,true,true,true,true,false,false,false,false,true,false,true,true,false,true,true,false,false,true,true,true,false,false,true,true,false,false,false,true,false,true,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,true,false,false,true,false,false,false,false,false,false,true,true,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,true,true,false,true,true,false,false,true,true,true,false,true,false,false,false,false,false,true,true,true,false,false,true,false,false,false,false,false,true,true,true,true,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,true,true,true,false,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,true,true,false,true,false,false,true,true,false,false,true,false,false,false,true,true,true,true,false,true,true,true,true,true,false,false,true,false,false,false,false,false,true,true,false,false,false,false,false,false,false,true,true,true,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,false,false,false,true,true,false,true,false,true,false,false,false,true,false,false,false,false,false,false,false,false,false,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,true,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,false,true,false,false,true,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,true,false,false,true,true,true,true,false,false,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
		boolean ignore = false;
		FindFluxModules ff = new FindFluxModules(load,ignore);
		SparseMatrix rctMetArr = ff.matrixBuild();
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		
		SparseMatrix spma = ff.computeMinModules(var, rctMetArr, solver);
		ff.dfs(var, spma);
		ff.printModules();
	}
	

	//test if optimize gets the correct parameter
	@Test
	public void FindFluxModules_optimizeParameterControl() throws LpSolveException{
		boolean ignore = false;
		FindFluxModules ff = new FindFluxModules(load,ignore);
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		
		SparseMatrix rctMetArr = ff.matrixBuild();
		double[] solved = ff.optimize(solver,rctMetArr);
		assertNotNull(solved);
		double opt = ff.getBiomassOptValue();
		assertTrue(opt>0);
	}
	
	//small example for minModules
	@Test	 
	public void testForMinModules(){
		SparseMatrix test = new SparseMatrix(6,6);
		test.set(0,0,-1);test.set(0,1,1);test.set(1,0,1);test.set(1,1,-1);test.set(2,2,-1);
		test.set(2,3,1);test.set(3,3,-1);test.set(3,4,1);test.set(4,2,1);test.set(4,4,-1);
		boolean[] isVariable = {true,true,true,true,true,false};
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		SBMLLoad load2 = new SBMLLoad();
		load2.numS = 5;
		load2.numR = 5;
		boolean ignore = true;
		FindFluxModules ff = new FindFluxModules(load2,ignore);
		ff.computeMinModules(isVariable, test, solver);
		ff.printModules();
	}
	
	
	//minMax
	@Test
	public void testFindingMinimalModules() throws LpSolveException{
		boolean ignore = true;
		FindFluxModules ff = new FindFluxModules(load,ignore);
		SparseMatrix rctMetArr = ff.matrixBuild();
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		
		double[] solved = ff.optimize(solver,rctMetArr);
		boolean[] isConstant = new boolean[load.getNumR()];
		ff.computeMinModules(isConstant, rctMetArr, solver);
		
	}
	
	//computes the optimal biomass and print value
	@Test
	public void getBiomass(){
		int c =  load.getBiomassOptValuePos();
		System.out.println("position Biomass: " +c);	
		Reaction r = load.getModel().getReaction(c);
		String id = r.getId();
		System.out.println("idBiomass: " + id);
		System.out.println();
	}
}
