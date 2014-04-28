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
import scpsolver.constraints.LinearConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;
import scpsolver.util.SparseMatrix;

public class TestComputeModules {

	SBMLLoad load;
	
	@Before
	public void setUp() throws XMLStreamException, IOException {
		load = new SBMLLoad();
		load.loadSBML("/home/guru/Downloads/S_aureus_iSB619.xml");
	}
	
	@Test
	public void testMinModules() throws LpSolveException{
		boolean[] var = {false,false,false,false,false,false,true,true,true,true,true,true,true,true,true,false,true,false,false,true,false,false,false,false,true,true,true,true,true,true,true,true,true,true,true,true,true,true,false,false,false,true,false,false,false,false,false,false,false,false,false,true,false,true,true,false,false,false,false,true,true,false,true,false,true,true,true,false,true,true,false,false,true,true,false,true,true,false,false,false,false,false,true,true,false,false,true,true,false,true,true,true,false,false,false,true,true,false,false,true,false,true,false,false,false,true,false,false,true,true,true,true,false,false,true,true,false,false,false,false,false,true,false,false,false,true,false,true,false,false,false,false,true,true,true,false,true,true,false,true,false,false,false,true,false,false,true,false,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,true,true,true,false,false,false,false,false,false,false,true,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,true,false,false,true,false,false,false,false,false,false,true,false,true,true,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,false,true,false,false,false,true,false,false,false,true,true,true,false,false,true,true,true,false,false,true,true,false,false,true,false,false,false,false,true,false,false,true,true,true,true,true,false,false,false,false,true,false,true,true,false,true,true,false,false,true,true,true,false,false,true,true,false,false,false,true,false,true,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,true,false,false,true,false,false,false,false,false,false,true,true,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,true,true,false,true,true,false,false,true,true,true,false,true,false,false,false,false,false,true,true,true,false,false,true,false,false,false,false,false,true,true,true,true,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,true,true,true,false,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,true,true,false,true,false,false,true,true,false,false,true,false,false,false,true,true,true,true,false,true,true,true,true,true,false,false,true,false,false,false,false,false,true,true,false,false,false,false,false,false,false,true,true,true,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,false,false,false,true,true,false,true,false,true,false,false,false,true,false,false,false,false,false,false,false,false,false,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,true,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,false,true,false,false,true,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,true,false,false,true,true,true,true,false,false,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
		
		
		FindFluxModules ff = new FindFluxModules(load);
		SparseMatrix rctMetArr = ff.matrixBuild();
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		
		ff.computeMinModules(var, rctMetArr, solver);
		
	}

	
	@Test
	public void testMinMax(){
		
	}
	

//optimize
	@Test
	public void FindFluxModules_optimizeParameterControl() throws LpSolveException{
		FindFluxModules ff = new FindFluxModules(load);
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		
		SparseMatrix rctMetArr = ff.matrixBuild();
		double[] solved = ff.optimize(solver,rctMetArr);
		assertNotNull(solved);
		double opt = ff.getBiomassOptValue();
		assertTrue(opt>0);
	}
	
	@Test	//with 
	public void testForMinModules(){
		SparseMatrix test = new SparseMatrix(6,6);
		test.set(0,0,-1);test.set(0,1,1);test.set(1,0,1);test.set(1,1,-1);test.set(2,2,-1);
		test.set(2,3,1);test.set(3,3,-1);test.set(3,4,1);test.set(4,2,1);test.set(4,4,-1);
		boolean[] isVariable = {true,true,true,true,true,false};
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		SBMLLoad load2 = new SBMLLoad();
		load2.numS = 5;
		load2.numR = 5;
		
		FindFluxModules ff = new FindFluxModules(load2);
		ff.computeMinModules(isVariable, test, solver);
		ff.printModules();
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
