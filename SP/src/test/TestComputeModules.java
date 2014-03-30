package test;

import static org.junit.Assert.*;

import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import org.junit.Before;
import org.junit.Test;

import project.SBMLLoad;

public class TestComputeModules {

	SBMLLoad load;
	
	@Before
	public void setUp() throws XMLStreamException, IOException {
		load = new SBMLLoad();
		load.loadSBML("SP/data/Ec_iJR904.xml");
	}
	
	@Test
	public void isSetUp() {
		assertTrue(load != null);
	}
	
	@Test
	public void testFBA() {
		fail("Not yet implemented");
	}
	
	@Test
	public void testFVA() {
		fail("Not yet implemented");
	}
	
	@Test
	public void testComputeModules() {
		fail("Not yet implemented");
	}

}
