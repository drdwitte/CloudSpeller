package be.iminds.testcloudspeller;

import static org.junit.Assert.*;
import org.junit.Test;

/* Tutorial http://www.vogella.com/articles/JUnit/article.html
@Test
public void testMultiply() {

   // MyClass is tested
   MyClass tester = new MyClass();
   
   // Check if multiply(10,5) returns 50
   assertEquals("10 x 5 must be 50", 50, tester.multiply(10, 5));
 }
*/

/*
Annotations: @Test, @Before,@BeforeClass,@Ignore,@Test(expected=Exception.class),@Test(timeout=100) 
Test Methods:
	fail, assertTrue, assertEquals, assertNull, assertNotNull, assertSame(refcheck), assertNotSame
*/


public class ExampleJUnit {


	/*
	@Test
	public void testFail(){
		fail("Example of a fail message");
	}
	*/
	
	@Test
	public void testAssertTrue(){
		assertTrue(true);
		//assertTrue(false);
	}
	
	/*@Test
	public void testAssertEquals(){
		assertEquals(0.005,0.0053,3);
	}*/

/*
	@BeforeClass
  public static void testSetup() {
  }

  @AfterClass
  public static void testCleanup() {
    // Teardown for data used by the unit tests
  }

  @Test(expected = IllegalArgumentException.class)
  public void testExceptionIsThrown() {
    MyClass tester = new MyClass();
    tester.multiply(1000, 5);
  }

  @Test
  public void testMultiply() {
    MyClass tester = new MyClass();
    assertEquals("10 x 5 must be 50", 50, tester.multiply(10, 5));
  } 
*/
}
