package be.iminds.testcloudspeller;


import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({TestBLSCalculator.class, TestDegSuffixTree.class, TestExactDiscovery.class
	, TestByteRepresentations.class})
public class AllTests {

}