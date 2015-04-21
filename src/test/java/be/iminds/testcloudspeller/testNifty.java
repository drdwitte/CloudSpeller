package be.iminds.testcloudspeller;

import be.iminds.cloudspeller.test.NiftyPigTest;
import be.iminds.cloudspeller.test.reporter.StringReporter;
import be.iminds.cloudspeller.test.result.DataSetReport;
import be.iminds.cloudspeller.test.validator.DataSetValidator;
import org.junit.Assert;
import org.junit.Test;

import static be.iminds.cloudspeller.test.validator.DataSetValidator.dataset;
import static be.iminds.cloudspeller.test.validator.FieldValidator.isNull;
import static be.iminds.cloudspeller.test.validator.FieldValidator.string;
import static be.iminds.cloudspeller.test.validator.TupleValidator.tuple;

/**
 * Created by dieter on 13.11.14.
 */
public class testNifty {


       private static final String PIG_SCRIPT = "prep_arbor_account_info.pig";

        @Test
        /**
         * Test 1 record
         */
        public void test_prep_arbor_account_info_positive_1() throws Exception{

            String[] raw_data = {
                    "21758;495987;herman.vanherp@telenet.be;Dhr;Van Herp;Herman;" +
                            "null;Landbouwstraat 39----;Landbouwstraat 39----;Mechelen;null;" +
                            "2800;Mechelen;2800;null;25111998;M01;" +
                            "3;0;3;e-UTB;1;" +
                            "06062014 10:29:36;null;46968;177392;171105"
            };

            NiftyPigTest test = new NiftyPigTest(PIG_SCRIPT);
            test.input("raw_data",raw_data,NiftyPigTest.STORAGE_PIG_CSV);

            test.execute();

            DataSetReport report = test.validate(dataset("cleaned_data").size(1).mode(DataSetValidator.ValidationMode.ByOrder)
                            .add(tuple()
                                            .field(string("21758"))
                                            .field(string("495987"))
                                            .field(string("herman.vanherp@telenet.be"))
                                            .field(string("1998-11-25T12:00:00.000Z"))
                                            .field(isNull())
                                            .field(string("e-UTB"))
                                            .field(string("46968"))
                                            .field(string("177392"))
                                            .field(string("171105"))
                                            .field(string("2014-06-06T10:29:36.000Z"))
                                            .field(string("arbor"))
                                            .field(string("account"))
                            )
            );
            System.out.println(StringReporter.format(report));

            // -- check the report was valid
            Assert.assertTrue(report.isValid());
        }

        @Test
        /**
         * Test 3 records
         */
        public void test_prep_arbor_account_info_positive_2() throws Exception{

            String[] raw_data = {
                    "21757;495987;herman.vanherp@telenet.be;Dhr;Van Herp;Herman;null;Landbouwstraat 39----;Landbouwstraat 39----;Mechelen;null;2800;Mechelen;2800;null;25111998;M01;3;0;3;e-UTB;1;06062014 10:29:36;null;46968;177392;171105",
                    "21758;495987;herman.vanherp@telenet.be;Dhr;Van Herp;Herman;null;Landbouwstraat 39----;Landbouwstraat 39----;Mechelen;null;2800;Mechelen;2800;null;25111998;M01;3;0;3;e-UTB;1;06062014 10:29:36;null;46968;177392;171105",
                    "21759;495987;herman.vanherp@telenet.be;Dhr;Van Herp;Herman;null;Landbouwstraat 39----;Landbouwstraat 39----;Mechelen;null;2800;Mechelen;2800;null;25111998;M01;3;0;3;e-UTB;1;06062014 10:29:36;null;46968;177392;171105"
            };

            NiftyPigTest test = new NiftyPigTest(PIG_SCRIPT);
            test.input("raw_data",raw_data,NiftyPigTest.STORAGE_PIG_CSV);

            test.execute();

            DataSetReport report = test.validate(dataset("cleaned_data").size(3).mode(DataSetValidator.ValidationMode.ByOrder)
                            .add(tuple()
                                            .field(string("21757"))
                                            .field(string("495987"))
                                            .field(string("herman.vanherp@telenet.be"))
                                            .field(string("1998-11-25T12:00:00.000Z"))
                                            .field(isNull())
                                            .field(string("e-UTB"))
                                            .field(string("46968"))
                                            .field(string("177392"))
                                            .field(string("171105"))
                                            .field(string("2014-06-06T10:29:36.000Z"))
                                            .field(string("arbor"))
                                            .field(string("account"))
                            ).add(tuple()
                                            .field(string("21758"))
                                            .field(string("495987"))
                                            .field(string("herman.vanherp@telenet.be"))
                                            .field(string("1998-11-25T12:00:00.000Z"))
                                            .field(isNull())
                                            .field(string("e-UTB"))
                                            .field(string("46968"))
                                            .field(string("177392"))
                                            .field(string("171105"))
                                            .field(string("2014-06-06T10:29:36.000Z"))
                                            .field(string("arbor"))
                                            .field(string("account"))
                            ).add(tuple()
                                            .field(string("21759"))
                                            .field(string("495987"))
                                            .field(string("herman.vanherp@telenet.be"))
                                            .field(string("1998-11-25T12:00:00.000Z"))
                                            .field(isNull())
                                            .field(string("e-UTB"))
                                            .field(string("46968"))
                                            .field(string("177392"))
                                            .field(string("171105"))
                                            .field(string("2014-06-06T10:29:36.000Z"))
                                            .field(string("arbor"))
                                            .field(string("account"))
                            )
            );
            System.out.println(StringReporter.format(report));

            // -- check the report was valid
            Assert.assertTrue(report.isValid());
        }

        @Test
        /**
         * Test Null input
         */
        public void test_prep_arbor_account_info_negative_1() throws Exception{

            String[] raw_data = {
                    "null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null;null"
            };

            NiftyPigTest test = new NiftyPigTest(PIG_SCRIPT);
            test.input("raw_data",raw_data,NiftyPigTest.STORAGE_PIG_CSV);

            test.execute();

            DataSetReport report = test.validate(dataset("cleaned_data").size(1).mode(DataSetValidator.ValidationMode.ByOrder)
                            .add(tuple()
                                            .field(string("null"))
                                            .field(string("null"))
                                            .field(string("null"))
                                            .field(isNull())
                                            .field(isNull())
                                            .field(string("null"))
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(string("arbor"))
                                            .field(string("account"))
                            )
            );
            System.out.println(StringReporter.format(report));

            // -- check the report was valid
            Assert.assertTrue(report.isValid());
        }

        @Test
        /**
         * Test blank input
         */
        public void test_prep_arbor_account_info_negative_2() throws Exception{
            String[] raw_data = {
                    ";;;;;;;;;;;;;;;;;;;;;;;;;;"
            };

            NiftyPigTest test = new NiftyPigTest(PIG_SCRIPT);
            test.input("raw_data",raw_data,NiftyPigTest.STORAGE_PIG_CSV);

            test.execute();

            DataSetReport report = test.validate(dataset("cleaned_data").size(1).mode(DataSetValidator.ValidationMode.ByOrder)
                            .add(tuple()
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(string("arbor"))
                                            .field(string("account"))
                            )
            );
            System.out.println(StringReporter.format(report));

            // -- check the report was valid
            Assert.assertTrue(report.isValid());
        }
        @Test
        /**
         * Test input with big numbers
         */
        public void test_prep_arbor_account_info_negative_3() throws Exception{

            String[] raw_data = {
                    ";;;;;;;;;;;;;;;;;;;;;;;;1000000000000;1000000000000;1000000000000"
            };

            NiftyPigTest test = new NiftyPigTest(PIG_SCRIPT);
            test.input("raw_data",raw_data,NiftyPigTest.STORAGE_PIG_CSV);

            test.execute();

            DataSetReport report = test.validate(dataset("cleaned_data").size(1).mode(DataSetValidator.ValidationMode.ByOrder)
                            .add(tuple()
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(isNull())
                                            .field(string("arbor"))
                                            .field(string("account"))
                            )
            );
            System.out.println(StringReporter.format(report));

            // -- check the report was valid
            Assert.assertTrue(report.isValid());
        }


}

