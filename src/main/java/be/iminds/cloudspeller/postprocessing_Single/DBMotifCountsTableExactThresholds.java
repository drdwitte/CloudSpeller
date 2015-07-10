package be.iminds.cloudspeller.postprocessing_Single;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.output.SimultaneousOccurrenceConfidenceAndBLSFiltering;
import be.iminds.cloudspeller.output.TripleThresholdFilterExactBLST;
import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by ddewitte on 08.12.14.
 */
public class DBMotifCountsTableExactThresholds extends Configured implements Tool {

    private static String input;
    private static String output;

    public static class MotifMapper extends Mapper<LongWritable,Text,Text,IntWritable> {

        private static final int [] t = {15,50,60,70,90,95};

        private static final int [] cThr = { 50,60,70,80,90,95};
        private static final int [] blsThr = {15,70,95};
        private static final int [] Fthr = {1,5,10,15,20,50,500};
        private static final int [] motifLength = {6,7,8,9,10,11,12};
        private static final int [] degeneracies = {1,2,4,8,16,32,64};

        private static String [] strKeys;
        private static int [] F;
        private static double [] p;


        private static int numRestrictions;
        private static Map<String,Integer> countMap = new HashMap<String,Integer>();
        private static Text mapOutputKey = new Text("");
        private static IntWritable mapOutputValue = new IntWritable(1);

        private static ConfidenceGraphRestrictions [] restrictions;

        @Override
        protected void setup(Context context) throws IOException,
                InterruptedException {
            BLS.initializeBLSConstants(t);
            FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());
            F = new int[FreqVec.getNumberOfIntervals()];
            p = new double[FreqVec.getNumberOfIntervals()];

            numRestrictions = cThr.length*blsThr.length*Fthr.length;


            restrictions = new ConfidenceGraphRestrictions[numRestrictions];


            strKeys = new String[numRestrictions];

            for (int i=0; i<restrictions.length; i++){
                int b = blsThr[i%blsThr.length];             //012012012012012012
                int c = cThr[(i/blsThr.length)%cThr.length]; //000111222000111222
                int f = Fthr[i/(blsThr.length*cThr.length)]; //0000000001111111115555555 0000011115555

                restrictions[i] = new TripleThresholdFilterExactBLST(c, f, b);

                strKeys[i] = "_"+b+"_"+c+"_"+f;

                for (int k : motifLength) {
                    for (int d: degeneracies) {
                        countMap.put(k + "_" + d + strKeys[i], 0);
                    }
                }
            }



        }

        @Override
        protected void map(LongWritable key, Text value, Context context)
                throws IOException, InterruptedException {


            String strValue = value.toString();
            if (strValue == null || strValue.length()==0)
                return;






            GeneralToolbox.parseConfidenceGraphValues(strValue,F,p);

            for (int i=0; i<restrictions.length; i++){

                if (restrictions[i].checkRestrictions(F,p)){
                    String ckey = getMotifLength(strValue)+"_"+getMotifDegeneracy(strValue)+strKeys[i];
                    countMap.put(ckey,countMap.get(ckey)+1);
                }
            }






        }



        private static void pr(String s){
            System.err.println(s);
        }

        private int getMotifLength(String confGraphRecord){

            return confGraphRecord.indexOf('\t');

        }

        @Override
        protected void cleanup(Context context) throws IOException, InterruptedException {
            //flush countmap

            for (Map.Entry<String,Integer> e : countMap.entrySet()){


                mapOutputKey.set(e.getKey());
                mapOutputValue.set(e.getValue());
                context.write(mapOutputKey,mapOutputValue);
            }

            countMap.clear();
        }

        private static final String baseChars = "ACGT";

        private int getMotifDegeneracy(String confGraphRecord) {

            int twofolds = 0;
            int Ns = 0;
            int length = getMotifLength(confGraphRecord);
            String group = confGraphRecord.substring(0,length);

            for (char c : group.toCharArray()){
                if (c == 'N'){
                    Ns++;
                } else if (baseChars.indexOf(c)>=0){

                } else {
                    twofolds++;
                }
            }

            int deg = (1 << twofolds) * (1 << 2*Ns);
            return deg;
        }


    }
;

    static class MotifReducer extends Reducer<Text,IntWritable,Text,IntWritable> {


        private static IntWritable outputValue = new IntWritable();

        @Override
        protected void reduce(Text key, Iterable<IntWritable> values,
                              Context context)
                throws IOException, InterruptedException {

            int sum = 0;
            for (IntWritable i : values){
                sum+=i.get();
            }

            outputValue.set(sum);
            context.write(key,outputValue);

            /*if (sum!=0){
                System.err.println(key.toString()+" : "+sum);
            }*/
        }

        private static void pr(String s){
            System.err.println(s);
        }

    }



        @Override
    public int run(String[] args) throws Exception {


        //NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
        Job job = Job.getInstance(getConf(), "MotifCountsTable");

        job.setJarByClass(DBMotifCountsTableExactThresholds.class);

        //SETTING MAPREDUCE CLASSES

        // IO
        job.setInputFormatClass(TextInputFormat.class);
        job.setOutputFormatClass(TextOutputFormat.class);

        // (K3,V3)
		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(IntWritable.class);

        // set hadoop methods
        job.setMapperClass(MotifMapper.class);
        job.setReducerClass(MotifReducer.class);

        processCommandLineArgs(args, job);

        // for multiple runs remove be.iminds.cloudspeller.output folder of previous run
        FileSystem fs = FileSystem.get(job.getConfiguration());
        fs.delete(new Path(output), true); // true stands for recursively deleting the folder you gave

        // setting IO

        try {
            FileInputFormat.setInputPaths(job, new Path(input));
            FileOutputFormat.setOutputPath(job, new Path(output));
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

        return job.waitForCompletion(true) ? 0 : 1;

    }

    private static void processCommandLineArgs(String[] args, Job job) {


       List<String> other_args = new ArrayList<String>();

        for (int i = 0; i < args.length; ++i) {
            try {

                    other_args.add(args[i]);

            } catch (NumberFormatException except) {
                System.out.println("ERROR: Integer expected instead of " + args[i]);
                printUsage();
            }
        }

        if (other_args.size()!=2){
            System.out.println("ERROR: Wrong number of parameters: " + args.length + " instead of 2.");
            printUsage();
        }

        input = other_args.get(0);
        output = other_args.get(1);

    }

    private static int printUsage() {
        System.out.println("bin/hadoop jar input/ output/");
        return -1;

    }

    public static void main(String[] args) throws Exception {

        int res = 0;
        try {
            res = ToolRunner
                    .run(new Configuration(), new DBMotifCountsTableExactThresholds(), args);
        } catch (Exception ex) {
            System.err.println("exception: " + ex.getMessage());
        }
        System.exit(res);
    }
}
