package be.iminds.cloudspeller.postprocessing_Single;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.output.SimultaneousOccurrenceConfidenceAndBLSFiltering;
import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
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
import java.util.List;

/**
 * Created by ddewitte on 08.12.14.
 */
public class DBMotifTranslator extends Configured implements Tool {

    private static String input;
    private static String output;

    public static class TranslatorMapper extends Mapper<LongWritable,Text,Text,Text> {

        private static final int [] t = {15,50,60,70,90,95};

        private static Text outputKey = new Text("");

        @Override
        protected void setup(Context context) throws IOException,
                InterruptedException {
            BLS.initializeBLSConstants(t);
            FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());

        }

        @Override
        protected void map(LongWritable key, Text value, Context context)
                throws IOException, InterruptedException {


            String strValue = value.toString();


            if (strValue == null || strValue.length()==0)
                return;


            String strKey = calculateKey(strValue.split("\t")[0]);


            context.write(outputKey,value);


        }

        private static String calculateKey(String permGroupKey){

            String translate = mapToDC(permGroupKey);

            String comp_trans = mapToComplement(translate);

            return comp_trans;

        }

        private static String mapToDC(String s){
            StringBuilder mappedStr = new StringBuilder("");
            for (int i=0; i<s.length(); i++){
                if (baseChars.indexOf(s.charAt(i))<0){
                    mappedStr.append(dc);
                } else {
                    mappedStr.append(s.charAt(i));
                }
            }
            return mappedStr.toString();
        }

        private static String mapToComplement(String s){
            String mappedStr = "";
        return mappedStr;
    }

        private static void pr(String s){
            System.err.println(s);
        }

        private static final String baseChars = "ACGT";

        private static final char dc = 'N';
    }


    static class TranslatorGrouper extends Reducer<Text,Text,Text,Text> {


        private static Text outputValue = new Text("");

        @Override
        protected void reduce(Text key, Iterable<Text> values,
                              Context context)
                throws IOException, InterruptedException {


            context.write(key,outputValue);

        }

        private static void pr(String s){
            System.err.println(s);
        }

    }




    @Override
    public int run(String[] args) throws Exception {


        //NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
        Job job = Job.getInstance(getConf(), "DBMotifTranslator");

        job.setJarByClass(DBMotifTranslator.class);

        //SETTING MAPREDUCE CLASSES

        // IO
        job.setInputFormatClass(TextInputFormat.class);
        job.setOutputFormatClass(TextOutputFormat.class);

        // (K3,V3)
        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(Text.class);

        // set hadoop methods (MAP ONLY)
        job.setMapperClass(TranslatorMapper.class);
        job.setReducerClass(TranslatorGrouper.class);


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
            System.out.println("ERROR: Wrong number of parameters: " + args.length);
            printUsage();
        }

        input = other_args.get(0);
        output = other_args.get(1);
    }


    private static int printUsage() {
        System.out.println("bin/hadoop jar MotifTranslator.jar input/ output/");
        return -1;

    }

    public static void main(String[] args) throws Exception {

        int res = 0;
        try {
            res = ToolRunner
                    .run(new Configuration(), new DBMotifTranslator(), args);
        } catch (Exception ex) {
            System.err.println("exception: " + ex.getMessage());
        }
        System.exit(res);
    }
}
