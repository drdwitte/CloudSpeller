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
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import java.io.IOException;
import java.util.*;

/**
 * Created by ddewitte on 08.12.14.
 */
public class DBMotifDeepFilterCFTKS extends Configured implements Tool {

    private static String input;
    private static String output;

    public static class DeepFilterMapperCFTKS extends Mapper<LongWritable,Text,NullWritable,Text> {

        private static final int [] t = {15,50,60,70,90,95};


        private static int motifLength = 12;
        private static int degeneracy = 64;

        private static int [] F;
        private static double [] p;


        private static IntWritable mapOutputValue = new IntWritable(1);

        private static ConfidenceGraphRestrictions restrictions;

        @Override
        protected void setup(Context context) throws IOException,
                InterruptedException {
            BLS.initializeBLSConstants(t);
            FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());

            F = new int[FreqVec.getNumberOfIntervals()];
            p = new double[FreqVec.getNumberOfIntervals()];

            int cThr = context.getConfiguration().getInt("C",90);
            int blsThr = context.getConfiguration().getInt("BLS",15);
            int FThr = context.getConfiguration().getInt("F",10);
            motifLength = context.getConfiguration().getInt("k",12);
            degeneracy = context.getConfiguration().getInt("s",64);

            restrictions = new SimultaneousOccurrenceConfidenceAndBLSFiltering
                    (cThr, FThr, blsThr );
        }

        @Override
        protected void map(LongWritable key, Text value, Context context)
                throws IOException, InterruptedException {


            String strValue = value.toString();


            if (strValue == null || strValue.length()==0)
                return;


            int k = getMotifLength(strValue);


            if (k!=motifLength)
                return;

            int s = getMotifDegeneracy(strValue);
            if (s!=degeneracy)
                return;



            GeneralToolbox.parseConfidenceGraphValues(strValue,F,p);



            if (restrictions.checkRestrictions(F,p)){
                    context.write(NullWritable.get(),value);
            }

        }



        private static void pr(String s){
            System.err.println(s);
        }

        private int getMotifLength(String confGraphRecord){

            return confGraphRecord.indexOf('\t');

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






    @Override
    public int run(String[] args) throws Exception {


        //NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
        Job job = Job.getInstance(getConf(), "DBMotifDeepFilter");

        job.setJarByClass(DBMotifDeepFilterCFTKS.class);

        //SETTING MAPREDUCE CLASSES

        // IO
        job.setInputFormatClass(TextInputFormat.class);
        job.setOutputFormatClass(TextOutputFormat.class);

        // (K3,V3)
        job.setMapOutputKeyClass(NullWritable.class);
        job.setMapOutputValueClass(Text.class);

        // set hadoop methods (MAP ONLY)
        job.setMapperClass(DeepFilterMapperCFTKS.class);
        job.setNumReduceTasks(0);

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


        int confidenceCutoff=90;
        int familyCutoff=10;
        int minBLS=15;
        int k = 12;
        int s = 64;


        List<String> other_args = new ArrayList<String>();

        for (int i = 0; i < args.length; ++i) {
            try {
                if ("-C".equals(args[i])) {
                    confidenceCutoff = Integer.parseInt(args[++i]);
                } else if ("-F".equals(args[i])){
                    familyCutoff = Integer.parseInt(args[++i]);
                } else if ("-BLS".equals(args[i])){
                    minBLS = Integer.parseInt(args[++i]);
                } else if ("-k".equals(args[i])){
                    k = Integer.parseInt(args[++i]);
                } else if ("-s".equals(args[i])){
                    s = Integer.parseInt(args[++i]);

                } else {
                    other_args.add(args[i]);
                }
            } catch (NumberFormatException except) {
                System.out.println("ERROR: Integer expected instead of " + args[i]);
                printUsage();
            }
        }



        job.getConfiguration().setInt("C",confidenceCutoff);
        job.getConfiguration().setInt("F",familyCutoff);
        job.getConfiguration().setInt("BLS",minBLS);
        job.getConfiguration().setInt("k",k);
        job.getConfiguration().setInt("s",s);



        if (other_args.size()!=2){
            System.out.println("ERROR: Wrong number of parameters: " + args.length);
            System.out.println("Current parameter values: C="+confidenceCutoff+" F="+ familyCutoff+" BLS="+minBLS+" k="+k+" s="+s);
            printUsage();
        }

        input = other_args.get(0);
        output = other_args.get(1);
    }


    private static int printUsage() {
        System.out.println("bin/hadoop jar DeepFilter.jar -C 90 -F 10 -BLS 15 -k 12 -s 64 input/ output/");
        return -1;

    }

    public static void main(String[] args) throws Exception {

        int res = 0;
        try {
            res = ToolRunner
                    .run(new Configuration(), new DBMotifDeepFilterCFTKS(), args);
        } catch (Exception ex) {
            System.err.println("exception: " + ex.getMessage());
        }
        System.exit(res);
    }
}
