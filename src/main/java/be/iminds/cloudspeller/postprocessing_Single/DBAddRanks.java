package be.iminds.cloudspeller.postprocessing_Single;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.output.ConfidenceGraphRestrictions;
import be.iminds.cloudspeller.output.MotifBLSRestrictionsWithCutoff;
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
public class DBAddRanks extends Configured implements Tool {

    private static String input;
    private static String output;

    public static class AddRanksMapper extends Mapper<LongWritable,Text,Text,Text> {

        private static final int [] t = {15,50,60,70,90,95};
        private static Text mapOutputKey = new Text("");

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

            mapOutputKey.set(strValue.split("\t")[0]);
            context.write(mapOutputKey,value);


            //output: motifPerm, (motif, F) => voor een gegeven BLS

       }



    }

    static class AddRanksReducer extends Reducer<Text,Text,NullWritable,Text> {

        private static final int [] t = {15,50,60,70,90,95};
        MotifBLSRestrictionsWithCutoff restr;
        private static Text outputValue = new Text("");
        private static final int maxPQSize = 1000;


        @Override
        protected void setup(Context context) throws IOException,
                InterruptedException {
            BLS.initializeBLSConstants(t);
            FreqVec.setNumberOfIntervals(BLS.getNumberOfIntervals());

            int[] F = new int[FreqVec.getNumberOfIntervals()];
            double[] p = new double[FreqVec.getNumberOfIntervals()];

            int FThr =   context.getConfiguration().getInt("F",10);
            int pThr =   context.getConfiguration().getInt("C",90);
            int blsThr = context.getConfiguration().getInt("BLS",15);


            restr = new MotifBLSRestrictionsWithCutoff(FThr, pThr, blsThr);

        }



        @Override
        protected void reduce(Text key, Iterable<Text> values,
                              Context context)
                throws IOException, InterruptedException {

            BestMotifContainer container = new BestMotifContainer(maxPQSize, restr);


            for (Text t : values){
                container.add(t.toString());
            }

            List<String> all = container.extractAll();
            int fams = 0;
            int currentRank = 1;
            int counter = 0;
            for (int i=all.size()-1; i>=0; i--){
                int currentFams = restr.getMaxFamilies(all.get(i));
                if (currentFams != fams){
                    currentRank=counter+1;
                    fams = currentFams;
                }

                outputValue.set(all.get(i)+"\t"+currentRank);
                context.write(NullWritable.get(), outputValue);
                counter++;

            }
        }


    }





    @Override
    public int run(String[] args) throws Exception {


        //NOTE: Job takes a deep copy of conf so to modify conf use job.getConfiguration().doSomething();
        Job job = Job.getInstance(getConf(), "DBAddRanks");

        job.setJarByClass(DBAddRanks.class);

        //SETTING MAPREDUCE CLASSES

        // IO
        job.setInputFormatClass(TextInputFormat.class);
        job.setOutputFormatClass(TextOutputFormat.class);

        // (K3,V3)
        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(Text.class);

        job.setOutputKeyClass(NullWritable.class);
        job.setOutputValueClass(Text.class);

        // set hadoop methods
        job.setMapperClass(AddRanksMapper.class);
        job.setReducerClass(AddRanksReducer.class);


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

        List<String> other_args = new ArrayList<String>();

        for (int i = 0; i < args.length; ++i) {
            try {
                if ("-C".equals(args[i])) {
                    confidenceCutoff = Integer.parseInt(args[++i]);
                } else if ("-F".equals(args[i])){
                    familyCutoff = Integer.parseInt(args[++i]);
                } else if ("-BLS".equals(args[i])){
                    minBLS = Integer.parseInt(args[++i]);
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




        if (other_args.size()!=2){
            System.out.println("ERROR: Wrong number of parameters: " + args.length);
            printUsage();
        }

        input = other_args.get(0);
        output = other_args.get(1);
    }


    private static int printUsage() {
        System.out.println("bin/hadoop jar AddRanks.jar -C 90 -F 10 -BLS 15 input/ output/");
        return -1;

    }

    public static void main(String[] args) throws Exception {

        int res = 0;
        try {
            res = ToolRunner
                    .run(new Configuration(), new DBAddRanks(), args);
        } catch (Exception ex) {
            System.err.println("exception: " + ex.getMessage());
        }
        System.exit(res);
    }
}
