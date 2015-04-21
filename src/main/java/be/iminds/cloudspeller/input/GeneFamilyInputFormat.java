package be.iminds.cloudspeller.input;

import java.io.IOException;

import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;


public class GeneFamilyInputFormat extends FileInputFormat<LongWritable,GeneFamily> {

	//FIXME this implies that every file will be treated by a different mapper so no optimization
	//is possible!! The hadoop reference refers to CombineFileFormat as a way out of this.
	//the reason for using one file per family is that the splits need to correspond the 
	//an integer number of records which is not predictable since every gene family can have
	//a different number of genes!

	@Override
	public RecordReader<LongWritable, GeneFamily> createRecordReader(
			InputSplit arg0, TaskAttemptContext arg1) throws IOException,
			InterruptedException {
		return new GeneFamilyRecordReader();
	}

}
