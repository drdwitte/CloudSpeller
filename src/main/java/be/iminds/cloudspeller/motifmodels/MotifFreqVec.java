package be.iminds.cloudspeller.motifmodels;

import org.apache.hadoop.io.BytesWritable;



public class MotifFreqVec {
	
	private static final String delimiter="-";
	private static MotifFactory factory;
	private Motif motif;
	private FreqVec vec;
	
		
	
	
	
	//CONSTRUCTORS
	
	public MotifFreqVec(){
		vec=new FreqVec();
	}
	
	public MotifFreqVec(Motif motif, FreqVec vec){
		setMotif(motif);
		setVec(vec);
		
	}
	
	//SETTERS
	
	public static void setMotifFactory(MotifFactory mFac){
		factory=mFac;
	}
	
	public void setMotif(Motif m){
		this.motif=m;
	}
	
	public void setVec(FreqVec v){
		this.vec=v;
	}
	
	
	//GETTERS
	public Motif getMotif(){
		return motif;
	}
	
	public FreqVec getVec(){
		return vec;
	} 
	
	
	
	
	
	 //METHODS
	
	public boolean isInitialized(){
		return motif!=null;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof MotifFreqVec){
			MotifFreqVec mFreq = (MotifFreqVec) obj;
			return mFreq.motif.equals(this.motif) && mFreq.vec.equals(this.vec);
		}
		return false;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(motif.toString());
		sb.append("\t");
		sb.append(vec.toString());
		return sb.toString();
	}
	
	
	public static MotifFreqVec createMotifFreqVecFromString(String s){
		String [] splits=s.split(delimiter);
		Motif motif=factory.createMotifFromString(splits[0]);
		FreqVec vec = FreqVec.createFreqVecFromString(splits[1]);
		return new MotifFreqVec(motif, vec);
	}
	
	public String createStringRepresentation(){
		return factory.createStringRepresentation(motif)+delimiter
						+vec.createStringRepresentation();
	}
	
	public BytesWritable createBytesRepresentation(){
		
		byte [] motifBytes = factory.createBytesRepresentation(motif);
		byte [] vecBytes   = vec.createBytesRepresentation();
		
		//join to byte arrays
		byte [] total = new byte[motifBytes.length+vecBytes.length];
		for (int i=0; i<motifBytes.length; i++){
			total[i]=motifBytes[i];
		}
		
		int c = motifBytes.length;
		for (int i=0; i<vecBytes.length; i++){
			total[c++]=vecBytes[i];
		}
		return new BytesWritable(total);
	}
	
	public static MotifFreqVec createMotifFreqVecFromBytes(BytesWritable b){
		
		byte [] total = b.getBytes();
		int offset=0;
		Motif motif = factory.createMotifFromBytes(offset,total);
		offset = factory.getNumberOfBytesForMotif();
		FreqVec vec = FreqVec.createFreqVecFromBytes(offset,total);
		return new MotifFreqVec(motif,vec);
	}
	

	


	

}
