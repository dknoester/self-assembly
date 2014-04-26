
class Datafile { 
  HashMap<Integer,String[]> data;
  int xdim;
  int ydim;
  int zdim; 
  Experiment expr;
  int end;


  Datafile(String filename, Experiment expr) {
    print("loading datafile: " + filename + "...");
    this.expr = expr;
    String[] file = loadStrings(filename);
    data = new HashMap<Integer,String[]>();
    int step = 0;
    Boolean init = false;


    for(int i=1; i<file.length; ++i) {
      String l = trim(file[i]);
      if((l.length()==0) || (l.charAt(0)=='#')) {
        continue;
      }
      String[] e = splitTokens(l);
      if (!init) {
         xdim = int(e[0]);
         ydim = int(e[1]);
         zdim = int(e[2]);
         init = true;
         println ("Initializing!");    
      } else {
          data.put(step, e);
          ++step;
      }
    }
    end = --step;
    println(" done.");
  }
  
  String[] getStep(int e) {
     return data.get(e); 
  }
  
  boolean atEnd(int update) {
    return update >= end;
  }
  
}


