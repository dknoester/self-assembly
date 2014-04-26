int[] INFO_AREA= {
  10, 430, 400, 400
};
PFont FONT;
int FONT_SIZE=16;

int cols;  
int rows;
int pages; 

class Experiment {
  int currentUpdate;
  boolean paused;

  Datafile df;

  Experiment() {    
    FONT = createFont("SansSerif", FONT_SIZE);
    currentUpdate = 0;
    background(255);
  }
  
  void addDatafile(Datafile d) {
    df = d;
    cols = df.xdim;
    rows = df.ydim;
    pages = df.zdim;
  }
  
  int getCurrentUpdate () {
    return currentUpdate; 
  }
  
  void reset() {
    currentUpdate=0;
    paused = true;
    redraw();
  }

  void step(int direction) {
    // are we doing anything?
    if((direction > 0) && atEnd()) {
      pause();
      return;
    }
    if((direction < 0) && atBegin()) {
      pause();
      return;
    }
    currentUpdate += direction;
    redraw();
  }
  
  void togglePause() {
    if(paused) {
      unpause();
    } 
    else {
      pause();
    }
  }

  void pause() {
    paused = true;
    noLoop();
  }

  void unpause() {
    paused = false;
    loop();
  }
  
  boolean atBegin() {
    return (currentUpdate <= 0);
  }

  boolean atEnd() {    
    return (df.atEnd(currentUpdate));
  }
  
 void draw() {
    String[] cMap = df.getStep(currentUpdate);

    hint(ENABLE_DEPTH_TEST);

    println("Current update... " + currentUpdate);
    int cellCount = 0;
    int xoffset = 150; 
    int yoffset = height/2;
    int zoffset = 0;
    int edge = 50;
    
    camera(mouseX, mouseY, (height/2) / tan(PI/6), width/2, height/2, 0, 0, 1, 0);

    background(255);

  
      for (int k = 0; k < pages; k++) { 
        for (int i = 0; i < cols; i++) {
          for (int j = 0; j < rows; j++) {
            int tempC = int(cMap[cellCount]);
            if (tempC == 1) {
              stroke(255); 
              fill(BLACK[0], BLACK[1], BLACK[2], 100);
              
            } else {
              //fill(WHITE[0], WHITE[1], WHITE[2], 50);
              stroke(180);
              noFill();
            }
            pushMatrix(); 
            int x = xoffset + i * edge; 
            int y = yoffset + j * edge;
            int z = zoffset + k * -edge;
            translate(y, x, z); 
            box(edge);
            popMatrix();
    
            ++cellCount;
          }
        }
      }
  }
  
}

