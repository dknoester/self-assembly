// 2D Array of objects
Experiment expr;


String filename = "/Users/heathergoldsby/Desktop/processing-projects/ca_1d.dat";


void setup() {
  size(640,360,P3D);
  frameRate(5);
  //lights();
  background(0);
  smooth();
  expr = new Experiment();
  expr.addDatafile(new Datafile(filename, expr));
  
}

void draw() {

  if (!expr.paused) {
    expr.step(1);
  }
  expr.draw();
} 

void keyPressed() {
  switch(key) {
 
  case 'P':
  case 'p': 
    {
      expr.togglePause();
      break;
    }
      case 'Q':
  case 'q': 
    {
      exit();
    }
  case 'r': 
    {
      expr.reset();
      break;
    }
  case 'S': 
    {
      expr.step(-1);
      break;
    }
  case 's': 
    {
      expr.step(1);
      break;
    }
  }
}
  

