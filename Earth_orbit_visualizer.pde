// Gilles Regamey 2024

int scrsize = 1000;
float camdist = 0.5;
int mwcnt = 0;
float cam_long = 0;
float cam_lat = 0;
float xoff = 0;
float yoff = 0;
int earthrad = 6371;
int sunrad = 696340;
int au = 149600000;
float mu = 398600.4418;
int loadedsat = 0;
float spdfactor = 2048; 
int currentsat = 0;
int maxsat = 1000;

PImage earth,suntext,cosmicbg;
PShape globe,sun,bg;
Satelite[] sats;

class Satelite{
  float a,e,i,raan,argp,t0,epoch,n,meananom;
  String name,catalog_nbr,designation,info;
  char classification;
  int epochyr;
  Satelite(){
    this.name = "Uninitialized-------------------";
    this.a = 0;
    this.e = 0;
    this.i = 0;
    this.raan = 0;
    this.argp = 0;
    this.t0 = 0;
    this.catalog_nbr = "00000000";
    this.designation = "00000000";
    this.classification = 'U';
    this.epoch = 0;
    this.epochyr = 0;
    this.n = 0;
  }
  Satelite(Satelite s){
    this.name = s.name;
    this.a = s.a;
    this.e = s.e;
    this.i = s.i;
    this.raan = s.raan;
    this.argp = s.argp;
    this.t0 = s.t0;
    this.catalog_nbr = s.catalog_nbr;
    this.classification = s.classification;
    this.epoch = s.epoch;
    this.epochyr = s.epochyr;
    this.n = s.n;
  }
  void printinfo(){
    print("______________________________\n");
    print("Name: ",this.name,"\n");
    print("Catalog Number: ",this.catalog_nbr,"\n");
    print("Classification: ",this.classification,"\n");
    print("Designation: ",this.designation,"\n");
    print("Epoch: ",this.epochyr,"-",this.epoch,"\n");
    print("a: ",this.a);
    print(" e: ",this.e);
    print(" i: ",this.i,"\n");
    print("raan: ",this.raan);
    print(" argp: ",this.argp);
    print(" mAn ",this.meananom,"\n");
    print("revs/day: ",this.n,"\n");
  }
  void draw(){
    drawOrbit(this.a/1000,this.e,this.i,this.raan,this.argp,0,100,int(this.catalog_nbr));
  }

  void drawsat(float secs){
    float x,y,z;
    strokeWeight(5/camdist);
    int i = int(this.catalog_nbr);
    stroke((i*1135) % 255, (i*2742) % 255, (i*3171) % 255,255);
    
    float meananom = this.meananom + this.n*(secs - this.t0);
    meananom = meananom % (2*PI);
    float trueanom = trueAnomaly(meananom,this.e,20,10);  // not optimal as it is valid only for low eccentricity orbits
    float doff = this.argp;
    float ang = trueanom + this.argp;
    float r = (this.a/1000)*(1-this.e*this.e)/(1+this.e*cos(ang-doff));
    float diffang = ang ;
    x = r*(cos(this.raan)*cos(diffang) - sin(this.raan)*sin(diffang)*cos(this.i));
    y = r*(sin(this.raan)*cos(diffang) + cos(this.raan)*sin(diffang)*cos(this.i));
    z = r*sin(diffang)*sin(this.i);
    line(0,0,0,x,y,z);
    }
}

int getTLEs(String filename){
  String[] lines = loadStrings(filename);
  int satcnt = 0;
  Satelite tempsat = new Satelite();
  print("Loading TLEs\n");
  for(int i = 0; i < lines.length; i++){
    String line = lines[i];
    if(line.charAt(0) == '1' && line.charAt(1) == ' '){
      tempsat.catalog_nbr = line.substring(2,7);
      tempsat.classification = line.charAt(7);
      tempsat.designation = line.substring(9,17);
      tempsat.epochyr = int(line.substring(18,20));
      tempsat.epoch = float(line.substring(20,32));
      tempsat.t0 = tempsat.epoch*24*60*60;
    }else if(line.charAt(0) == '2' && line.charAt(1) == ' '){
      tempsat.i = radians(float(line.substring(8,16)));
      tempsat.raan = radians(float(line.substring(17,25)));
      tempsat.e = float("0."+line.substring(26,33));
      tempsat.argp = radians(float(line.substring(34,42)));
      tempsat.meananom = radians(float(line.substring(43,51)));
      tempsat.n = float(line.substring(52,63))*2*PI;  // rads per day
      tempsat.n = tempsat.n/(24*60*60);                   // rads per sec
      tempsat.t0 = tempsat.t0 % (2*PI/tempsat.n);
      tempsat.a = pow(mu*(1/(tempsat.n*tempsat.n)),1.0/3.0);
      if(satcnt >= sats.length){
        i = lines.length;
        print("\nToo many satellites\n");
      }else{
        sats[satcnt] = new Satelite(tempsat);
        print("\r",satcnt,"/",lines.length," - ",tempsat.name);
        float degprad = 180/PI;
        sats[satcnt].info = (tempsat.name + "\n" + tempsat.catalog_nbr + " " + tempsat.classification + " " + tempsat.designation + "\nepoch:" + tempsat.epochyr + "Y " + tempsat.epoch + "days\na:" + tempsat.a + " km, e:" + tempsat.e + " , i:" + degprad*tempsat.i + "deg\nRAAN:" + degprad*tempsat.raan + " deg, w:" + degprad*tempsat.argp + " deg, M0:" + degprad*tempsat.meananom + "deg\nn:" + (tempsat.n*12*60*60)/PI + " revs/day");
        satcnt++;
      }
    }else{
      tempsat.name = line;
    }
  }
  print("Loaded ",satcnt," satellites\n");
  return satcnt;
}

float besselJ(int n, float x) {
  int numSteps = 1000; // Number of steps in numerical integration
  float step = PI / numSteps; // Integration step size
  float sum = 0.0;

  // Numerical integration
  for (int i = 0; i < numSteps; i++) {
    float tau = i * step;
    float cosTerm = cos(n * tau - x * sin(tau));
    sum += cosTerm * step;
  }

  return sum / PI;
}

float trueAnomaly(float M, float e, int maxK, int maxN) {
  if (e <= 0 || e >= 1) {
    println("Eccentricity e must be between 0 and 1.");
    return Float.NaN;
  }

  float beta = (1 - sqrt(1 - e * e)) / e;
  float nu = M;

  // Sum over k
  for (int k = 1; k <= maxK; k++) {
    float sumN = 0;
    
    // Sum over n
    for (int n = -maxN; n <= maxN; n++) {
      float besselTerm = besselJ(n, -k * e);
      float betaTerm = pow(beta, abs(k + n));
      sumN += besselTerm * betaTerm;
    }
    
    // Add to the true anomaly
    nu += (2.0 / k) * sumN * sin(k * M);
  }
  
  return nu;
}

void drawOrbit(float a, float e, float i, float raan, float argp, float man,int opacity,int id){
  float x,y,z;
  float px = 0,py = 0,pz = 0;
  strokeWeight(1/camdist);
  if(id == 0){
    stroke(255,255,255,opacity);
  }else{
    stroke((id*1135) % 255, (id*2742) % 255, (id*3171) % 255,opacity);
  }
  noFill();
  beginShape();
  vertex(0,0,0);
  float doff = raan+argp;
  for(float ang = doff; ang <= doff+TWO_PI; ang += 0.05){
    float r = a*(1-e*e)/(1+e*cos(ang-doff));
    float diffang = ang - raan;
    x = (r )*(cos(raan)*cos(diffang) - sin(raan)*sin(diffang)*cos(i));
    y = r*(sin(raan)*cos(diffang) + cos(raan)*sin(diffang)*cos(i));
    z = r*sin(diffang)*sin(i);
    if(ang == doff){
      px = x;
      py = y;
      pz = z;
    }
    vertex(x,y,z);
  }
  vertex(px,py,pz);
  endShape();
} 

void drawAxis(){
  strokeWeight(10/camdist);
  stroke(255,0,0);
  noFill();
  line(0,0,0,42.164,0,0);
  stroke(0,255,0);
  line(0,0,0,0,42.164,0);
  stroke(0,100,255);
  line(0,0,0,0,0,42.164);
}

void mouseWheel(MouseEvent event){
  if(keyPressed && key == 's'){
    currentsat += event.getCount();
    if(currentsat < 0){
      currentsat = 0;
    }else if(currentsat >= loadedsat){
      currentsat = loadedsat-1;
    }
  }else if(keyPressed && key == 't'){
    if(event.getCount() > 0){
      spdfactor *= 2;
    }else{  
      spdfactor /= 2;
    }
    if(spdfactor < 1){
      spdfactor = 1;
    }
  }else{
    mwcnt -= event.getCount();
    camdist = exp(mwcnt);
    if(camdist < 0.1){
      camdist = 0.1;
      mwcnt = int(log(camdist));
    }else if(camdist > 10*au){
      camdist = 10*au;
      mwcnt = int(log(camdist));
    }
  }
}

void mouseDragged(){
  if(mouseButton == LEFT){
    cam_long += radians(mouseX - pmouseX);
    cam_lat += radians(mouseY - pmouseY);
  }else{
    xoff += (mouseX - pmouseX)*(1/camdist);
    yoff -= (mouseY - pmouseY)*(1/camdist);
  }
}

void displayinfo(){
  textAlign(LEFT);
  textSize(32);
  text("entry n°" + currentsat + " - " + sats[currentsat].info, 64, 64);
  text("scroll: zoom, +s: change satellite, +t: change timespeed", 64, height-64);
  text("t * " + spdfactor + "days: " + times/(24*60*60), width-32*10, height-64);
}

void setup(){
  size(1940,1080,P3D);
  noStroke();
  noFill();
  earth = loadImage("earth.jpg");
  cosmicbg = loadImage("stars2.jpg");
  suntext = loadImage("sun.jpg");
  globe = createShape(SPHERE, earthrad/1000);
  sun = createShape(SPHERE, sunrad/1000);
  bg = createShape(SPHERE, 3*au);
  textureMode(NORMAL);
  globe.setTexture(earth);
  sun.setTexture(suntext);
  bg.setTexture(cosmicbg);
  sats = new Satelite[maxsat];
  //loadedsat = getTLEs("tles.txt");
  loadedsat = getTLEs("https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle");
}

float times = 0;
int lastmillis = millis();
void draw(){
  background(0);
  displayinfo();
  float dtsec = (millis()-lastmillis)/1000.0;
  lastmillis = millis();
  times += dtsec*spdfactor;
  perspective(PI/3,width*1.0/height,0.1,au);
  scale(1,-1,1);
  translate(width/2,-height/2,0);
  scale(camdist);

  rotateZ(cam_long);
  //translate(xoff,yoff,0);
  translate(cos(cam_long)*xoff + sin(cam_long)*yoff,-sin(cam_long)*xoff + cos(cam_long)*yoff,0);
  rotate(cam_lat,cos(cam_long),-sin(cam_long),0);

  // in solar ref frame
  rotateX(-PI/2);
  stroke(255,100,0,50);
  shape(sun,au/1000,0);
  shape(bg);
  rotateX(PI/2);
  pointLight(255,255,255, -au/1000, 0, 0);

  strokeWeight(1/camdist);
  stroke(255,255,255,100);
  ellipse(au/1000, 0, au/500, au/500);

  rotateX(radians(23.5));
  // in earth ref frame
  drawAxis();

  // GEO for reference
  drawOrbit(42.164,0,radians(0),radians(0),radians(0),0,50,0);
  //draw the satellites
  sats[currentsat].draw();
  sats[currentsat].drawsat(times);

  // draw the earth
  float radpsec = 2*PI/(24*60*60);
  rotateZ(times*radpsec);
  rotateX(-PI/2);
  shape(globe);
}
