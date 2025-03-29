const int Button_pin = 3;
const int TTL_in = 8;
const int TTL_out1 = 6; //orange
const int TTL_out2 = 5; //yellow
const int LED_pin = 13;

int buttonState =  0;
int sweepState =  0;
int pre_sweepState = 0;
int pre_buttonState =  0;
bool sweepStarts = false;
bool sweepEnds = false;
int countStart = 0;
unsigned long avg_sweep_time = 0;
unsigned long avg_sweep_period = 0;
unsigned long LastStartTime = 0;
unsigned long LastEndTime = 0;
unsigned long curtime = 0;

bool toggleButton = false;
bool changerequest = false;
bool isON = false;

unsigned long lastToggle = 0;
const unsigned long timeout = 2000; // 2s in ms
const unsigned long timeavg = 30000; // 2s in ms

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(Button_pin,INPUT);
  pinMode(TTL_in,INPUT);
  pinMode(TTL_out1,OUTPUT);
  pinMode(TTL_out2,OUTPUT);
  pinMode(LED_pin,OUTPUT);

  digitalWrite(LED_pin, HIGH);
  digitalWrite(TTL_out1, LOW);
  digitalWrite(TTL_out2, LOW);
}

void loop() {
  // put your main code here, to run repeatedly:
  
  buttonState = digitalRead(Button_pin);
  sweepState = digitalRead(TTL_in);
  curtime = millis();

  sweepStarts = (sweepState - pre_sweepState)==1;
  sweepEnds = (sweepState - pre_sweepState)==-1;

  toggleButton = (buttonState - pre_buttonState)==1;
  if (toggleButton and (millis() - lastToggle >= timeout) and !(isON) ) {
    changerequest = true;
    lastToggle = millis();
  }

  if ((changerequest) and (sweepStarts)) {
    isON = true;
    countStart = 0;
    avg_sweep_time = 0;
    avg_sweep_period = 0;
    changerequest = false;
  }

  if ((millis() - lastToggle >= timeavg) and isON) {
    isON = false;
    Serial.println("END OF AVERAGE");
  }

  if (isON) {
    if (sweepStarts) {
      countStart = countStart+1;
      if (countStart==1) {
        LastStartTime = curtime;
      } else {
        avg_sweep_period = ( avg_sweep_period*(countStart-2) + (curtime-LastStartTime) ) / (countStart-1);
        LastStartTime = curtime;
        Serial.print("Average sweep period (ms) : ");
        Serial.println(avg_sweep_period);
      }
    }
    if (sweepEnds) {
      avg_sweep_time = ( avg_sweep_time*(countStart-1) + (curtime-LastStartTime) ) / (countStart);
      Serial.print("Average sweep time (ms) : ");
      Serial.println(avg_sweep_time);
    }
  }

  

  pre_sweepState = sweepState;
  pre_buttonState = buttonState;
}
