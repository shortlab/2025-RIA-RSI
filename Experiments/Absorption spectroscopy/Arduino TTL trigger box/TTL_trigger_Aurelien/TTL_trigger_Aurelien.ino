const int Button_pin = 3;
const int TTL_in = 8;
const int TTL_out1 = 5; //yellow
const int TTL_out2 = 6; //orange
const int LED_pin = 13;

int buttonState =  0;
int sweepState =  0;
int pre_sweepState = 0;
int pre_buttonState =  0;
bool sweepStarts = false;
bool toggleButton = false;
bool changerequest = false;
bool isON = false;

unsigned long lastToggle = 0;
const unsigned long timeout = 2000; // 2s in ms

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(Button_pin,INPUT);
  pinMode(TTL_in,INPUT);
  pinMode(TTL_out1,OUTPUT);
  pinMode(TTL_out2,OUTPUT);
  pinMode(LED_pin,OUTPUT);

  digitalWrite(LED_pin, LOW);
  digitalWrite(TTL_out1, LOW);
  digitalWrite(TTL_out2, LOW);
}

void loop() {
  // put your main code here, to run repeatedly:
  buttonState = digitalRead(Button_pin);
  sweepState = digitalRead(TTL_in);

  sweepStarts = (sweepState - pre_sweepState)==1;

  toggleButton = (buttonState - pre_buttonState)==1;
  if (toggleButton and (millis() - lastToggle >= timeout) ) {
    changerequest = true;
    lastToggle = millis();
  }

  if (sweepStarts && changerequest){
    if (!isON){
      digitalWrite(LED_pin, HIGH);
      digitalWrite(TTL_out1, HIGH);
      digitalWrite(TTL_out2, HIGH);
    }
    if (isON){
      digitalWrite(LED_pin, LOW);
      digitalWrite(TTL_out1, LOW);
      digitalWrite(TTL_out2, LOW);
    }
    isON = !isON;
    Serial.println(isON);
    changerequest = false;
  }

  pre_sweepState = sweepState;
  pre_buttonState = buttonState;
}
