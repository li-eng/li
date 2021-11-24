import RPi.GPIO as GPIO
import time

# w1,w2,w3,w4,w5,w6 = 0,1,0,0,0,1,0   细分1600，电流3A，电压12V


dir = 21#dianji
pul =20


#Gnd=19;
# 初始化树莓派引脚，设置树莓派的引脚为输出状态
def setup():
    
    GPIO.setwarnings(False)
    GPIO.setmode(GPIO.BCM)
    
    GPIO.setup(dir, GPIO.OUT)
    GPIO.setup(pul, GPIO.OUT)
    #GPIO.setup(Gnd, GPIO.OUT)
    #GPIO.output(Gnd, GPIO.LOW)

def setStep(w1, w2):
    GPIO.output(dir, w1)
    GPIO.output(pul, w2)
    
def forward():
    #setup()
    #while 1:
    t = 0.001
    for i in range(1, 800):
        setStep(1, 1)
        time.sleep(t)
        setStep(1, 0)
        t = t- 0.000000625
        #time.sleep(0.0005)
    #t = 0.005-0.000001*800
    for i in range(1, 800):
        setStep(1, 1)
        time.sleep(t)
        setStep(1, 0)
    for i in range(1, 800):
        setStep(1, 1)
        time.sleep(t)
        setStep(1, 0)
        t = t + 0.000000625
def backward():
    #setup()
    #while 1:
    t = 0.001
    for i in range(1, 800):
        setStep(0, 0)
        time.sleep(t)
        setStep(0, 1)
        t = t- 0.000000625
    #t = 0.002- 0.000001*800
    for i in range(1, 800):
        setStep(0, 0)
        time.sleep(t)
        setStep(0, 1)
    for i in range(1, 800):
        setStep(0, 0)
        time.sleep(t)
        setStep(0, 1)
        t = t + 0.000000625
#         time.sleep(t)
i=1    
while i==1:
    setup()
    forward()
    time.sleep(1)
    backward()
    i = i -1
