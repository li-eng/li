import RPi.GPIO as GPIO
import time

der1=21#4988-1
step1=20

key1=24


# 初始化树莓派引脚，设置树莓派的引脚为输出状态
def setup():
    GPIO.setwarnings(False)
    GPIO.setmode(GPIO.BCM)
    GPIO.setup(der1, GPIO.OUT)#setup_motor1
    GPIO.setup(step1, GPIO.OUT)

    GPIO.setup(key1, GPIO.IN, pull_up_down=GPIO.PUD_UP)

#gan_down
def down_loop():
    setup()
    GPIO.output(der1, GPIO.LOW)#motor1_down
    for i in range(0,800):
        GPIO.output(step1, GPIO.HIGH)
        time.sleep(0.001)
        GPIO.output(step1, GPIO.LOW)
        time.sleep(0.001)
    

#gan_up
def up_loop():
    setup()


    GPIO.output(der1, GPIO.HIGH) #motor2_up  
    for i in range(0,800):
        GPIO.output(step1, GPIO.HIGH)
        time.sleep(0.001)
        GPIO.output(step1, GPIO.LOW)
        time.sleep(0.001)
       
def stop():
    setup()
    GPIO.output(der1, GPIO.HIGH)
    GPIO.output(step1, GPIO.HIGH)
    
def loop():
    while True:
        down_loop()
        if GPIO.input(key1) == 0:
            up_loop()
            break
        

up_loop()

#down_loop()
