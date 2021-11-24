import RPi.GPIO as GPIO
import time
import threading

dir1 = 21#1_dianji
step1 =20

dir2 = 16#2_dianji
step2 =12

dir3 = 26#3_dianji
step3 =19

dir4 = 6#4_dianji
step4 =5

dir5 = 23#5_dianji
step5 =24

dir6 = 7
step6 = 8#6_dianji

motor1_delay = 0.003
motor1_steps = int((395-385)*200)
motor2_steps = int((395-385)*200)
motor3_steps = int((395-390)*200)
motor4_steps = int((413-395)*200)
motor5_steps = int((413-395)*200)
motor6_steps = int((395-390)*200)

class Motor:
    def __init__(self, d1, s1):
        self.d1= d1
        self.s1 = s1
        
        GPIO.setwarnings(False)
        GPIO.setmode(GPIO.BCM)
        GPIO.setup(self.d1, GPIO.OUT)
        GPIO.setup(self.s1, GPIO.OUT)
     
    def loop(self, delay_loop, steps_loop):
        print("forward...")
        Motor.forward(self, delay_loop, steps_loop)
                
        print("stop...")
        Motor.stop(self)
        time.sleep(2)
        
        print("backward...")
        Motor.backward(self, delay_loop, steps_loop)
        
        print("stop...")
        Motor.stop(self)
        time.sleep(2)
        
    def setStep(self, w1, w2):
        GPIO.output(self.d1, w1)
        GPIO.output(self.s1, w2)
        
    def forward(self, delay, steps):
        for i in range(0, steps):
            Motor.setStep(self, 1, 1)
            time.sleep(delay)
            Motor.setStep(self, 1, 0)
            '''delay -= 0.00000125
        for i in range(0, steps):
            Motor.setStep(self, 1, 1)
            time.sleep(delay)
            Motor.setStep(self, 1, 0)
        for i in range(0, steps):
            Motor.setStep(self, 1, 1)
            time.sleep(delay)
            Motor.setStep(self, 1, 0)
            delay += 0.00000125 '''   
    def backward(self, delay, steps):
        for i in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)
            '''delay -= 0.00000125
        for i in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)
        for i in range(0, steps):
             Motor.setStep(self, 0, 1)
             time.sleep(delay)
             Motor.setStep(self, 0, 0)
             time.sleep(delay)
             delay. += 0.00000125 '''    
            
    def stop(self):
        Motor.setStep(self, 0, 0)
       
def destroy():
    GPIO.cleanup()
    
if __name__ == '__main__':
    motot1 = Motor(dir1, step1)
    motot2 = Motor(dir2, step2)
    motot3 = Motor(dir3, step3)
    motot4 = Motor(dir4, step4)
    motot5 = Motor(dir5, step5)
    motot6 = Motor(dir6, step6)
    try:
        while True:
            shuru = int(input('input:'))
            if shuru == 1:
                thread1 = threading.Thread(target=motot1.backward, args=(motor1_delay, motor1_steps))
                thread2 = threading.Thread(target=motot2.backward, args=(motor1_delay, motor2_steps))
                thread3 = threading.Thread(target=motot3.backward, args=(motor1_delay, motor3_steps))
                thread4 = threading.Thread(target=motot4.forward, args=(motor1_delay, motor4_steps))
                thread5 = threading.Thread(target=motot5.forward, args=(motor1_delay, motor5_steps))
                thread6 = threading.Thread(target=motot6.backward, args=(motor1_delay, motor6_steps))
                
                thread1.start()
                thread2.start()
                thread3.start()
                thread4.start()
                thread5.start()
                thread6.start()
    except KeyboardInterrupt:
        destroy()
    