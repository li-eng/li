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
step4 = 5

dir5 = 23
step5 = 24#5_dianji

dir6 = 7
step6 = 8#6_dianji


key1 = 27
key2 = 22
key3 = 25
key4 = 13
key5 = 4
key6 = 17

motorback_delay = 0.001
motorup_delay = 0.002
motor1_steps = 400*22

class Motor:
    def __init__(self, d1, s1):
        self.d1= d1
        self.s1 = s1
        
        #self.m1 = m1
        #self.m2 = m2
        #self.m3 = m3
        
        GPIO.setwarnings(False)
        GPIO.setmode(GPIO.BCM)
        GPIO.setup(self.d1, GPIO.OUT)
        GPIO.setup(self.s1, GPIO.OUT)
        #GPIO.setup(self.m1, GPIO.OUT)
        #GPIO.setup(self.m2, GPIO.OUT)
        #GPIO.setup(self.m3, GPIO.OUT)
        GPIO.setup(key1, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key2, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key3, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key4, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key5, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key6, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        
 #zihanshu       
    def setStep(self, w1, w2):
        GPIO.output(self.d1, w1)
        GPIO.output(self.s1, w2)
        
    def backward1(self, delay, steps):
        #delay1=delay
        while True:
            
        #for i in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
           
            if GPIO.input(key1) == 0:
                Motor.forward(self, 0.005, 800)
                break                            
    
    def backward2(self, delay, steps):
        #delay1=delay
        while True:
        #for i in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
       
            if GPIO.input(key2) == 0:
                Motor.forward(self, 0.005, 800)
                break
    def backward3(self, delay, steps):
        #delay1=delay
        while True:
        #for iorward( in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
        #step = steps
       
            if GPIO.input(key3) == 0:
                Motor.forward(self, 0.005, 800)
                break
    def backward4(self, delay, steps):
        #delay1=delay
        while True:
        #for i in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
        #step = steps
       
            if GPIO.input(key4) == 0:
                Motor.forward(self, 0.005, 800)
                break        
    def backward5(self, delay, steps):
        #delay1=delay
        while True:
        #for i in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
        #step = steps
       
            if GPIO.input(key5) == 0:
                Motor.forward(self, 0.005, 800)
                break
    def backward6(self, delay, steps):
        #delay1=delay
        while True:
        #for i in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
        #step = steps
       
            if GPIO.input(key6) == 0:
                Motor.forward(self, 0.005, 800)
                break      
    def forward(self, delay, steps):
        #while True:
        #step = 400#one circle
        for i in range(0, steps):
            Motor.setStep(self, 1, 1)
            time.sleep(delay)
            Motor.setStep(self, 1, 0)              
            '''delay -= 0.0000025
        for i in range(0, steps):
            Motor.setStep(self, 1, 1)
            time.sleep(delay)
            Motor.setStep(self, 1, 0)
        for i in range(0,2*step):
             Motor.setStep(self, 1, 1)
             time.sleep(delay)
             Motor.setStep(self, 1, 0)
             
             delay += 0.00000125'''  
    def stop(self):
        Motor.setStep(self, 0, 0)
        
def destroy():
    GPIO.cleanup()
    
if __name__ == '__main__':
    
    #motot2 = Motor(dir1, step1)
    #motot3 = Motor(dir2, step2)
    
    motot1 = Motor(dir1, step1)
    motot2 = Motor(dir2, step2)
    motot3 = Motor(dir3, step3)
    motot4 = Motor(dir4, step4)
    motot5 = Motor(dir5, step5)
    motot6 = Motor(dir6, step6)
                
    thread1 = threading.Thread(target=motot1.backward1, args=(motorback_delay, motor1_steps))
    thread2 = threading.Thread(target=motot2.backward2, args=(motorback_delay, motor1_steps))
    thread3 = threading.Thread(target=motot3.backward3, args=(motorback_delay, motor1_steps))
    thread4 = threading.Thread(target=motot4.backward4, args=(motorback_delay, motor1_steps))
    thread5 = threading.Thread(target=motot5.backward5, args=(motorback_delay, motor1_steps))
    thread6 = threading.Thread(target=motot6.backward6, args=(motorback_delay, motor1_steps))
                
    thread1.start()
    thread2.start()
    thread3.start()
    thread4.start()
    thread5.start()
    thread6.start()
    
    thread1.join()
    thread2.join()
    thread3.join()
    thread4.join()
    thread5.join()
    thread6.join()
    
    t1_up = threading.Thread(target=motot1.forward, args=(motorup_delay, motor1_steps))
    t2_up = threading.Thread(target=motot2.forward, args=(motorup_delay, motor1_steps))    
    t3_up = threading.Thread(target=motot3.forward, args=(motorup_delay, motor1_steps))
    t4_up = threading.Thread(target=motot4.forward, args=(motorup_delay, motor1_steps))
    t5_up = threading.Thread(target=motot5.forward, args=(motorup_delay, motor1_steps))
    t6_up = threading.Thread(target=motot6.forward, args=(motorup_delay, motor1_steps))
    
    t1_up.start()
    t2_up.start()
    t3_up.start()
    t4_up.start()
    t5_up.start()
    t6_up.start()
    
    '''motot5 = Motor(dir1, step1)
    thread5 = threading.Thread(target=motot5.backward5, args=(motor1_delay, motor1_steps))
    thread5.start()'''
    

