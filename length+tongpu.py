import numpy as np
import math

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

key1 = 27
key2 = 22
key3 = 25
key4 = 13
key5 = 4
key6 = 17

def length(a,b,c,tet):
#上平台初始位置
    p = np.array([0,0,440.4])
#姿态参数
    sx, sy, sz = a, b, c#旋转轴
    theta = tet#旋转角
    
    #上平台铰链位置
    b1 = np.array([73.20508076,-73.20508076,-36.5])
    b2 = np.array([-73.20508076,-73.20508076,-36.5])
    b3 = np.array([-100,-26.79491924,-36.5])
    b4 = np.array([-26.79491924,100,-36.5])
    b5 = np.array([26.79491924,100,-36.5])
    b6 = np.array([100,-26.79491924,-36.5])


    #下平台铰链位置
    a1 = np.array([45,-167.94228634,21.5])
    a2 = np.array([-45,-167.94228634,21.5])
    a3 = np.array([-167.94228634,45,21.5])
    a4 = np.array([-125.94229,125.94229,21.5])
    a5 = np.array([125.94229,125.94229,21.5])
    a6 = np.array([167.94228634,45,21.5])
    
    def ll(a1,b1):
        sth = math.sin(theta)
        cth = math.cos(theta)
        vth = 1 - cth
        r_11 = sx*sx*vth+cth
        r_12 = sx*sy*vth-sz*sth
        r_13 = sx*sz*vth+sy*sth
        r_21 = sy*sx*vth+sz*sth
        r_22 = sy*sy*vth+cth
        r_23 = sy*sz*vth-sx*sth
        r_31 = sz*sx*vth-sy*sth
        r_32 = sz*sy*vth+sx*sth
        r_33 = sz*sz*vth+cth
        r = np.array([[r_11, r_12, r_13],[r_21, r_22, r_23],[r_31, r_32, r_33]])#rotation matrix

        p_p = np.dot(p,p)

        r_b1 = np.dot(r,b1)
        r_b1_2 = np.dot(r_b1,r_b1)
        a1_a1 = np.dot(a1,a1)
        p_a1 = np.dot(p,a1)
        p_rb1 = np.dot(p,r_b1)
        rb1_a1 = np.dot(r_b1,a1)
        l= (p_p+r_b1_2+a1_a1-2*p_a1+2*p_rb1-2*rb1_a1)**0.5#length of rod-1
        return l
    #qiu gan length
    l_1 = ll(a1,b1)
    l_2 = ll(a2,b2)
    l_3 = ll(a3,b3)
    l_4 = ll(a4,b4)
    l_5 = ll(a5,b5)
    l_6 = ll(a6,b6)
    l_l=[l_1,l_2,l_3,l_4,l_5,l_6]
    return l_l

class Motor:
    def __init__(self, d1, s1):
        self.d1= d1
        self.s1 = s1
        
        GPIO.setwarnings(False)
        GPIO.setmode(GPIO.BCM)
        GPIO.setup(self.d1, GPIO.OUT)#define direction pin output
        GPIO.setup(self.s1, GPIO.OUT)#define pulse pin output

        GPIO.setup(key1, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key2, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key3, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key4, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key5, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        GPIO.setup(key6, GPIO.IN, pull_up_down=GPIO.PUD_UP)
        
    def setStep(self, w1, w2):
        GPIO.output(self.d1, w1)
        GPIO.output(self.s1, w2)
        
    def forward(self, delay, steps):
        for i in range(0, steps):
            Motor.setStep(self, 1, 1)
            time.sleep(delay)
            Motor.setStep(self, 1, 0)

    #translation or rotation backward
    def backward(self, delay, steps):
        for i in range(0, steps):
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)
    
    #fuwei-backward
    def backward1(self, delay, steps):
        while True:
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
           
            if GPIO.input(key1) == 0:
                Motor.forward(self, 0.005, 800)
                break              
    
    def backward2(self, delay, steps):
        while True:
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
       
            if GPIO.input(key2) == 0:
                Motor.forward(self, 0.005, 800)
                break

    def backward3(self, delay, steps):
        while True:
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
       
            if GPIO.input(key3) == 0:
                Motor.forward(self, 0.005, 800)
                break

    def backward4(self, delay, steps):
        while True:
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
       
            if GPIO.input(key4) == 0:
                Motor.forward(self, 0.005, 800)
                break
            
    def backward5(self, delay, steps):
        while True:
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
       
            if GPIO.input(key5) == 0:
                Motor.forward(self, 0.005, 800)
                break
            
    def backward6(self, delay, steps):
        while True:
            Motor.setStep(self, 0, 1)
            time.sleep(delay)
            Motor.setStep(self, 0, 0)            
       
            if GPIO.input(key6) == 0:
                Motor.forward(self, 0.005, 800)
                break            
            
    def stop(self):
        Motor.setStep(self, 0, 0)            

       
def destroy():
    GPIO.cleanup()



motor1_delay = 0.003
motorback_delay = 0.001
motorup_delay = 0.002

fuwei_steps = 400*22
if __name__ == '__main__':
    motor1 = Motor(dir1, step1)
    motor2 = Motor(dir2, step2)
    motor3 = Motor(dir3, step3)
    motor4 = Motor(dir4, step4)
    motor5 = Motor(dir5, step5)
    motor6 = Motor(dir6, step6)
    a = int(input('please input x:'))
    b = int(input('please input y:'))
    c = int(input('please input z:'))
    d = int(input('please input theta:'))
    
    l=length(a,b,c,d*math.pi/180)
    
    l1=l[0]
    l2=l[1]
    l3=l[2]
    l4=l[3]
    l5=l[4]
    l6=l[5]
    #position limit
    if (l1<=345 or l1>=445 or l2<=345 or l2>=445 or l3<=345 or l3>=445 or l4<=345 or l4>=445 or l5<=345 or l5>=445 or l6<=345 or l6>=445):
        print(l1,l2,l3,l4,l5,l6)
        print('length is beyong the limit')
    #decide motor1 is turning forward or backward
    else:
        if l1 <= 395:
            motor1_steps = int((395-l1)*200)
            target1 = motor1.backward
        else:
            motor1_steps = int((l1-395)*200)
            target1 = motor1.forward
    
    #decide motor2 is turning forward or backward
        if l2 <= 395:
            motor2_steps = int((395-l2)*200)
            target2 = motor2.backward
        else:
            motor2_steps  = int((l2-395)*200)
            target2 = motor2.forward
    
    #decide motor3 is turning forward or backward
        if l3 <= 395:
            motor3_steps = int((395-l3)*200)
            target3 = motor3.backward
        else:
            motor3_steps = int((l3-395)*200)
            target3 = motor3.forward

    #decide motor4 is turning forward or backward
        if l4 <= 395:
            motor4_steps = int((395-l4)*200)
            target4 = motor4.backward
        else:
            motor4_steps = int((l4-395)*200)
            target4 = motor4.forward

    #decide motor5 is turning forward or backward
        if l5 <= 395:            
            motor5_steps = int((395-l5)*200)
            target5 = motor5.backward
        else:
            motor5_steps = int((l5-395)*200)
            target5 = motor5.forward

    #decide motor6 is turning forward or backward
        if l6 <= 395:
            motor6_steps = int((395-l6)*200)
            target6 = motor6.backward
        else:
            motor6_steps = int((l6-395)*200)
            target6 = motor6.forward

    try:
        while True:
            shuru = int(input('input:'))
            if shuru == 1:#translation or rotation function
                thread1 = threading.Thread(target=target1, args=(motor1_delay, motor1_steps))
                thread2 = threading.Thread(target=target2, args=(motor1_delay, motor2_steps))
                thread3 = threading.Thread(target=target3, args=(motor1_delay, motor3_steps))
                thread4 = threading.Thread(target=target4, args=(motor1_delay, motor4_steps))
                thread5 = threading.Thread(target=target5, args=(motor1_delay, motor5_steps))
                thread6 = threading.Thread(target=target6, args=(motor1_delay, motor6_steps))
                
                thread1.start()
                thread2.start()
                thread3.start()
                thread4.start()
                thread5.start()
                thread6.start()
            if shuru == 2:#fuwei-function
                #print(2)
                thread1 = threading.Thread(target=motor1.backward1, args=(motorback_delay, motor1_steps))
                thread2 = threading.Thread(target=motor2.backward2, args=(motorback_delay, motor1_steps))
                thread3 = threading.Thread(target=motor3.backward3, args=(motorback_delay, motor1_steps))
                thread4 = threading.Thread(target=motor4.backward4, args=(motorback_delay, motor1_steps))
                thread5 = threading.Thread(target=motor5.backward5, args=(motorback_delay, motor1_steps))
                thread6 = threading.Thread(target=motor6.backward6, args=(motorback_delay, motor1_steps))
                
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
    
                t1_up = threading.Thread(target=motor1.forward, args=(motorup_delay, fuwei_steps))
                t2_up = threading.Thread(target=motor2.forward, args=(motorup_delay, fuwei_steps))    
                t3_up = threading.Thread(target=motor3.forward, args=(motorup_delay, fuwei_steps))
                t4_up = threading.Thread(target=motor4.forward, args=(motorup_delay, fuwei_steps))
                t5_up = threading.Thread(target=motor5.forward, args=(motorup_delay, fuwei_steps))
                t6_up = threading.Thread(target=motor6.forward, args=(motorup_delay, fuwei_steps))
    
                t1_up.start()
                t2_up.start()
                t3_up.start()
                t4_up.start()
                t5_up.start()
                t6_up.start()
    except KeyboardInterrupt:
        destroy()