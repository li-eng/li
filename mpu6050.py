import smbus
import math
import time
import matplotlib.pyplot as plt
import numpy as np

# 电源管理寄存器地址
power_mgmt_1 = 0x6b
power_mgmt_2 = 0x6c

def read_byte(adr):
    return bus.read_byte_data(address, adr)

def read_word(adr):
    high = bus.read_byte_data(address, adr)
    low = bus.read_byte_data(address, adr+1)
    val = (high << 8) + low
    return val

def read_word_2c(adr):
    val = read_word(adr)
    if (val >= 0x8000):
        return -((65535 - val) + 1)
    else:
        return val

def dist(a,b):
    return math.sqrt((a*a)+(b*b))
    #math.sqrt(x) 方法返回数字x的平方根。

def get_y_rotation(x,y,z):
    radians = math.atan2(x, dist(y,z))
    #math.atan2(y, x) 返回给定的 X 及 Y 坐标值的反正切值。
    return -math.degrees(radians)
    #math.degrees(x) 将弧度x转换为角度。

def get_x_rotation(x,y,z):
    radians = math.atan2(y, dist(x,z))
    return math.degrees(radians)

#display data
fig = plt.figure()
ax = fig.add_subplot(1,2,1)
ay = fig.add_subplot(1,2,2)

ax.set_xlabel('Time')
ax.set_ylabel('xrotation')
ax.set_title('')

ay.set_xlabel('Time')
ay.set_ylabel('yrotation')

line = None
plt.grid(True)
plt.ion()
obsX = []
obsY = []

obsY1 = []

bus = smbus.SMBus(1) # or bus = smbus.SMBus(1) for Revision 2 boards
address = 0x68       # This is the address value read via the i2cdetect command

# Now wake the 6050 up as it starts in sleep mode
bus.write_byte_data(address, power_mgmt_1, 0)
t = 0

while True:
    #time.sleep(0.1)
    gyro_xout = read_word_2c(0x43)
    gyro_yout = read_word_2c(0x45)
    gyro_zout = read_word_2c(0x47)
    
    #print ' '
    #print ("gyro_xout : ", (gyro_xout / 131))#倍率：±250°/s
    #print ("gyro_yout : ", (gyro_yout / 131))
    #print ("gyro_zout : ", (gyro_zout / 131))

    accel_xout = read_word_2c(0x3b)
    accel_yout = read_word_2c(0x3d)
    accel_zout = read_word_2c(0x3f)

    accel_xout_scaled = accel_xout / 16384.0 #倍率：±2g
    accel_yout_scaled = accel_yout / 16384.0
    accel_zout_scaled = accel_zout / 16384.0

    #print ("accel_xout scaled: ", accel_xout_scaled)
    #print ("accel_yout scaled: ", accel_yout_scaled)
    #print ("accel_zout scaled: ", accel_zout_scaled)
    x_ro = get_x_rotation(accel_xout_scaled, accel_yout_scaled, accel_zout_scaled)-1.48
    y_ro = get_y_rotation(accel_xout_scaled, accel_yout_scaled, accel_zout_scaled)+1.97
    
    #print ("x rotation: " , x_ro)
    #print ("y rotation: " , y_ro)
    
    t = t + 0.2
    obsX.append(t)
    obsY.append(x_ro)
    obsY1.append(y_ro) 
    print('x_ro:',x_ro)
    print('y_ro:',y_ro)
    #print(np.mean(obsY))
    
    #if line is None:
    line = ax.plot(obsX, obsY, '-g', marker = '*')[0]
    #line1 = ay.plot(obsX, obsY1, '-g', marker = '*')[0]
    
    line.set_xdata(obsX)
    line.set_ydata(obsY)
    
    #line1.set_xdata(obsX)
    #line1.set_ydata(obsY1)    

    
    ax.set_xlim([t-20,t+10])
    ax.set_ylim([-2,2])
    
    ay.set_xlim([t-20,t+10])
    ay.set_ylim([-2,2])
       
    
    plt.pause(0.2)
