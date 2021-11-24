import numpy as np
import math

def length(a1,b1,):
#上平台初始位置
    p = np.array([0,0,440.4])
#姿态参数
    sx, sy, sz = 1, 0, 0#旋转轴
    theta = math.pi/18#旋转角

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
    l = (p_p+r_b1_2+a1_a1-2*p_a1+2*p_rb1-2*rb1_a1)**0.5#length of rod
    return l

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

l1 = length(a1,b1)
l2 = length(a2,b2)
l3 = length(a3,b3)
l4 = length(a4,b4)
l5 = length(a5,b5)
l6 = length(a6,b6)

print(l1, l2, l3, l4, l5, l6)