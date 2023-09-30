import math

def mph2ms(x):
    return x * 0.44704

def in2m(x):
    return x * 0.0254

rho = 1.2 # kg/m^3 density of air
w = in2m(19.2)
h = in2m(5.3)
A = w * h
v = mph2ms(60)

pwind = 0.5 * rho * A * v**3

print("Wind Power at Intake: %.2f W" % pwind)

eff = 0.33

print("Wind Power at 33pct Efficiency: %.2f W" % (pwind*eff))

rad = in2m(8) # m
circ = 2 * math.pi * rad

print("Disc Circumference: %.2f m" % circ)

rps = v / circ
rpm = rps * 60

print("Disc RPS at Equal Tangential Windspeed: %d" % rps)
print("Disc RPM at Equal Tangential Windspeed: %d" % rpm)

pulley1_rad = in2m(3)
pulley1_circ = 2 * math.pi * pulley1_rad

pulley2_rad = in2m(1)
pulley2_circ = 2 * math.pi * pulley2_rad

pulley1_rpm = rpm
pulley2_rpm = pulley1_rpm * pulley1_circ/pulley2_circ

print("RPM at Motor Shaft: %d" % pulley2_rpm)
