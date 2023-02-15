# https://www.govinfo.gov/content/pkg/CFR-2022-title14-vol4/pdf/CFR-2022-title14-vol4-part420-appC.pdf


from math import sqrt, pi, exp

π = pi


# ––––––––––––––––––––––––––––––––––––––––––––––––––– #
# TABLE C–2—IIP RANGE RATE VS. IIP RANGE              #
# ––––––––––––––––––––––––––––––––––––––––––––––––––– #
#
#    IIP RANGE (NM)    |     IIP RANGE RATE (NM/S)
# ––––––––––––––––––––––––––––––––––––––––––––––––––– #
#    0-75              |     0.75
#    76-300            |     1.73
#    301-900           |     4.25
#    901-1700          |     8.85
#    1701-2600         |     19.75
#    2601-3500         |     42.45
#    3501-4500         |     84.85
#    4501-5250         |     154.95






# (Equation C1)

# closest and farthest downrange distance (nm) along the flight corridor centerline to the populated area (see figure C–1)
x1 = 61.1
x2 = 61.8

# closest and farthest cross range distance (nm) to the populated area measured from the flight corridor centerline (see figure C–1)
y1 = 40.7
y2 = 41.2

# one-third of the cross range distance from the centerline to the flight corridor boundary (see figure C–1)
σy = 44.5/3

# Probability of failure = 0.10
Pf = 0.1

# IIP range rate (nm/sec) (see table C–2)
R = 0.75

# C = 643 seconds (constant)
C = 643

q1 = (abs(y2-y1)/σy)/(6*sqrt(2*π))
q2 = exp( -(y1/σy)**2/2 ) + 4 * exp( -((y1+y2)/(2*σy))**2/2 ) + exp( -(y1/σy)**2/2 )
q3 = (Pf/C) * ((x2-x1)/R)
Pi =  q1 * q2 * q3

print("Probability of Impact: %s" % str(Pi))

Ac = 3.14 * 10e-2
Ak = 2.5
Nk = 98

Eck = Pi * (Ac/Ak) * Nk

print("Casualty Expectancy: %s" % str(Eck))

# ×
# Δ
# ℎ
# Φ
# λ
# φ
# α
# ε
# Ω
# β
# Ø
# Σ
# σ
# π
# δ
