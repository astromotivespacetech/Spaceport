# Testing 

# https://www.govinfo.gov/content/pkg/CFR-2012-title14-vol4/pdf/CFR-2012-title14-vol4-part420-appA.pdf


# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Table 1 of § 420.19  Orbital Expendable Launch Vehicle Classes by Payload Weight (lbs)
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
#
# 100 nm orbit           |	                             Weight class
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
#                              Small       |      Medium      |     Medium large     |	     Large
#
# 28 degrees inclination *	   ≤4400       |  >4400 to ≤11100 |	  >11100 to ≤18500	 |      >18500
#
# 90 degrees inclination	   ≤3300	   |  >3300 to ≤8400  |	  >8400 to ≤15000	 |      >15000




# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# TABLE A–1—DEBRIS DISPERSION RADIUS (Dmax) (IN)                             #
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
#                                                                            #
#                    ORBITAL LAUNCH VEHICLES               |  SUBORBITAL LV  #
#                                                                            #
#   |  SMALL  |   MEDIUM   |   MEDIUM LARGE   |   LARGE    |  GUIDED         #
#                                                                            #
#      87,600     111,600      127,200            156,000     96,000         #
#      (1.2 nm)   (1.53 nm)    (1.74 nm)          (2.14 nm)   (1.32 nm)      #

DEBRIS_DISPERSION_RADUIS = {
    'SMALL': 1.2,
    'MEDIUM': 1.53,
    'MEDIUM-LARGE': 1.74,
    'LARGE': 2.14,
    'GUIDED': 1.32
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# TABLE A–2—OVERFLIGHT EXCLUSION ZONE DOWNRANGE DISTANCE (Doez) (IN)         #
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
#                                                                            #
#                    ORBITAL LAUNCH VEHICLES               |  SUBORBITAL LV  #
#                                                                            #
#   |  SMALL  |   MEDIUM   |   MEDIUM-LARGE   |   LARGE    |  GUIDED         #
#                                                                            #
#      240,500    253,000      310,300            937,700     232,100        #
#      (3.3 nm)   (3.47 nm)    (4.26 nm)         (12.86 nm)   (3.18 nm)      #


OVERFLIGHT_EXCLUSION_ZONE_DOWNRANGE_DISTANCE = {
    'SMALL': 3.3,
    'MEDIUM': 3.47,
    'MEDIUM-LARGE': 4.26,
    'LARGE': 12.86,
    'GUIDED': 3.18
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# TABLE A-3: FLIGHT CORRIDOR LINE SEGMENTS LENGTH                            #
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
#                                                                            #
# Dmax(in)                  | LINE SEGMENT LENGTHS (x 10^6 inches)           #
#                           |                                                #
# ORBITAL LAUNCH VEHICLES   |      CF      |      DE      |      HI      |   #
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# SMALL         | 87,600    |    2.87620   |    8.59452   |    128.566   |   #
#                               (39.45 nm)    (117.87 nm)   (1,763.27 nm)
# MEDIUM        | 111,600   |    2.97220   |    8.64252   |    128.566   |   #
#                               (40.76 nm)    (118.53 nm)   (1,763.27 nm)
# MEDIUM-LARGE  | 127,200   |    3.03460   |    8.67372   |    128.566   |   #
#                               (41.62 nm)    (118.96 nm)   (1,763.27 nm)
# LARGE         | 156,000   |    3.14979   |    8.73131   |    128.566   |   #
#                               (43.20 nm)    (119.75 nm)   (1,763.27 nm)
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# SUBORBITAL LV | 96,000    |    2.90980   |    8.61132   |    N/A       |   #
#                               (39.91 nm)    (118.10 nm)                    #
#                                                                            #


FLIGHT_CORRIDOR_LINE_SEGMENTS_LENGTH = {
    'SMALL': {
        'CF': 39.45,
        'DE': 117.87,
        'HI': 1763.27
    },
    'MEDIUM': {
        'CF': 40.76,
        'DE': 118.53,
        'HI': 1763.27
    },
    'MEDIUM-LARGE': {
        'CF': 41.62,
        'DE': 118.96,
        'HI': 1763.27
    },
    'LARGE': {
        'CF': 43.2,
        'DE': 119.75,
        'HI': 1763.27
    },
    'GUIDED': {
        'CF': 39.91,
        'DE': 118.1,
        'HI': None
    }
}


from math import pi, sin, cos, asin, acos, tan, atan, degrees, radians
π = pi

# Appendix A to Part 420 - Method for Defining a Flight Corridor

# (a) Introduction

# (1) This appendix provides a method for constructing a flight corridor from a launch point for a guided suborbital launch vehicle
# or any one of the four classes of guided orbital launch vehicles from table 1, § 420.19, without the use of local meteorological data
# or a launch vehicle trajectory.

# (2) A flight corridor includes an overflight exclusion zone in a launch area and, for a guided suborbital launch vehicle, an impact
# dispersion area in a downrange area. A flight corridor for a guided suborbital launch vehicle ends with the impact dispersion area,
# and, for the four classes of guided orbital launch vehicles, 5000 nautical miles (nm) from the launch point.

# (b) Data requirements

# (1) Maps. An applicant shall use any map for the launch site region with a scale not less than 1:250,000 inches per inch in the launch
# area and 1:20,000,000 inches per inch in the downrange area. As described in paragraph (b)(2), an applicant shall use a mechanical method,
# a semi-automated method, or a fully-automated method to plot a flight corridor on maps. A source for paper maps acceptable to the FAA is
# the U.S. Dept. of Commerce, National Oceanic and Atmospheric Administration, National Ocean Service.

# (i) Projections for mechanical plotting method. An applicant shall use a conic projection. The FAA will accept a “Lambert-Conformal” conic
# projection. A polar aspect of a plane-azimuthal projection may also be used for far northern launch sites.

# (ii) Projections for semi-automated plotting method. An applicant shall use cylindrical, conic, or plane projections for semi-automated
# plotting. The FAA will accept “Mercator” and “Oblique Mercator” cylindrical projections. The FAA will accept “Lambert-Conformal” and
# “Albers Equal-Area” conic projections.The FAA will accept “Lambert Azimuthal Equal-Area” and “Azimuthal Equidistant” plane projections.

# (iii) Projections for fully-automated plotting method. The FAA will accept map projections used by geographical information system software
# scaleable pursuant to the requirements of paragraph (b)(1).

# (2) Plotting Methods.

# (i) Mechanical method. An applicant may use mechanical drafting equipment such as pencil, straight edge, ruler, protractor, and compass to
# plot the location of a flight corridor on a map. The FAA will accept straight lines for distances less than or equal to 7.5 times the map
# scale on map scales greater than or equal to 1:1,000,000 inches per inch (in/in); or straight lines representing 100 nm or less on map scales
# less than 1:1,000,000 in/in.

# (ii) Semi-automated method. An applicant may employ the range and bearing techniques in paragraph (b)(3) to create latitude and longitude
# points on a map. The FAA will accept straight lines for distances less than or equal to 7.5 times the map scale on map scales greater than or
# equal to 1:1,000,000 inches per inch (in/in); or straight lines representing 100 nm or less on map scales less than 1:1,000,000 in/in.

# (iii) Fully-automated method. An applicant may use geographical information system software with global mapping data scaleable in accordance
# with paragraph (b)(1).

# (3) Range and bearing computations on an ellipsoidal Earth model.

# (i) To create latitude and longitude pairs on an ellipsoidal Earth model, an applicant shall use the following equations to calculate geodetic
# latitude (+N) and longitude (+E) given the launch point geodetic latitude (+N), longitude (+E), range (nm), and bearing (degrees, positive
# clockwise from North).

# (A) Input. An applicant shall use the following input in making range and bearing computations. Angle units must be in radians.

# Φ1 = Geodetic latitude of launch point (radians) = Φ1 (DDD)·π/180 (radians per degree)

# λ1 = Longitude of launch point (DDD) = λ(DDD)·π/180 (radians per degree)

# S = Range from launch point (nm) = S(DDD)·π/180 (radians per degree)

# α12 = Azimuth bearing from launch point (deg) = α12(DDD)·π/180 (radians per degree)


#(B) Computations. An applicant shall use the following equations to determine the latitude (φ2) and longitude (λ2) of a target point situated “S”
# nm from the launch point on an azimuth bearing (α12) degrees.


# where:
# a = WGS-84 semi-major axis (3443.91846652 nmi)
# b = WGS-84 semi-minor axis (3432.37165994 nmi)
a = 3443.91846652
b = 3432.37165994


def lat_lon_of_target(Φ, λ, α, s):

    # print(Φ, λ, α, s)

    # (Equation A1)
    f = 1 - b/a

    # (Equation A2)
    ε2 = (a**2-b**2)/b**2

    # (Equation A3)
    Ø = s/b

    # (Equation A4)
    β1 = atan( (b*sin(Φ)) / (a*cos(Φ)) )

    # (Equation A5)
    g = cos(β1)*cos(α)

    # (Equation A6)
    h = cos(β1)*sin(α)

    # (Equation A7)
    m = ( (1 + (ε2/2) * sin(β1)**2) * (1 - h**2) ) / 2

    # (Equation A8)
    n = ( (1 + (ε2/2) * sin(β1)**2) * (sin(β1)**2 * cos(Ø) + g * sin(β1) * sin(Ø)) ) / 2

    # (Equation A9)
    L = h * (-f * Ø + 3 * f**2 * n * sin(Ø) + ((3 * f**2 * m * (Ø-sin(Ø)*cos(Ø))) / 2)  )

    # (Equation A10)
    M = m * ε2

    # (Equation A11)
    N = n * ε2

    # (Equation A12)
    A1 = N * sin(Ø)

    # (Equation A13)
    A2 = (M/2) * (sin(Ø) * cos(Ø) - Ø)

    # (Equation A14)
    A3 = (5/2) * (N**2 * sin(Ø) * cos(Ø))

    # (Equation A15)
    A4 = (M**2/16) * (11 * Ø - 13 * sin(Ø) * cos(Ø) - 8 * Ø * cos(Ø)**2 + 10 * sin(Ø) * cos(Ø)**3)

    # (Equation A16)
    A5 = (M*N/2) * (3 * sin(Ø) + 2 * Ø * cos(Ø) - 5 * sin(Ø) * cos(Ø)**2)

    # (Equation A17)
    δ = Ø - A1 + A2 + A3 + A4 + A5

    # (Equation A18)
    sinβ2 = sin(β1) * cos(δ) + g * sin(δ)

    # Equation A19)
    cosβ2 = (h**2 + (g * cos(δ) - sin(β1) * sin(δ))**2 )**0.5

    # (Equation A20) (geodetic latitude of target point, DDD)
    Φ2 = atan( (a*sinβ2)/(b*cosβ2) ) * 180/π

    # (Equation A21)
    Δ = atan( (sin(δ)*sin(α))/(cos(β1)*cos(δ)-sin(β1)*sin(δ)*cos(α)) )

    # (Equation A22) (longitude of target point,DDD)
    λ2 = (λ+Δ+L) * 180/π

    # print("Geodetic Latitde of Target Point: %s deg" % str(Φ2))
    # print("Geodetic Longitude of Target Point: %s deg" % str(λ2))

    return Φ2, λ2





# (ii) To create latitude and longitude pairs on an ellipsoidal Earth model, an applicant shall use the following equations to
# calculate the distance (S) of the geodesic between two points (P1 and P2), the forward azimuth (a12) of the geodesic at P1, and the
# back azimuth (a21) of the geodesic at P2, given the geodetic latitude (+N), longitude (+E) of P1 and P2. Azimuth is measured
# positively clockwise from North. (A) Input. An applicant shall use the following input. Units must be in radians.

# (B) Computations. An applicant shall use the following equations to determine the distance (S), the forward azimuth (a12) of the
# geodesic at P1, and the back azimuth (a12) of the geodesic at P2.

# f = 1 − b/a (Equation A23)

# where:
# a = WGS-84 semi-major axis (3443.91846652 nmi)
# b = WGS-84 semi-minor axis (3432.37165994 nmi)
# a = 3443.91846652
# b = 3432.37165994
#
# # (Equation A23)
# f = 1 - b/a
#
# # (Equation A24)
# L = λ2 - λ1
#
# # (Equation A25)
# β1 = atan( (b*sin(Φ1)) / (α12*cos(Φ1)) )
#
# # (Equation A26)
# β2 = atan( (b*sin(Φ2)) / (α12*cos(Φ2)) )
#
# # (Equation A27)
# A = sin(β1) * sin(β2)
#
# # (Equation A28)
# B = cos(β1) * cos(β2)
#
# # (Equation A29)
# cosδ = A + B * cos(L)
#
# # (Equation A30)
# n = (a-b)/(a+b)
#
# # (Equation A31)
# β2_β1 = (Φ2-Φ1) + 2 * (A * (n+n**2+n**3) - B * (n-n**2+n**3) ) * sin(Φ2-Φ1)
#
# # (Equation A32)
# sinδ = ( (sin(L)*cos(β2))**2 + ( sin(β2_β1) + 2 * cos(β2) * sin(β1) * sin(L/2)**2 )**2  )**0.5
#
# # (Equation A33)
# δ = atan(sinδ/cosδ) # evaluated in positive radians <= π
#
# # (Equation A34)
# c = B * sin(L) / sinδ
#
# # (Equation A35)
# m = 1-c*2
#
# # (Equation A36)
# one = δ * (1+f+f**2) + A * ( (f+f**2) * sin(δ) - (f**2 * δ**2)/(2*sin(δ)) )
# two = -(m/2) * ( (f+f**2) * (δ+sin(δ)*cos(δ)) - (f**2 * δ**2)/(tan(δ)) )
# three = -(A**2 * f**2/2) * sin(δ) * cos(δ)
# four = (f**2 * m**2/16) * (δ + sin(δ) * cos(δ) - 2 * sin(δ) * cos(δ)**3 - 8 * δ**2 / (tan(δ)) )
# five = (A**2 * m * f**2/2) * (sin(δ) * cos(δ) + δ + δ**2 / (sin(δ)) )
# S = b * (one + two + three + four + five)
#
# # (Equation A37)
# Λ = L + c * (δ * (f+f**2) - (A*f**2/2) * (sin(δ)+2*δ**2/(sin(δ))) + (m*f**2/4) * (sin(δ)*cos(δ)-5*δ+4*δ**2/(tan(δ)))  )
#
# # (Equation A38)
# α12 = atan( (cos(β2)*sin(Λ)) / (sin(β2_β1)+2*cos(β2)*sin(β1*sin(Λ/2)**2)) ) * (180/π)
#
# # (Equation A39)
# α21 = atan( (-cos(β1)*sin(Λ)) / (2*cos(β1)*sin(β2)*sin(Λ/2)**2-sin(β2_β1)) ) * (180/π)










# (c) Creation of a Flight Corridor
#
# (1) To define a flight corridor, an applicant shall:
#
# (i) Select a guided suborbital or orbital launch vehicle, and, for an orbital launch vehicle, select from table 1 of §420.19 a launch
# vehicle weight class that best represents the launch vehicle the applicant plans to support at its launch point;

LV = 'SMALL'

#
# (ii) Select a debris dispersion radius (Dmax) from table A–1 corresponding to the guided suborbital launch vehicle or orbital launch
# vehicle class selected in paragraph (c)(1)(i);

Dmax = DEBRIS_DISPERSION_RADUIS[LV]

#
# (iii) Select a launch point geodetic latitude and longitude; and
#


# Φ1 = Geodetic latitude of launch point (radians) = Φ1 (DDD)·π/180 (radians per degree)
ΦLP = 43.07961 * π/180

# λ1 = Longitude of launch point (DDD) = λ(DDD)·π/180 (radians per degree)
λLP = -119.94511 * π/180


# (iv) Select a flight azimuth.

# α12 = Azimuth bearing from launch point (deg) = α12(DDD)·π/180 (radians per degree)
αft = 105 * π/180

αdmax = αft + π

ΦDmax, λDmax = lat_lon_of_target(ΦLP, λLP, αdmax, Dmax)

print(ΦDmax, λDmax)


# (2) An applicant shall define and map an overflight exclusion zone using the following method:
#
# (i) Select a debris dispersion radius (Dmax) from table A–1 and a downrange distance (DOEZ) from table A–2 to define an overflight
# exclusion zone for the guided suborbital launch vehicle or orbital launch vehicle class selected in paragraph (c)(1)(i).
#

Doez = OVERFLIGHT_EXCLUSION_ZONE_DOWNRANGE_DISTANCE[LV]

ΦDoez, λDoez = lat_lon_of_target(ΦLP, λLP, αft, Doez)

ΦCF, λCF = lat_lon_of_target(ΦLP, λLP, αft, 10)

αt = αft - π/2
αb = αft + π/2

S = FLIGHT_CORRIDOR_LINE_SEGMENTS_LENGTH[LV]['CF']*0.5

ΦC, λC = lat_lon_of_target(ΦCF * π/180, λCF * π/180, αt, S)
ΦF, λF = lat_lon_of_target(ΦCF * π/180, λCF * π/180, αb, S)


ΦDE, λDE = lat_lon_of_target(ΦLP, λLP, αft, 100)

S = FLIGHT_CORRIDOR_LINE_SEGMENTS_LENGTH[LV]['DE']*0.5

ΦD, λD = lat_lon_of_target(ΦDE * π/180, λDE * π/180, αt, S)
ΦE, λE = lat_lon_of_target(ΦDE * π/180, λDE * π/180, αb, S)


# Flight trajectory
ΦHI, λHI = lat_lon_of_target(ΦLP, λLP, αft, 5000)

print(ΦHI, λHI)

S = FLIGHT_CORRIDOR_LINE_SEGMENTS_LENGTH[LV]['HI']*0.5

ΦH, λH = lat_lon_of_target(ΦHI * π/180, λHI * π/180, αt, S)
ΦI, λI = lat_lon_of_target(ΦHI * π/180, λHI * π/180, αb, S)

print(ΦH, λH)
print(ΦI, λI)



# (ii) An overflight exclusion zone is described by the intersection of the following boundaries, which are depicted in figure A–1:
#
# (A) An applicant shall define an uprange boundary with a half-circle arc of radius Dmax and a chord of length twice Dmax connecting
# the half-circle arc endpoints. The uprange boundary placement on a map has the chord midpoint positioned on the launch point
# with the chord oriented along an azimuth ±90°from the launch azimuth and the halfcircle arc located uprange from the launch point.



# (B) An applicant shall define the downrange boundary with a half-circle arc of radius Dmax and a chord of length twice Dmax connecting
# the half-circle arc endpoints. The downrange boundary placement on a map has the chord midpoint intersecting the nominal flight azimuth
# line at a distance DOEZ inches downrange with the chord oriented along an azimuth ±90°from the launch azimuth and the half-circle arc
# located downrange from the intersection of the chord and the flight azimuth line.
#
# (C) Crossrange boundaries of an overflight exclusion zone are defined by two lines segments. Each is parallel to the flight azimuth
# with one to the left side and one to the right side of the flight azimuth line. Each line connects an uprange half-circle arc endpoint
# to a downrange half-circle arc endpoint as shown in figure A–1.
#
# (iii) An applicant shall identify the overflight exclusion zone on a map that meets the requirements of paragraph (b).
#
# (3) An applicant shall define and map a flight corridor using the following method:
#
# (i) In accordance with paragraph (b), an applicant shall draw a flight corridor on one or more maps with the Dmax origin centered on
# the intended launch point and the flight corridor centerline (in the downrange direction) aligned with the initial flight azimuth. The
# flight corridor is depicted in figure A–2 and its line segment lengths are tabulated in table A–3.
#
# (ii) An applicant shall define the flight corridor using the following boundary definitions:
#
# (A) An applicant shall draw an uprange boundary, which is defined by an arc-line GB (figure A–2), directly uprange from and centered on
# the intended launch point with radius Dmax.
#
# (B) An applicant shall draw line CF perpendicular to and centered on the flight azimuth line, and positioned 10 nm downrange from
# the launch point. The applicant shall use the length of line CF provided in table A–3 corresponding to the guided suborbital launch
# vehicle or orbital launch vehicle class selected in paragraph (c)(1)(i).
#
# (C) An applicant shall draw line DE perpendicular to and centered on the flight azimuth line, and positioned 100 nm downrange from
# the launch point. The applicant shall use the length of line DE provided in table A–3 corresponding to the guided suborbital launch
# vehicle or orbital launch vehicle class selected in paragraph (c)(1)(i).
#
# (D) Except for a guided suborbital launch vehicle, an applicant shall draw a downrange boundary, which is defined by line HI and is
# drawn perpendicular to and centered on the flight azimuth line, and positioned 5,000 nm downrange from the launch point.
# The applicant shall use the length of line HI provided in table A–3 corresponding to the orbital launch vehicle class selected in
# paragraph (c)(1)(i).
#
# (E) An applicant shall draw crossrange boundaries, which are defined by three lines on the left side and three lines on the right side of
# the flight azimuth. An applicant shall construct the left flight corridor boundary according to the following, and as depicted in figure A–3 :
#
# (1) The first line (line BC in figure A–3) is tangent to the uprange boundary arc, and ends at endpoint C of line CF, as depicted in figure A–3;
#
# (2) The second line (line CD in figure A–3) begins at endpoint C of line BC and ends at endpoint D of line DH, as depicted in figure A–3;
#
# (3) For all orbital launch vehicles, the third line (line DH in figure A–3) begins at endpoint D of line CD and ends at endpoint H of line HI,
# as depicted in figure A–3; and
#
# (4) For a guided suborbital launch vehicle, the line DH begins at endpoint D of line CD and ends at a point tangent to the impact dispersion
# area drawn in accordance with paragraph (c)(4) and as depicted in figure A–4.
#
# (F) An applicant shall repeat the procedure in paragraph (c)(3)(ii)(E) for the right side boundary.
#
# (iii) An applicant shall identify the flight corridor on a map that meets the requirements of paragraph (b).
#
# (4) For a guided suborbital launch vehicle, an applicant shall define a final stage impact dispersion area as part of the flight corridor
# and show the impact dispersion area on a map, as depicted in figure A–4, in accordance with the following:
#
# (i) An applicant shall select an apogee altitude (Hap) for the launch vehicle final stage. The apogee altitude should equal the highest
# altitude intended to be reached by a guided suborbital launch vehicle launched from the launch point.
#
# (ii) An applicant shall define the impact dispersion area by using an impact range factor [IP(Hap)] and a dispersion factor [DISP(Hap)] as
# shown below:
#
# (A) An applicant shall calculate the impact range (D) for the final launch vehicle stage. An applicant shall set D equal to the maximum
# apogee altitude (Hap) multiplied by the impact range factor as shown below:

# (Equation A40)
# D = Hap * IP(Hap)
#
# where: IP(Hap) = 0.4 for an apogee less than 100 km; and IP(Hap) = 0.7 for an apogee 100 km or greater.
#
# (B) An applicant shall calculate the impact dispersion radius (R) for the final launch vehicle stage. An applicant shall set R equal to
# the maximum apogee altitude (Hap) multiplied by the dispersion factor as shown below:

# (Equation A41)
# R = Hap * DISP(Hap)

# where: DISP(Hap) = 0.05


# (iii) An applicant shall draw the impact dispersion area on a map with its center on the predicted impact point. An applicant
# shall then draw line DH in accordance with paragraph (c)(3)(ii)(E)(4).
#
# (d) Evaluate the Flight Corridor
#
# (1) An applicant shall evaluate the flight corridor for the presence of any populated areas. If an applicant determines that no
# populated area is located within the flight corridor, then no additional steps are necessary.
#
# (2) If a populated area is located in an overflight exclusion zone, an applicant may modify its proposal or demonstrate that there
# are times when no people are present or that the applicant has an agreement in place to evacuate the public from the overflight
# exclusion zone during a launch.
#
# (3) If a populated area is located within the flight corridor, an applicant may modify its proposal and create another flight corridor
# pursuant to appendix A, use appendix B to narrow the flight corridor, or complete a risk analysis in accordance with appendix C.





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







#
