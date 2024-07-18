extends Node2D

""" Summary
	-------
This node simulates the displacement of points flexible structure
defined by the parameters below. Continuous time version of
the model can be written as
dx/dt = Ax(t)+Bu(t)      	(1)
y(t) = Cx(t)				(2)
u: input (scalar)
y: output (vector)
x: state vector

When an impulse is applied to the flexible structure, it starts
to vibrate meaning that the points on the structure starts to oscillate.
Some frequencies are dominant in this oscillatory behavior. These dominant
frequencies are called natural modes of the system and depends on the
system parameters.

These natural modes are embedded in the state matrix A. One of the modes can be
represented as
A_i = [
	[  0,     om_i ],
	[ -om_i, -2*zeta_i*om_i ]
]
om_i: frequency in radians
zeta_i: damping cofficient of the natural mode

If zeta is larger, oscillatory behavior vanishes more quickly.

The state matrix A, is block diagonal composition of N modes.
A = diag{ A_1, A_2, ... A_N } which is an 2Nx2N matrix.
Each mode is excited by input u(t), which leads to the 1st equation above.

The amount of displacement of equidistant points on the structure is a
linear combination of natural modes. It is represented by 2nd equation.

The system is discretized in time to iteratively calculate y(t)
in `_process` function. The discrete time version can be written as
x[n+1] = Ax[n]+Bu[n], y[n] = Cx[n]
The system is discretized by using zero order hold method.

The system matrices are generated once when `_ready` function is called.
The calculation steps are given in `_generate_system_matrices`.
Simulation steps are evaluated in `_process` function.
Oscillations are plotted in `_draw` function.
"""

const Lb: float = 20 # Beam length
const Eb: float = 1.0 # Young modulus
const Ib: float = 1.0 # Inertia
const rhob: float = 1.0 # density
const Sb: float = 20.0 # Cross sectional area
const Nb: int = 256 # Number of natural mods
const Ts: float = 1.0 / 40.0 # Sampling period
const zeta: float = 0.01 # Damping
const u1_pos: float = 0.5 * Lb # Position of input on the beam
const u1_mag: float = 300 # Magnitude of input on the beam
const u1_period: float = 30 # Frequency of input on the beam

const NUM_POINTS: int = Nb
const FIRST_POINT: Vector2 = Vector2(100, 300)
const LAST_POINT: Vector2 = Vector2(1000, 300)

var sysA = null
var sysB = null
var sysC = null
var sysx = null
var sysu = null
var show_hand: bool = false

@onready var beam: Path2D = get_node("%Path2D")
@onready var hand: Sprite2D = get_node("%Sprite2D")
var curve: Curve2D = Curve2D.new()
var curve_points: Array = []

func _resonance_freq(i):
	var c1 = pow(i * PI, 2)
	var c2 = sqrt(Eb * Ib / (rhob * Sb * pow(Lb, 4)))

	return c1 * c2

func _generate_system_matrix_block(gbl, om, p_zeta):
	var Ai_cont = gbl.new_mat()
	Ai_cont.from_array([ [ 0, om ], [ -om, -2*p_zeta*om ] ])
	var U = gbl.new_complex_mat(2)
	var V = gbl.new_complex_mat(2)
	var c1 = gbl.new_complex_mat(1)
	c1.set(p_zeta - 1, 0, 0)
	var c2 = gbl.new_complex_mat(1)
	c2.set(p_zeta + 1, 0, 0)
	c1.mul(c2)
	c1.sqrt()
	c2 = c1.copy()
	c1.add(-p_zeta)
	c2.add(p_zeta)
	c2.mul(-1)

	U.set(c1.get(0, 0), 0, 0)
	U.set(c2.get(0, 0), 0, 1)
	U.set(1, 1, 0)
	U.set(1, 1, 1)
	c1.mul(Ts * om)
	c2.mul(Ts * om)
	V.set(c2.get(0, 0), 0, 0)
	V.set(c1.get(0, 0), 1, 1)
	V.exp()
	V.set(0, 0, 1)
	V.set(0, 1, 0)

	var Uinv = U.inv()
	var Ai = U.prod(V).prod(Uinv)

	return [ Ai, Ai_cont ]

func _psi(i, s = null):
	var c1 = sqrt(2.0 / Lb)
	if s == null:
		s = i * Lb / Nb
	var c2 = sin(PI * s / Lb)

	return c1 * c2

func _generate_initial_state(gbl):
	sysx = gbl.new_mat(2 * Nb, 1)

func _generate_system_matrices():

	var gbl = GDBlas.new()
	var Acont = gbl.new_mat(2 * Nb) # Continuous time state matrix
	var A = gbl.new_mat(2 * Nb) # Discrete time state matrix
	var C = gbl.new_mat(Nb, 2 * Nb) # Output matrix is not effected by discretization
	var B = gbl.new_mat(2 * Nb, 1) # Continuous time input matrix
	for i in range(1, Nb + 1):
		var om_i = _resonance_freq(i)
		var Ai_list = _generate_system_matrix_block(gbl, om_i, zeta)
		A.set(Ai_list[0].real(), 2 * (i - 1), 2 * (i - 1)) # Set diagonal blocks
		Acont.set(Ai_list[1], 2 * (i - 1), 2 * (i - 1))

		for j in range(1, Nb + 1):
			C.set(_psi(j) * om_i, j - 1, 2 * (i - 1) + 1)

		B.set(_psi(1, u1_pos) / om_i, 2 * (i - 1) + 1, 0)

	var I = gbl.new_mat(2 * Nb)
	I.eye(-1.0)
	I.add(A)

	sysA = A
	sysC = C
	sysB = Acont.inv().prod(I).prod(B) # Discrete time input matrix
	sysu = gbl.new_mat(1)

	_generate_initial_state(gbl)

func _system_output_to_curve(y):
	curve.clear_points()
	curve.add_point(FIRST_POINT)
	for i in range(1, Nb + 1):
		var yi = y.get(i - 1, 0)
		curve_points[i].y = yi + FIRST_POINT.y
		if i == floor(Nb * u1_pos / Lb):
			hand.position = curve_points[i]

	_create_beam(curve_points)
	beam.curve = curve

func _cubic_bezier(p0: Vector2, p1: Vector2, p2: Vector2, p3: Vector2, t: float):
	var q0 = p0.lerp(p1, t)
	var q1 = p1.lerp(p2, t)
	var q2 = p2.lerp(p3, t)

	var r0 = q0.lerp(q1, t)
	var r1 = q1.lerp(q2, t)

	var s = r0.lerp(r1, t)
	return s

func _create_segment(curve: Curve2D, q0: Vector2, q1: Vector2, n: int = 5):
	var t = 0;
	var dt = 1.0 / float(n);
	var m = 0.5 * (q1 + q0)
	var c1 = Vector2(m.x, q0.y)
	var c2 = Vector2(m.x, q1.y)

	while t < 1:
		var q = _cubic_bezier(q0, c1, c2, q1, t)
		curve.add_point(q)
		t += dt

func _create_beam(points):
	for i in range(len(points) - 1):
		_create_segment(curve, points[i], points[i+1])

func init_path():
	beam.curve = curve

	var dp = (LAST_POINT - FIRST_POINT) / NUM_POINTS
	var point = FIRST_POINT
	for n in range(NUM_POINTS + 1):
		curve_points.append(point)
		point += dp

	_create_beam(curve_points)

# Called when the node enters the scene tree for the first time.
func _ready():
	init_path()
	_generate_system_matrices()

func _draw():
	hand.visible = show_hand
	draw_set_transform(Vector2(550, 300), 0, Vector2(1.0, 0.5))
	draw_circle(Vector2(0, 0), 500, Color.SEA_GREEN)
	draw_set_transform(Vector2(0, 0), 0, Vector2(1.0, 1.0))
	_system_output_to_curve(sysC.prod(sysx))
	draw_polyline(beam.curve.get_baked_points(), Color.TAN, Sb, true)

# Called every frame. 'delta' is the elapsed time since the previous frame.
var accum_delta: float = 0
func _process(delta):
	var u1 = 0
	show_hand = false

	if accum_delta < u1_period / 30:
		u1 = u1_mag * accum_delta
		show_hand = true
	elif accum_delta > u1_period:
		accum_delta = 0.0

	sysu.set(u1, 0, 0)
	sysx = sysA.prod(sysx)
	sysx.add(sysB.prod(sysu))
	queue_redraw()

	accum_delta += delta
