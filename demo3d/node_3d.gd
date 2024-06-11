extends Node3D

const CAPTURE_PERIOD: float = 1.0 / 3e1
const DRAW_PERIOD: float = 1.0 / 3e1
var GBL: GDBlas = GDBlas.new()
var R = GBL.new_mat()
var G = GBL.new_mat()
var B = GBL.new_mat()
var _img_size: Vector2i = Vector2i()
var _edge_filter = null # Edge detection filter
var _motion_filt = null # Motion blur filter
var _filtered_mat1 = null
var _filtered_mat2 = null
@onready var filtered_view1: Sprite2D = get_node("%Sprite2D")
@onready var filtered_view2: Sprite2D = get_node("%Sprite2D2")

func _capture_and_process():
	var vp: Viewport = get_viewport()
	var img: Image = vp.get_texture().get_image()
	var img_size: Vector2i = img.get_size() / 3
	img.resize(img_size.x, img_size.y)
	if img_size != _img_size:
		R.resize(img_size.y, img_size.x)
		G.resize(img_size.y, img_size.x)
		B.resize(img_size.y, img_size.x)
		_img_size = img_size

	var format = img.get_format()
	if format == Image.FORMAT_RGB8:
		var img_data: PackedByteArray = img.get_data()
		# Get image channels into GDBlasMat matrix
		R.unpack(img_data, GDBlas.REAL_COMPONENT, 3, 0)
		G.unpack(img_data, GDBlas.REAL_COMPONENT, 3, 1)
		B.unpack(img_data, GDBlas.REAL_COMPONENT, 3, 2)
		# Convert to gray scale
		R.mul(0.299)
		G.mul(0.587)
		B.mul(0.114)
		R.add(G)
		R.add(B)
		_filtered_mat1 = R.conv(_edge_filter, true)
		_filtered_mat2 = R.conv(_motion_filt, true)
	else:
		printerr("Image must be in RGB8 format")

func _draw_filtered_viewport():
	if _filtered_mat1 == null or _filtered_mat2 == null:
		return

	var vp: Viewport = get_viewport()
	var sprite_pos = vp.size - _img_size;
	filtered_view2.position = sprite_pos

	var small_size = _filtered_mat1.dimension()
	var small_img_data: PackedByteArray = GBL.mat_to_image_data([ _filtered_mat1 ])
	var small_img: Image = Image.create_from_data(small_size.y, small_size.x, false,
												  Image.FORMAT_L8, small_img_data)
	var texture = ImageTexture.create_from_image(small_img)
	filtered_view1.set_texture(texture)

	small_img_data = GBL.mat_to_image_data([ _filtered_mat2 ])
	small_img = Image.create_from_data(small_size.y, small_size.x, false, Image.FORMAT_L8, small_img_data)
	texture = ImageTexture.create_from_image(small_img)
	filtered_view2.set_texture(texture)

func _ready():
	var _hp_filt = GBL.new_mat(3)
	_hp_filt.from_array([
		[0.1667,  0.6667, 0.1667],
		[0.6667, -3.3333, 0.6667],
		[0.1667,  0.6667, 0.1667],
	])
	var _lp_filt = GBL.new_mat(3)
	_lp_filt.from_array([
		[0, 1, 0],
		[1, 1, 1],
		[0, 1, 0],
	])
	_lp_filt.div(_lp_filt.integrate()) # Normalize
	_edge_filter = _hp_filt.conv(_lp_filt)
	_motion_filt = GBL.new_mat(5, 5)
	_motion_filt.from_array([
		[ 0,           0,      0, 0.0501, 0.0304 ],
		[ 0,           0, 0.0519, 0.1771, 0.0501 ],
		[ 0,      0.0519, 0.1771, 0.0519,      0 ],
		[ 0.0501, 0.1771, 0.0519,      0,      0 ],
		[ 0.0304, 0.0501,      0,      0,      0 ],
	])

func _input(event):
	if event.is_action_pressed("toggle_occlusion_culling"):
		get_viewport().use_occlusion_culling = not get_viewport().use_occlusion_culling
		update_labels()
	if event.is_action_pressed("toggle_mesh_lod"):
		get_viewport().mesh_lod_threshold = 1.0 if is_zero_approx(get_viewport().mesh_lod_threshold) else 0.0
		update_labels()
	if event.is_action_pressed("cycle_draw_mode"):
		get_viewport().debug_draw = wrapi(get_viewport().debug_draw + 1, 0, 5)
		update_labels()

var time_accum: float = 0.0
func _process(delta):
	$Performance.text = """%d FPS (%.2f mspf)

Currently rendering:
%d objects
%dK primitive indices
%d draw calls
""" % [
	Engine.get_frames_per_second(),
	1000.0 / Engine.get_frames_per_second(),
	RenderingServer.get_rendering_info(RenderingServer.RENDERING_INFO_TOTAL_OBJECTS_IN_FRAME),
	RenderingServer.get_rendering_info(RenderingServer.RENDERING_INFO_TOTAL_PRIMITIVES_IN_FRAME) * 0.001,
	RenderingServer.get_rendering_info(RenderingServer.RENDERING_INFO_TOTAL_DRAW_CALLS_IN_FRAME),
]

	if time_accum > DRAW_PERIOD:
		_draw_filtered_viewport()
		time_accum = 0.0

	time_accum += delta

var physics_time_accum: float = 0.0
func _physics_process(delta):

	if physics_time_accum > CAPTURE_PERIOD:
		_capture_and_process()
		physics_time_accum = 0.0

	physics_time_accum += delta

func update_labels():
	$OcclusionCulling.text = "Occlusion culling: %s" % ("Enabled" if get_viewport().use_occlusion_culling else "Disabled")
	$MeshLOD.text = "Mesh LOD: %s" % ("Enabled" if not is_zero_approx(get_viewport().mesh_lod_threshold) else "Disabled")
	$DrawMode.text = "Draw mode: %s" % get_draw_mode_string(get_viewport().debug_draw)

func get_draw_mode_string(draw_mode):
	match draw_mode:
		0:
			return "Normal"
		1:
			return "Unshaded"
		2:
			return "Lighting"
		3:
			return "Overdraw"
		4:
			return "Wireframe"
