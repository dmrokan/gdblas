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

var _thread: Thread
var _continue_thread: bool = true
var _mutex: Mutex
var _viewport_image: Image
var _vp_image_mutex: Mutex

func _capture_and_process():
	var img: Image = Image.new()

	_vp_image_mutex.lock()
	if _viewport_image == null:
		_vp_image_mutex.unlock()
		return

	img.copy_from(_viewport_image)
	_vp_image_mutex.unlock()

	var img_size: Vector2i = img.get_size() / 3
	img.resize(img_size.x, img_size.y)
	if img_size != _img_size:
		R.resize(img_size.y, img_size.x)
		G.resize(img_size.y, img_size.x)
		B.resize(img_size.y, img_size.x)
		_img_size = img_size

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

	_mutex.lock()
	_filtered_mat1 = R.conv(_edge_filter, true)
	_filtered_mat2 = R.conv(_motion_filt, true)
	_mutex.unlock()

func _draw_filtered_viewport():
	_mutex.lock()

	if _filtered_mat1 == null or _filtered_mat2 == null:
		_mutex.unlock()
		return

	var small_img_data1: PackedByteArray = GBL.mat_to_image_data([ _filtered_mat1 ])
	var small_size1 = _filtered_mat1.dimension()
	var small_img_data2: PackedByteArray = GBL.mat_to_image_data([ _filtered_mat2 ])
	var small_size2 = _filtered_mat2.dimension()
	_mutex.unlock()

	var vp: Viewport = get_viewport()
	var sprite_pos = vp.size - _img_size;
	filtered_view2.position = sprite_pos

	var small_img: Image = Image.create_from_data(small_size1.y, small_size1.x, false, Image.FORMAT_L8, small_img_data1)
	var texture = ImageTexture.create_from_image(small_img)
	filtered_view1.set_texture(texture)

	small_img = Image.create_from_data(small_size2.y, small_size2.x, false, Image.FORMAT_L8, small_img_data2)
	texture = ImageTexture.create_from_image(small_img)
	filtered_view2.set_texture(texture)

func _image_processing_thread():
	const expected_loop_time: int = round(CAPTURE_PERIOD * 1000.0)

	while _continue_thread:
		var loop_start = Time.get_ticks_msec()

		_capture_and_process()

		var loop_end = Time.get_ticks_msec()
		var dt = loop_end - loop_start
		if dt > 0:
			OS.delay_msec(dt)

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

	_thread = Thread.new()
	_mutex = Mutex.new()
	_vp_image_mutex = Mutex.new()

	_continue_thread = true
	_thread.start(_image_processing_thread)

func _notification(what):
	if what == NOTIFICATION_WM_CLOSE_REQUEST:
		_continue_thread = false
		_thread.wait_to_finish()
		get_tree().quit() # default behavior

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
var total_time: float = 0.0
var proc_counter: int = 0
var avg_fps: float = 0.0
func _process(delta):
	if proc_counter == 0:
		avg_fps = Engine.get_frames_per_second()
		proc_counter += 1
	else:
		avg_fps = avg_fps * proc_counter + Engine.get_frames_per_second()
		proc_counter += 1
		avg_fps /= proc_counter

		if total_time > 2.0:
			total_time = 0.0
			proc_counter = 0

			$Performance.text = """%d FPS (%.2f mspf)

		Currently rendering:
		%d objects
		%dK primitive indices
		%d draw calls
		""" % [
			avg_fps,
			1000.0 / avg_fps,
			RenderingServer.get_rendering_info(RenderingServer.RENDERING_INFO_TOTAL_OBJECTS_IN_FRAME),
			RenderingServer.get_rendering_info(RenderingServer.RENDERING_INFO_TOTAL_PRIMITIVES_IN_FRAME) * 0.001,
			RenderingServer.get_rendering_info(RenderingServer.RENDERING_INFO_TOTAL_DRAW_CALLS_IN_FRAME),
		]

	total_time += delta

	if time_accum > DRAW_PERIOD:
		var vp: Viewport = get_viewport()

		_vp_image_mutex.lock()
		_viewport_image = vp.get_texture().get_image()
		_vp_image_mutex.unlock()

		_draw_filtered_viewport()
		time_accum = 0.0

	time_accum += delta

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
