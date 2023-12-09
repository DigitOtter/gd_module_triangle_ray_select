# gd\_module\_triangle\_ray\_select

Triangle selection from ray. `TriangleRaySelect` takes an array of `MeshInstance3D`'s and a ray origin and normal to select the triangle that's closest to the ray's origin. The class checks the mesh's vertice positions after skeleton and blend shape transforms have been applied, so we can use up-to-date information for triangle selection.

Note:
- The class needs access to Godot's `MeshStorage` private members, so we have to apply a patch to the source code (see Installation below)
- Only implemented for Vulkan, not for OpenGL.


## Installation

- Clone Godot Engine from [https://github.com/godotengine/godot](https://github.com/godotengine/godot)
- Clone this repository into the `modules` subdirectory
- Run `git apply modules/gd_module_triangle_ray_select/mesh_storage.patch` from godot's root directory
- Compile godot [as usual](https://docs.godotengine.org/en/stable/contributing/development/compiling/index.html)


## Usage

An example project is available at [https://github.com/DigitOtter/gd_triangle_ray_select_example](https://github.com/DigitOtter/gd_triangle_ray_select_example).

https://github.com/DigitOtter/gd_triangle_ray_select_example/assets/102542997/ff05f6cc-ca81-4dab-b132-58879b1b54c2

TriangleRaySelect example to get the triangle at mouse position on click:
``` gdscript
extends Node

@onready
var triangle_ray_select: TriangleRaySelect = TriangleRaySelect.new()

var selectable_meshes: Array[MeshInstance3D]

func _ready():
	self.selectable_meshes = [ <ARRAY OF INTERSECTION MESHES> ]

func _input(event):
	if event.is_pressed() and event is InputEventMouseButton and event.button_index == MOUSE_BUTTON_LEFT:
		var camera: Camera3D = get_viewport().get_camera_3d()
		var pixel: Vector2i  = event.position
		
		# Select the triangle that's closest to the camera's origin
		var mesh_triangle_point: MeshTrianglePoint = self.triangle_ray_select.select_triangle_from_meshes_cam(self.selectable_meshes, camera, pixel)
		print(mesh_triangle_point)
		print("Point on triangle:        ", mesh_triangle_point.point_on_triangle)
		print("Distance to origin:       ", mesh_triangle_point.ray_origin_dist)
		print("Triangle vertice indices: ", mesh_triangle_point.vertex_ids)
```

TriangleTransform example to locate the global pose of a selected triangle:
``` gdscript
extends Node

var triangle_ray_select: TriangleRaySelect = null
var mesh_triangle_point: MeshTrianglePoint = null
var triangle_transform: TriangleTransform = null

var selectable_meshes: Array[MeshInstance3D]

func set_mesh_triangle_point(_triangle_ray_select: TriangleRaySelect, _mesh_triangle_point: MeshTrianglePoint):
	self.triangle_ray_select = _triangle_ray_select
	self.mesh_triangle_point = _mesh_triangle_point
	
	# Create transform for point_on_triangle from triangle vertices
	self.triangle_transform  = self.triangle_ray_select.get_triangle_transform_msi(self.mesh_triangle_point, Transform3D(Basis(), self.mesh_triangle_point.point_on_triangle))

func _process(_delta):
	if self.mesh_triangle_point and self.mesh_triangle_point.mesh_instance:
		# Get current triangle vertice positions
		var updated_vertice_positions: PackedVector3Array = self.triangle_ray_select.get_triangle_vertices(self.mesh_triangle_point)
		
		# Get point pose on triangle
		var updated_point_pose: Transform3D = self.triangle_transform.adjust_transform(updated_vertice_positions, 0.0)
		print(updated_point_pose)
```
