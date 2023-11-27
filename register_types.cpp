#include "register_types.h"

#include "triangle_ray_select.h"

#include <core/object/class_db.h>

void initialize_triangle_ray_select_module(ModuleInitializationLevel p_level)
{
	if(p_level == MODULE_INITIALIZATION_LEVEL_CORE)
	{
		TriangleRaySelect::vk_extensions_request_atomic();
	}
	else if(p_level == MODULE_INITIALIZATION_LEVEL_SCENE)
	{
		ClassDB::register_class<MeshSurfaceIndex>();
		ClassDB::register_class<TriangleRaySelect>();
	}
}

void uninitialize_triangle_ray_select_module(ModuleInitializationLevel p_level) {}
