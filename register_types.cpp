#include "register_types.h"

#include "triangle_ray_select.h"
#include "triangle_transform.h"

#include <core/object/class_db.h>

void initialize_gd_module_triangle_ray_select_module(ModuleInitializationLevel p_level)
{
	if(p_level == MODULE_INITIALIZATION_LEVEL_CORE)
	{
		TriangleRaySelect::vk_extensions_request_atomic();
	}
	else if(p_level == MODULE_INITIALIZATION_LEVEL_SCENE)
	{
		ClassDB::register_class<TriangleTransform>();
		ClassDB::register_class<MeshTrianglePoint>();
		ClassDB::register_class<TriangleRaySelect>();
	}
}

void uninitialize_gd_module_triangle_ray_select_module(ModuleInitializationLevel p_level) {}
