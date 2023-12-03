#[compute]
#version 450

#VERSION_DEFINES
//#extension GL_EXT_shader_atomic_float2 : enable
#extension GL_EXT_shader_atomic_int64 : enable

layout(local_size_x = 64, local_size_y = 1, local_size_z = 1) in;

layout(set=0, binding=0, std430) buffer restrict readonly IndexData {
	uint data[];
} index_data;

layout(set=0, binding=1, std430) buffer restrict readonly VertexData {
	uint data[];
} vertex_data;

layout(set=1, binding=0, std430) buffer restrict SelectedVertex {
	// Triangle closest to the origin (triangle_index is the first index_buffer id of the triangle)
	uint triangle_index;
	uint origin_dist;

	uint vertex_ids[3];
	uint pad_0;

	float point_on_triangle[3];
	uint pad_1;
} selected_triangle_data;

layout(push_constant, std430) uniform Params {
	uint index_count;
	uint index_stride;
	uint vertex_stride;
	uint pad0;

	vec3 ray_origin;
	uint pad1;

	vec3 ray_normal;
	uint pad2;
} params;


// Triangle-Ray Intersection
// (code from https://stackoverflow.com/questions/59257678/intersect-a-ray-with-a-triangle-in-glsl-c)
bool PointInOrOn( vec3 P1, vec3 P2, vec3 A, vec3 B ) {
	vec3 CP1 = cross(B - A, P1 - A);
	vec3 CP2 = cross(B - A, P2 - A);
	return dot(CP1, CP2) >= 0.0;
}

bool PointInTriangle( vec3 px, vec3 p0, vec3 p1, vec3 p2 ) {
	return PointInOrOn(px, p0, p1, p2) &&
	    PointInOrOn(px, p1, p2, p0) &&
	    PointInOrOn(px, p2, p0, p1);
}

vec3 PlaneNormal(vec3 p0, vec3 p1, vec3 p2) {
	return cross(p1-p0, p2-p0);
}

vec3 IntersectPlane(vec3 ray_origin, vec3 ray_normal, vec3 p0, vec3 plane_normal, float dotp) {
	vec3 X = ray_origin + ray_normal * dot(p0 - ray_origin, plane_normal) / dotp;

	return X;
}

bool IntersectTriangle(vec3 point_on_plane, vec3 p0, vec3 p1, vec3 p2) {
	return PointInTriangle(point_on_plane, p0, p1, p2);
}

void main() {
	const uint index = gl_GlobalInvocationID.x * params.index_stride;
	if (index >= params.index_count) {
		return;
	}

	// Get vertex data
	uint vert_index = index_data.data[index + 0] * params.vertex_stride;
	const vec3 vertex_a = uintBitsToFloat(uvec3(vertex_data.data[vert_index + 0],
	                                            vertex_data.data[vert_index + 1],
	                                            vertex_data.data[vert_index + 2]));

	vert_index = index_data.data[index + 1] * params.vertex_stride;
	const vec3 vertex_b = uintBitsToFloat(uvec3(vertex_data.data[vert_index + 0],
	                                            vertex_data.data[vert_index + 1],
	                                            vertex_data.data[vert_index + 2]));

	vert_index = index_data.data[index + 2] * params.vertex_stride;
	const vec3 vertex_c = uintBitsToFloat(uvec3(vertex_data.data[vert_index + 0],
	                                            vertex_data.data[vert_index + 1],
	                                            vertex_data.data[vert_index + 2]));

	// Check if ray intersects triangle
	const vec3 face_normal = PlaneNormal(vertex_a, vertex_b, vertex_c);
	const float dotp = dot(params.ray_normal, face_normal);
	if(dotp <= 0.0) {
		return;
	}

	const vec3 point_on_plane = IntersectPlane(params.ray_origin, params.ray_normal, vertex_a, face_normal, dotp);
	if(IntersectTriangle(point_on_plane, vertex_a, vertex_b, vertex_c))	{
		// TODO: Replace uint comparison with float comparison once float atomics are added
		const vec3 origin_dist_vec = point_on_plane - params.ray_origin;
		const uint origin_dist     = uint(dot(origin_dist_vec, origin_dist_vec) * 1000.0);

		// Check for triangle that's closest to ray origin
		atomicMin(selected_triangle_data.origin_dist, origin_dist);
		barrier();
		if(origin_dist == selected_triangle_data.origin_dist) {
			selected_triangle_data.triangle_index = index;

			selected_triangle_data.vertex_ids[0] = index_data.data[index + 0];
			selected_triangle_data.vertex_ids[1] = index_data.data[index + 1];
			selected_triangle_data.vertex_ids[2] = index_data.data[index + 2];

			selected_triangle_data.point_on_triangle[0] = point_on_plane[0];
			selected_triangle_data.point_on_triangle[1] = point_on_plane[1];
			selected_triangle_data.point_on_triangle[2] = point_on_plane[2];
		}
	}
}
