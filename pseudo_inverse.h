#pragma once

#include "core/math/vector3.h"
#include "core/variant/variant.h"

Vector3 pseudo_inverse_mult(const PackedVector3Array &array, const Vector3 &val);

Vector3 pseudo_inverse_mult(const Vector3 &vec0, const Vector3 &vec1, const Vector3 &vec2, const Vector3 &val);
