#include "register_vk_extensions.h"

#include <drivers/vulkan/vulkan_context.h>

void vk_extensions_request_atomic()
{
	// TODO: Update when float extension is added
	// VulkanContext::register_requested_device_extension(VK_KHR_SHADER_ATOMIC_FLOAT_EXTENSION_NAME, true);
	VulkanContext::register_requested_device_extension(VK_KHR_SHADER_ATOMIC_INT64_EXTENSION_NAME, true);
}
