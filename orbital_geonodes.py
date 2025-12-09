"""
Blender 5.0 script to create hydrogen orbital visualizations using Geometry Nodes.
Organizes orbitals by quantum numbers:
- X axis: angular momentum (l) - 0, 1, 2, 3... (s, p, d, f...)
- Y axis: magnetic quantum number (m) - varies by l
- Z axis: principal quantum number (n) - 1, 2, 3...
"""

import bpy
import math

def clear_scene():
    """Remove all mesh objects from the scene."""
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)

def create_orbital_point(n, l, m, spacing=3.0):
    """
    Create a point in 3D space for an orbital configuration.
    
    Args:
        n: principal quantum number (1, 2, 3...)
        l: angular momentum quantum number (0 to n-1)
        m: magnetic quantum number (-l to +l)
        spacing: distance between orbitals
    """
    x = l * spacing
    y = m * spacing
    z = (n - 1) * - spacing + 9
    return (x, y, z)

def create_base_mesh_with_points(max_n=4):
    """
    Create a mesh with vertices representing all valid orbital configurations.
    
    Args:
        max_n: maximum principal quantum number to generate
    """
    vertices = []
    
    for n in range(1, max_n + 1):
        for l in range(0, n):  # l goes from 0 to n-1
            for m in range(0, l + 1):  # m goes from -l to +l
                pos = create_orbital_point(n, l, m)
                vertices.append(pos)
    
    # Create mesh and object
    mesh = bpy.data.meshes.new("OrbitalPoints")
    mesh.from_pydata(vertices, [], [])
    mesh.update()
    
    obj = bpy.data.objects.new("OrbitalGrid", mesh)
    bpy.context.collection.objects.link(obj)
    
    # Store quantum number data as custom attributes
    n_attr = mesh.attributes.new(name="n", type='INT', domain='POINT')
    l_attr = mesh.attributes.new(name="l", type='INT', domain='POINT')
    m_attr = mesh.attributes.new(name="m", type='INT', domain='POINT')
    
    idx = 0
    for n in range(1, max_n + 1):
        for l in range(0, n):
            for m in range(0, l + 1):
                n_attr.data[idx].value = n
                l_attr.data[idx].value = l
                m_attr.data[idx].value = m
                idx += 1
    
    return obj

def create_instance_sphere():
    """Create a sphere to be instanced at each orbital point."""
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.125, location=(0, 0, 0))
    sphere = bpy.context.active_object
    sphere.name = "OrbitalSphere"
    
    # Add material with color variation
    mat = bpy.data.materials.new(name="OrbitalMaterial")
    mat.use_nodes = True
    sphere.data.materials.append(mat)
    
    return sphere

def setup_geometry_nodes(base_obj, instance_obj):
    """
    Set up geometry nodes for instancing on points.
    
    Args:
        base_obj: the object with points (orbital positions)
        instance_obj: the object to instance at each point
    """
    # Add geometry nodes modifier
    modifier = base_obj.modifiers.new(name="OrbitalInstances", type='NODES')
    
    # Create new node tree
    node_tree = bpy.data.node_groups.new(name="OrbitalGeometry", type='GeometryNodeTree')
    modifier.node_group = node_tree
    
    # Create nodes
    nodes = node_tree.nodes
    links = node_tree.links
    
    # Input and Output nodes
    input_node = nodes.new(type='NodeGroupInput')
    output_node = nodes.new(type='NodeGroupOutput')
    input_node.location = (-600, 0)
    output_node.location = (600, 0)
    
    # Create input and output sockets
    node_tree.interface.new_socket(name="Geometry", in_out='INPUT', socket_type='NodeSocketGeometry')
    node_tree.interface.new_socket(name="Geometry", in_out='OUTPUT', socket_type='NodeSocketGeometry')
    
    # Mesh to Points node
    mesh_to_points = nodes.new(type='GeometryNodeMeshToPoints')
    mesh_to_points.location = (-400, 0)
    
    # Instance on Points node
    instance_on_points = nodes.new(type='GeometryNodeInstanceOnPoints')
    instance_on_points.location = (-200, 0)
    
    # Object Info node for the instance object
    object_info = nodes.new(type='GeometryNodeObjectInfo')
    object_info.location = (-400, -200)
    object_info.inputs['Object'].default_value = instance_obj
    object_info.transform_space = 'RELATIVE'
    
    # Named Attribute nodes for quantum numbers
    n_attribute = nodes.new(type='GeometryNodeInputNamedAttribute')
    n_attribute.location = (-400, -350)
    n_attribute.inputs['Name'].default_value = "n"
    n_attribute.data_type = 'INT'
    
    l_attribute = nodes.new(type='GeometryNodeInputNamedAttribute')
    l_attribute.location = (-400, -500)
    l_attribute.inputs['Name'].default_value = "l"
    l_attribute.data_type = 'INT'
    
    # Map Range for scaling based on n
    scale_map = nodes.new(type='ShaderNodeMapRange')
    scale_map.location = (-200, -350)
    scale_map.inputs['From Min'].default_value = 1.0
    scale_map.inputs['From Max'].default_value = 4.0
    scale_map.inputs['To Min'].default_value = 0.5
    scale_map.inputs['To Max'].default_value = 1.5
    
    # Mix RGB node for color variation based on l
    # Create color for low l (blue)
    color_low = nodes.new(type='FunctionNodeCombineColor')
    color_low.location = (-200, -600)
    color_low.inputs['Red'].default_value = 0.1
    color_low.inputs['Green'].default_value = 0.3
    color_low.inputs['Blue'].default_value = 1.0
    
    # Create color for high l (red)
    color_high = nodes.new(type='FunctionNodeCombineColor')
    color_high.location = (-200, -750)
    color_high.inputs['Red'].default_value = 1.0
    color_high.inputs['Green'].default_value = 0.3
    color_high.inputs['Blue'].default_value = 0.1
    
    # Mix colors based on l
    color_mix = nodes.new(type='ShaderNodeMix')
    color_mix.location = (0, -600)
    color_mix.data_type = 'RGBA'
    color_mix.blend_type = 'MIX'
    
    # Map Range to normalize l for color mixing
    normalize_l = nodes.new(type='ShaderNodeMapRange')
    normalize_l.location = (-200, -500)
    normalize_l.inputs['From Min'].default_value = 0.0
    normalize_l.inputs['From Max'].default_value = 3.0
    normalize_l.inputs['To Min'].default_value = 0.0
    normalize_l.inputs['To Max'].default_value = 1.0
    
    # Set Material node
    set_material = nodes.new(type='GeometryNodeSetMaterial')
    set_material.location = (400, 0)
    if instance_obj.data.materials:
        set_material.inputs['Material'].default_value = instance_obj.data.materials[0]
    
    # Store Named Attribute node for color
    store_color = nodes.new(type='GeometryNodeStoreNamedAttribute')
    store_color.location = (200, 0)
    store_color.inputs['Name'].default_value = "orbital_color"
    store_color.data_type = 'FLOAT_COLOR'
    
    # Link nodes
    links.new(input_node.outputs['Geometry'], mesh_to_points.inputs['Mesh'])
    links.new(mesh_to_points.outputs['Points'], instance_on_points.inputs['Points'])
    links.new(object_info.outputs['Geometry'], instance_on_points.inputs['Instance'])
    
    # Link quantum number attributes
    links.new(n_attribute.outputs['Attribute'], scale_map.inputs['Value'])
    links.new(scale_map.outputs['Result'], instance_on_points.inputs['Scale'])
    
    # Link color mixing
    links.new(l_attribute.outputs['Attribute'], normalize_l.inputs['Value'])
    links.new(normalize_l.outputs['Result'], color_mix.inputs['Factor'])
    links.new(color_low.outputs['Color'], color_mix.inputs[6])  # A input
    links.new(color_high.outputs['Color'], color_mix.inputs[7])  # B input
    
    # Link color to store attribute
    links.new(instance_on_points.outputs['Instances'], store_color.inputs['Geometry'])
    links.new(color_mix.outputs[2], store_color.inputs['Value'])  # Result output
    
    # Link to material and output
    links.new(store_color.outputs['Geometry'], set_material.inputs['Geometry'])
    links.new(set_material.outputs['Geometry'], output_node.inputs['Geometry'])

def add_text_labels():
    """Add text labels for axis identification."""
    label_data = [
        ("X: Angular Momentum (l)", (6, -2, 0)),
        ("Y: Magnetic Number (m)", (0, 6, 0)),
        ("Z: Principal Number (n)", (0, 0, 10))
    ]
    
    for text, location in label_data:
        bpy.ops.object.text_add(location=location)
        text_obj = bpy.context.active_object
        text_obj.data.body = text
        text_obj.scale = (0.5, 0.5, 0.5)

def setup_camera_and_lighting():
    """Set up camera and lighting for the scene."""
    # Add camera
    bpy.ops.object.camera_add(location=(15, -15, 12))
    camera = bpy.context.active_object
    camera.rotation_euler = (math.radians(60), 0, math.radians(45))
    bpy.context.scene.camera = camera
    
    # Add sun light
    bpy.ops.object.light_add(type='SUN', location=(10, 10, 10))
    sun = bpy.context.active_object
    sun.data.energy = 2.0
    
    # Add area light for fill
    bpy.ops.object.light_add(type='AREA', location=(-10, -10, 5))
    area = bpy.context.active_object
    area.data.energy = 500

def main():
    """Main function to create the orbital visualization."""
    # Clear existing scene
    clear_scene()
    
    # Create the base mesh with orbital points
    max_n = 4  # Generate orbitals up to n=4
    base_obj = create_base_mesh_with_points(max_n=max_n)
    
    # Create instance object (sphere)
    instance_obj = create_instance_sphere()
    
    # Hide the instance object from viewport and render
    instance_obj.hide_viewport = True
    instance_obj.hide_render = True
    
    # Set up geometry nodes
    setup_geometry_nodes(base_obj, instance_obj)
    
    # Add text labels
    add_text_labels()
    
    # Set up camera and lighting
    setup_camera_and_lighting()
    
    # Select the main orbital grid object
    bpy.context.view_layer.objects.active = base_obj
    base_obj.select_set(True)
    
    print(f"Orbital visualization created with n=1 to n={max_n}")
    print("X axis: Angular momentum (l) - 0=s, 1=p, 2=d, 3=f")
    print("Y axis: Magnetic quantum number (m)")
    print("Z axis: Principal quantum number (n)")

if __name__ == "__main__":
    main()
