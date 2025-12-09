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
    Create point clouds for all orbital configurations.
    
    Args:
        max_n: maximum principal quantum number to generate
    """
    all_vertices = []
    n_values = []
    l_values = []
    m_values = []
    phase_values = []
    
    print(f"Generating orbital clouds for n=1 to n={max_n}...")
    
    for n in range(1, max_n + 1):
        for l in range(0, n):  # l goes from 0 to n-1
            for m in range(0, l + 1):  # m goes from -l to +l
                pos = create_orbital_point(n, l, m)
                
                # Create point cloud for this orbital
                print(f"Creating orbital n={n}, l={l}, m={m}...")
                orbital_points, orbital_phases = create_orbital_cloud_mesh(n, l, m, pos, num_points=300)
                
                if len(orbital_points) == 0:
                    print(f"Warning: No points generated for n={n}, l={l}, m={m}")
                    continue
                    
                all_vertices.extend(orbital_points)
                phase_values.extend(orbital_phases)
                
                # Store quantum numbers for each point
                for _ in orbital_points:
                    n_values.append(n)
                    l_values.append(l)
                    m_values.append(m)
    
    if len(all_vertices) == 0:
        print("ERROR: No vertices generated!")
        return None
    
    print(f"Total points generated: {len(all_vertices)}")
    
    # Create mesh and object
    mesh = bpy.data.meshes.new("OrbitalClouds")
    mesh.from_pydata(all_vertices, [], [])
    mesh.update()
    
    obj = bpy.data.objects.new("OrbitalGrid", mesh)
    bpy.context.collection.objects.link(obj)
    
    # Store quantum number data as custom attributes
    n_attr = mesh.attributes.new(name="n", type='INT', domain='POINT')
    l_attr = mesh.attributes.new(name="l", type='INT', domain='POINT')
    m_attr = mesh.attributes.new(name="m", type='INT', domain='POINT')
    sign_attr = mesh.attributes.new(name="wf_sign", type='FLOAT', domain='POINT')
    
    for idx, (n, l, m, sign) in enumerate(zip(n_values, l_values, m_values, phase_values)):
        n_attr.data[idx].value = n
        l_attr.data[idx].value = l
        m_attr.data[idx].value = m
        sign_attr.data[idx].value = sign
    
    return obj

def create_orbital_cloud_mesh(n, l, m, position, num_points=500):
    """Create a point cloud mesh for a specific orbital."""
    import numpy as np
    from scipy.special import sph_harm
    from scipy.special import assoc_laguerre
    
    points = []
    phases = []
    radius_scale = 0.5
    
    # Generate points using rejection sampling
    attempts = 0
    max_attempts = num_points * 50
    
    while len(points) < num_points and attempts < max_attempts:
        attempts += 1
        
        # Sample radius with appropriate distribution
        r = np.random.exponential(n * radius_scale)
        if r > n * 3 * radius_scale:  # cutoff
            continue
            
        theta = np.arccos(2 * np.random.random() - 1)
        phi = 2 * np.pi * np.random.random()
        
        # Radial wave function
        rho = 2 * r / (n * radius_scale)
        if rho < 0.01:  # avoid singularity
            continue
            
        try:
            L = assoc_laguerre(rho, n - l - 1, 2 * l + 1)
            R = np.exp(-rho / 2) * (rho ** l) * L
            
            # Angular wave function (complex valued)
            Y = sph_harm(abs(m), l, phi, theta)
            
            # Probability density
            prob = (R ** 2) * (abs(Y) ** 2) * (r ** 2)
            
            # Rejection sampling
            if prob > 0 and np.random.random() < min(prob * 100, 1.0):
                x = r * np.sin(theta) * np.cos(phi) + position[0]
                y = r * np.sin(theta) * np.sin(phi) + position[1]
                z = r * np.cos(theta) + position[2]
                points.append((x, y, z))
                
                # Calculate sign of the wave function (real part)
                psi = R * Y.real  # Use real part for sign
                sign = 1.0 if psi >= 0 else -1.0
                phases.append(float(sign))
        except:
            continue
    
    # If we didn't get enough points, fill with simple sphere
    while len(points) < num_points:
        r = np.random.exponential(n * radius_scale)
        theta = np.arccos(2 * np.random.random() - 1)
        phi = 2 * np.pi * np.random.random()
        x = r * np.sin(theta) * np.cos(phi) + position[0]
        y = r * np.sin(theta) * np.sin(phi) + position[1]
        z = r * np.cos(theta) + position[2]
        points.append((x, y, z))
        phases.append(0.0)  # default phase for fill points
    
    return points[:num_points], phases[:num_points]

def create_instance_sphere():
    """Create a tiny sphere to be instanced at each point in the orbital cloud."""
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.03, location=(0, 0, 0))
    sphere = bpy.context.active_object
    sphere.name = "OrbitalParticle"
    
    # Add material with sign-based color emission (no diffuse)
    mat = bpy.data.materials.new(name="OrbitalMaterial")
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.clear()
    
    # Output node
    output = nodes.new(type='ShaderNodeOutputMaterial')
    output.location = (400, 0)
    
    # Emission shader (only emission, no diffuse)
    emission = nodes.new(type='ShaderNodeEmission')
    emission.location = (200, 0)
    emission.inputs['Strength'].default_value = 3.0
    
    # Named Attribute node for wave function sign
    sign_attr = nodes.new(type='ShaderNodeAttribute')
    sign_attr.location = (-400, 0)
    sign_attr.attribute_name = "wf_sign"
    
    # Map Range to convert sign from [-1, 1] to [0, 1] for color mixing
    # -1 (negative) -> 0 (cyan), +1 (positive) -> 1 (mint)
    map_range = nodes.new(type='ShaderNodeMapRange')
    map_range.location = (-200, 0)
    map_range.inputs['From Min'].default_value = -1.0
    map_range.inputs['From Max'].default_value = 1.0
    map_range.inputs['To Min'].default_value = 0.0
    map_range.inputs['To Max'].default_value = 1.0
    
    # Cyan color (RGB: 0, 255/255, 255/255) for negative wave function
    cyan_color = nodes.new(type='ShaderNodeRGB')
    cyan_color.location = (-200, -200)
    cyan_color.outputs['Color'].default_value = (0.0, 1.0, 1.0, 1.0)
    
    # Mint color (RGB: 152/255, 255/255, 152/255) for positive wave function
    mint_color = nodes.new(type='ShaderNodeRGB')
    mint_color.location = (-200, -350)
    mint_color.outputs['Color'].default_value = (0.596, 1.0, 0.596, 1.0)
    
    # Mix colors based on sign
    mix_color = nodes.new(type='ShaderNodeMix')
    mix_color.location = (0, 0)
    mix_color.data_type = 'RGBA'
    
    # Link nodes
    links.new(sign_attr.outputs['Fac'], map_range.inputs['Value'])
    links.new(map_range.outputs['Result'], mix_color.inputs['Factor'])
    links.new(cyan_color.outputs['Color'], mix_color.inputs[6])  # A (negative/cyan)
    links.new(mint_color.outputs['Color'], mix_color.inputs[7])  # B (positive/mint)
    links.new(mix_color.outputs[2], emission.inputs['Color'])
    links.new(emission.outputs['Emission'], output.inputs['Surface'])
    
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
    mesh_to_points.inputs['Radius'].default_value = 0.02
    
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
    
    # Named Attribute node for wave function sign
    sign_attribute = nodes.new(type='GeometryNodeInputNamedAttribute')
    sign_attribute.location = (-400, -650)
    sign_attribute.inputs['Name'].default_value = "wf_sign"
    sign_attribute.data_type = 'FLOAT'
    
    # Store sign attribute to pass through to instances
    store_sign = nodes.new(type='GeometryNodeStoreNamedAttribute')
    store_sign.location = (0, 0)
    store_sign.inputs['Name'].default_value = "wf_sign"
    store_sign.data_type = 'FLOAT'
    
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
    
    # Link sign attribute to store
    links.new(instance_on_points.outputs['Instances'], store_sign.inputs['Geometry'])
    links.new(sign_attribute.outputs['Attribute'], store_sign.inputs['Value'])
    
    # Link color to store attribute
    links.new(store_sign.outputs['Geometry'], store_color.inputs['Geometry'])
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
    
    if base_obj is None:
        print("Failed to create orbital mesh!")
        return
    
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
