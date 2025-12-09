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
    # Use direct removal instead of operators to handle unselectable objects
    objects_to_remove = [obj for obj in bpy.data.objects]
    for obj in objects_to_remove:
        bpy.data.objects.remove(obj, do_unlink=True)

def create_orbital_point(n, l, m, spacing=12.0):
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
    z = (4 - n) * spacing
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
    
    for n in range(1, max_n + 1):
        for l in range(0, n):  # l goes from 0 to n-1
            for m in range(0, l + 1):  # m goes from -l to +l
                pos = create_orbital_point(n, l, m)
                
                # Create point cloud for this orbital
                orbital_points, orbital_phases = create_orbital_cloud_mesh(n, l, m, pos, num_points=300)
                all_vertices.extend(orbital_points)
                phase_values.extend(orbital_phases)
                
                # Store quantum numbers for each point
                for _ in orbital_points:
                    n_values.append(n)
                    l_values.append(l)
                    m_values.append(m)
    
    # Check if we have any vertices before proceeding
    if len(all_vertices) == 0:
        print("ERROR: No vertices generated!")
        print("Check if scipy and numpy are installed in Blender's Python")
        return None
    
    print(f"Total points generated: {len(all_vertices)}")
    
    # Create mesh and object
    mesh = bpy.data.meshes.new("OrbitalClouds")
    mesh.from_pydata(all_vertices, [], [])
    mesh.update()
    
    # Verify mesh has vertices
    if len(mesh.vertices) == 0:
        print(f"ERROR: Mesh has no vertices even though all_vertices had {len(all_vertices)} points!")
        print("This shouldn't happen - mesh.from_pydata failed")
        return None
    
    print(f"Mesh created successfully with {len(mesh.vertices)} vertices")
    
    # Create and link object FIRST - attributes need the object in the scene to initialize
    obj = bpy.data.objects.new("OrbitalGrid", mesh)
    bpy.context.collection.objects.link(obj)
    print(f"✓ Object created and linked: {obj.name}")
    
    # NOW create attributes (they will initialize properly because object is in scene)
    print(f"Creating attributes for {len(mesh.vertices)} vertices")
    print(f"n_values: {len(n_values)}, l_values: {len(l_values)}, m_values: {len(m_values)}, phase_values: {len(phase_values)}")
    
    n_attr = mesh.attributes.new(name="n", type='FLOAT', domain='POINT')
    l_attr = mesh.attributes.new(name="l", type='FLOAT', domain='POINT')
    m_attr = mesh.attributes.new(name="m", type='FLOAT', domain='POINT')
    wf_sign_attr = mesh.attributes.new(name="wf_sign", type='FLOAT', domain='POINT')
    
    # Force update after creating attributes
    mesh.update()
    obj.update_from_editmode()
    bpy.context.view_layer.update()
    
    print(f"✓ Attributes created and scene updated")
    print(f"  Attribute sizes: n={len(n_attr.data)}, l={len(l_attr.data)}, m={len(m_attr.data)}, wf_sign={len(wf_sign_attr.data)}")
    print(f"  Value list sizes: n={len(n_values)}, l={len(l_values)}, m={len(m_values)}, phase={len(phase_values)}")
    
    # If attributes still have size 0, try reacquiring them
    if len(n_attr.data) == 0:
        print("  Attributes have size 0, reacquiring from mesh...")
        n_attr = mesh.attributes.get("n")
        l_attr = mesh.attributes.get("l")
        m_attr = mesh.attributes.get("m")
        wf_sign_attr = mesh.attributes.get("wf_sign")
        print(f"  After reacquiring: n={len(n_attr.data)}, l={len(l_attr.data)}, m={len(m_attr.data)}, wf_sign={len(wf_sign_attr.data)}")
    
    # Use foreach_set for efficient bulk assignment
    try:
        n_attr.data.foreach_set("value", n_values)
        l_attr.data.foreach_set("value", l_values)
        m_attr.data.foreach_set("value", m_values)
        wf_sign_attr.data.foreach_set("value", phase_values)
        print(f"✓ Attributes assigned using foreach_set")
    except Exception as e:
        print(f"foreach_set failed, using loop: {e}")
        for idx, (n, l, m, sign) in enumerate(zip(n_values, l_values, m_values, phase_values)):
            n_attr.data[idx].value = n
            l_attr.data[idx].value = l
            m_attr.data[idx].value = m
            wf_sign_attr.data[idx].value = sign
        print(f"✓ Attributes assigned using loop")
    
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
                
                # Calculate phase from wave function
                # Phase is the argument of the complex wave function value
                psi = R * Y
                phase = np.angle(psi) / np.pi  # normalize to [-1, 1]
                phases.append(float(phase))
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
    """Create a tiny icosphere to be instanced at each point in the orbital cloud."""
    # Calculate radius for 1.618x volume: r_new = r_old * (1.618)^(1/3)
    radius = 0.03 * (1.618 ** (1/3))  # ~0.0352
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=1, radius=radius, location=(0, 0, 0))
    sphere = bpy.context.active_object
    sphere.name = "OrbitalParticle"
    
    # Add material with phase-based color emission
    # Remove existing material if it exists to ensure we get the updated version
    if "OrbitalMaterial" in bpy.data.materials:
        bpy.data.materials.remove(bpy.data.materials["OrbitalMaterial"])
    
    mat = bpy.data.materials.new(name="OrbitalMaterial")
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.clear()
    
    # Output node
    output = nodes.new(type='ShaderNodeOutputMaterial')
    output.location = (400, 0)
    
    # Emission shader
    emission = nodes.new(type='ShaderNodeEmission')
    emission.location = (200, 0)
    emission.inputs['Strength'].default_value = 3.0
    
    # Named Attribute node for phase
    phase_attr = nodes.new(type='ShaderNodeAttribute')
    phase_attr.location = (-400, 0)
    phase_attr.attribute_name = "phase"
    
    # Greater Than comparison - true (1) if phase > 0.5, false (0) otherwise
    # Since phase is now in unorm range [0, 1], compare to 0.5 to distinguish positive/negative
    greater_than = nodes.new(type='ShaderNodeMath')
    greater_than.location = (-200, 0)
    greater_than.operation = 'GREATER_THAN'
    greater_than.inputs[1].default_value = 0.5  # threshold
    
    # Mint color (RGB: 152/255, 255/255, 152/255)
    mint_color = nodes.new(type='ShaderNodeRGB')
    mint_color.location = (-200, -200)
    mint_color.outputs['Color'].default_value = (0.596, 1.0, 0.596, 1.0)
    
    # Magenta color (RGB: 255/255, 0/255, 255/255)
    magenta_color = nodes.new(type='ShaderNodeRGB')
    magenta_color.location = (-200, -350)
    magenta_color.outputs['Color'].default_value = (1.0, 0.0, 1.0, 1.0)
    
    # Mix colors based on phase
    mix_color = nodes.new(type='ShaderNodeMix')
    mix_color.location = (0, 0)
    mix_color.data_type = 'RGBA'
    
    # Link nodes
    links.new(phase_attr.outputs['Fac'], greater_than.inputs[0])
    links.new(greater_than.outputs['Value'], mix_color.inputs['Factor'])
    links.new(mint_color.outputs['Color'], mix_color.inputs[6])  # A
    links.new(magenta_color.outputs['Color'], mix_color.inputs[7])  # B
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
    print(f"Setting up geometry nodes...")
    print(f"  Base object: {base_obj.name}")
    print(f"  Instance object: {instance_obj.name}")
    
    # Clear all materials from the file
    for mat in list(bpy.data.materials):
        bpy.data.materials.remove(mat)
    print("  ✓ Cleared all materials")
    
    # Create new material with magenta/mint emission
    mat = bpy.data.materials.new(name="OrbitalMaterial")
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.clear()
    
    # Output node
    output = nodes.new(type='ShaderNodeOutputMaterial')
    output.location = (400, 0)
    
    # Emission shader
    emission = nodes.new(type='ShaderNodeEmission')
    emission.location = (200, 0)
    emission.inputs['Strength'].default_value = 3.0
    
    # Named Attribute node for wave function sign
    sign_attr = nodes.new(type='ShaderNodeAttribute')
    sign_attr.location = (-400, 0)
    sign_attr.attribute_name = "wf_sign"
    
    # Map Range
    map_range = nodes.new(type='ShaderNodeMapRange')
    map_range.location = (-200, 0)
    map_range.inputs['From Min'].default_value = -1.0
    map_range.inputs['From Max'].default_value = 1.0
    
    # Magenta color
    magenta = nodes.new(type='ShaderNodeRGB')
    magenta.location = (-200, -200)
    magenta.outputs['Color'].default_value = (1.0, 0.0, 1.0, 1.0)
    
    # Mint color
    mint = nodes.new(type='ShaderNodeRGB')
    mint.location = (-200, -350)
    mint.outputs['Color'].default_value = (0.596, 1.0, 0.596, 1.0)
    
    # Mix colors
    mix_color = nodes.new(type='ShaderNodeMix')
    mix_color.location = (0, 0)
    mix_color.data_type = 'RGBA'
    
    # Link material nodes
    links.new(sign_attr.outputs['Fac'], map_range.inputs['Value'])
    links.new(map_range.outputs['Result'], mix_color.inputs['Factor'])
    links.new(magenta.outputs['Color'], mix_color.inputs[6])
    links.new(mint.outputs['Color'], mix_color.inputs[7])
    links.new(mix_color.outputs[2], emission.inputs['Color'])
    links.new(emission.outputs['Emission'], output.inputs['Surface'])
    
    # Apply material to instance object
    if len(instance_obj.data.materials) == 0:
        instance_obj.data.materials.append(mat)
    else:
        instance_obj.data.materials[0] = mat
    print("  ✓ Created and applied OrbitalMaterial")
    
    # Assign material to base object so geometry nodes can use it
    base_obj.data.materials.clear()
    base_obj.data.materials.append(mat)
    print(f"  ✓ Assigned {mat.name} to {base_obj.name}")
    
    # Add geometry nodes modifier
    modifier = base_obj.modifiers.new(name="OrbitalInstances", type='NODES')
    print(f"  Created modifier: {modifier.name}")
    
    # Create new node tree with unique name
    import time
    node_tree_name = f"OrbitalGeometry_{int(time.time())}"
    node_tree = bpy.data.node_groups.new(name=node_tree_name, type='GeometryNodeTree')
    print(f"  Created node group: {node_tree.name}")
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
    
    # Named Attribute node for phase
    phase_attribute = nodes.new(type='GeometryNodeInputNamedAttribute')
    phase_attribute.location = (-400, -650)
    phase_attribute.inputs['Name'].default_value = "phase"
    phase_attribute.data_type = 'FLOAT'
    
    # Store phase attribute to pass through to instances
    store_phase = nodes.new(type='GeometryNodeStoreNamedAttribute')
    store_phase.location = (0, 0)
    store_phase.inputs['Name'].default_value = "phase"
    store_phase.data_type = 'FLOAT'
    
    # Set Material node - use the material we just created
    set_material = nodes.new(type='GeometryNodeSetMaterial')
    set_material.location = (400, 0)
    set_material.inputs['Material'].default_value = mat
    print(f"  ✓ Set Material node configured with {mat.name}")
    
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
    
    # Link phase attribute to store
    links.new(instance_on_points.outputs['Instances'], store_phase.inputs['Geometry'])
    links.new(phase_attribute.outputs['Attribute'], store_phase.inputs['Value'])
    
    # Link color to store attribute
    links.new(store_phase.outputs['Geometry'], store_color.inputs['Geometry'])
    links.new(color_mix.outputs[2], store_color.inputs['Value'])  # Result output
    
    # Link to material and output
    links.new(store_color.outputs['Geometry'], set_material.inputs['Geometry'])
    links.new(set_material.outputs['Geometry'], output_node.inputs['Geometry'])
    
    print(f"✓ Geometry nodes setup complete")
    print(f"  Modifier: {modifier.name}")
    print(f"  Instance object reference: {instance_obj.name}")

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
    print("="*60)
    print("ORBITAL VISUALIZATION SCRIPT STARTING")
    print("="*60)
    
    # Clear existing scene
    clear_scene()
    print("✓ Scene cleared")
    
    # Create the base mesh with orbital points
    max_n = 4  # Generate orbitals up to n=4
    base_obj = create_base_mesh_with_points(max_n=max_n)
    
    if base_obj is None:
        print("="*60)
        print("✗ ERROR: Failed to create orbital mesh - no vertices generated!")
        print("This usually means scipy/numpy is not installed in Blender's Python")
        print("="*60)
        return
    
    print(f"✓ Base mesh created: {base_obj.name} with {len(base_obj.data.vertices)} vertices")
    
    # Create instance object (sphere)
    instance_obj = create_instance_sphere()
    print(f"✓ Instance sphere created: {instance_obj.name}")
    
    # Hide the instance object from viewport and render
    instance_obj.hide_viewport = True
    instance_obj.hide_render = True
    
    # Set up geometry nodes
    print("Setting up geometry nodes...")
    setup_geometry_nodes(base_obj, instance_obj)
    print("✓ Geometry nodes setup attempted")
    
    # If geometry nodes already exist from previous run, update the instance object reference
    for mod in base_obj.modifiers:
        if mod.type == 'NODES' and mod.node_group:
            for node in mod.node_group.nodes:
                if node.type == 'OBJECT_INFO':
                    node.inputs['Object'].default_value = instance_obj
                    print(f"✓ Updated Object Info reference in existing modifier: {mod.name}")
    
    # Add text labels
    add_text_labels()
    
    # Set up camera and lighting
    setup_camera_and_lighting()
    
    # Select the main orbital grid object
    bpy.context.view_layer.objects.active = base_obj
    base_obj.select_set(True)
    
    print("="*60)
    print(f"✓ ORBITAL VISUALIZATION COMPLETE")
    print(f"  Created n=1 to n={max_n}")
    print("  X axis: Angular momentum (l) - 0=s, 1=p, 2=d, 3=f")
    print("  Y axis: Magnetic quantum number (m)")
    print("  Z axis: Principal quantum number (n)")
    print(f"  Total vertices: {len(base_obj.data.vertices)}")
    print(f"  Modifiers on {base_obj.name}: {[m.name for m in base_obj.modifiers]}")
    print("="*60)

if __name__ == "__main__":
    main()
