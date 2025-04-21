import math
import bpy

class InvoluteSpurGear:
    def __init__(self, module: float, num_teeth: int, pressure_angle_deg: float, fillet_radius_factor: float = 0.38):
        """
        Initializes the gear generator with basic parameters.

        Args:
            module (float): Gear module.
            num_teeth (int): Number of teeth. Must be >= 3.
            pressure_angle_deg (float): Pressure angle in degrees.
            fillet_radius_factor (float, optional): Factor of module for approximate root fillet radius. Defaults to 0.38.

        Raises:
            ValueError: If input parameters are invalid.
        """
        # --- Validate Inputs ---
        if num_teeth < 3:
            raise ValueError("Number of teeth must be 3 or greater.")
        if module <= 0:
            raise ValueError("Module must be positive.")
        if pressure_angle_deg <= 0 or pressure_angle_deg >= 90:
             raise ValueError("Pressure angle must be between 0 and 90 degrees.")

        # --- Store Base Parameters ---
        self.module: float = module
        self.num_teeth: int = num_teeth
        self.pressure_angle_deg: float = pressure_angle_deg
        self.pressure_angle_rad: float = math.radians(pressure_angle_deg)
        self.fillet_radius_factor: float = fillet_radius_factor

        # --- Calculate and Store Derived Geometric Parameters ---
        self.pitch_diameter: float = self.module * self.num_teeth
        self.pitch_radius: float = self.pitch_diameter / 2.0
        self.base_radius: float = self.pitch_radius * math.cos(self.pressure_angle_rad)

        if self.base_radius <= 1e-9:
             raise ValueError(
                 f"Base radius ({self.base_radius:.4f}) is too small or zero. "
                 f"Check parameters (num_teeth={num_teeth}, pressure_angle={pressure_angle_deg} deg)."
             )

        # Standard Addendum/Dedendum
        self.addendum: float = self.module
        self.dedendum: float = 1.25 * self.module # Example standard value
        self.tip_radius: float = self.pitch_radius + self.addendum
        self.root_radius: float = self.pitch_radius - self.dedendum

        # Effective root radius for fillet calculation (cannot be below base radius)
        self.effective_root_radius: float = max(self.root_radius, self.base_radius)
        if self.root_radius < self.base_radius:
            print(f"Warning: Calculated root radius {self.root_radius:.4f} is below base radius {self.base_radius:.4f}.")
            # print("         Involute profile will start from base circle. Root fillet calculation needs adjustment.") # Reduced verbosity

        # Approximate fillet radius
        self.fillet_radius: float = self.module * self.fillet_radius_factor

        # Angular parameters (radians)
        self.tooth_pitch_angle: float = 2 * math.pi / self.num_teeth # Angle covering one tooth and one space
        self.tooth_angle_on_pitch: float = self.tooth_pitch_angle / 2.0 # Approx angle of tooth thickness at pitch circle
        self.space_angle_on_pitch: float = self.tooth_pitch_angle / 2.0 # Approx angle of space width at pitch circle

        # Max involute roll angle to reach the tip radius
        self.theta_max: float = 0.0
        if self.tip_radius > self.base_radius:
            argument = max(0.0, (self.tip_radius / self.base_radius)**2 - 1)
            self.theta_max = math.sqrt(argument)

        print(f"InvoluteSpurGearGenerator initialized:")
        print(f"  Rb={self.base_radius:.4f}, Rp={self.pitch_radius:.4f}, Ra={self.tip_radius:.4f}, Eff.Rf={self.effective_root_radius:.4f}")
        # print(f"  Approx Fillet Radius={self.fillet_radius:.4f}, Theta Max={math.degrees(self.theta_max):.2f} deg")
        # print(f"  Tooth Pitch Angle={math.degrees(self.tooth_pitch_angle):.2f} deg")


    def _calculate_involute_point(self, theta: float) -> tuple[float, float, float]:
        """Calculates a single point (vertex) on the involute curve for roll angle theta."""
        if theta < -1e-9:
             # print(f"Warning: Negative theta ({theta}) encountered in involute calculation.")
             theta = 0.0 # Clamp

        x = self.base_radius * (math.cos(theta) + theta * math.sin(theta))
        y = self.base_radius * (math.sin(theta) - theta * math.cos(theta))
        
        return (x, y, 0.0)

    def _get_involute_profile_points(self, num_points: int) -> list[tuple[float, float, float]]:
        """Generates points for one standard involute flank (base to tip)."""
        if num_points < 1:
             return [(self.base_radius, 0.0, 0.0)] if self.base_radius > 0 else []
        if self.theta_max <= 1e-9:
             return [(self.base_radius, 0.0, 0.0)] if self.base_radius > 0 else []

        involute_points = []
        for i in range(num_points + 1):
            theta = (i / num_points) * self.theta_max
            point = self._calculate_involute_point(theta)
            involute_points.append(point)
            
        return involute_points

    def _rotate_point(self, point: tuple[float, float, float], angle_rad: float) -> tuple[float, float, float]:
        """Rotates a 2D point (x, y, z) around the Z axis by angle_rad."""
        x, y, z = point
        cos_a = math.cos(angle_rad)
        sin_a = math.sin(angle_rad)
        new_x = x * cos_a - y * sin_a
        new_y = x * sin_a + y * cos_a
        
        return (new_x, new_y, z)

    def _calculate_circular_root_fillet_points(
        self,
        start_point_approx: tuple[float, float, float],
        end_point_approx: tuple[float, float, float],
        slot_center_angle_rad: float,
        num_points: int
    ) -> list[tuple[float, float, float]]:
        """
        (PLACEHOLDER - Highly Simplified Approximation)
        Generates points for a circular arc fillet between two approximate points.
        NOTE: This is NOT mathematically correct for tangent blending.
        """
        fillet_points = []
        # --- Very Basic Placeholder Logic ---
        center_x = self.effective_root_radius * math.cos(slot_center_angle_rad)
        center_y = self.effective_root_radius * math.sin(slot_center_angle_rad)
        start_vec_x = start_point_approx[0] - center_x
        start_vec_y = start_point_approx[1] - center_y
        start_angle = math.atan2(start_vec_y, start_vec_x)
        end_vec_x = end_point_approx[0] - center_x
        end_vec_y = end_point_approx[1] - center_y
        end_angle = math.atan2(end_vec_y, end_vec_x)

        if end_angle <= start_angle: end_angle += 2 * math.pi

        if num_points > 0:
            for i in range(1, num_points + 1): # Exclude start point, include end point
                interp_factor = i / num_points
                current_angle = start_angle + (end_angle - start_angle) * interp_factor
                x = center_x + self.fillet_radius * math.cos(current_angle)
                y = center_y + self.fillet_radius * math.sin(current_angle)
                fillet_points.append((x, y, 0.0))

        if not fillet_points: # Ensure at least end point is returned
             fillet_points.append(end_point_approx)

        return fillet_points

    # --- Public Interface Methods ---

    def generate_single_tooth_slot_profile(
        self,
        slot_center_angle_rad: float = 0.0, # Center angle of the space BETWEEN teeth
        num_profile_points: int = 20,
        num_fillet_points: int = 10
    ) -> list[tuple[float, float, float]]:
        """
        Generates vertices for a single tooth slot profile (one fillet + two flanks).
        Points are ordered roughly from Tip of preceding flank -> Root Fillet -> Tip of succeeding flank.
        The profile is generated locally and then rotated by slot_center_angle_rad.
        """
        # 1. Calculate base involute profile (Flank A - local coords, starts on X-axis)
        involute_A_local = self._get_involute_profile_points(num_points=num_profile_points)
        if not involute_A_local or len(involute_A_local) < 2:
            print("Error: Failed to generate involute profile A.")
            
            return []

        # 2. Determine rotation angles for flanks relative to slot center
        #    Calculation based on involute properties at base circle
        inv_alpha = math.tan(self.pressure_angle_rad) - self.pressure_angle_rad
        base_angle_shift = self.tooth_angle_on_pitch / 2.0 - inv_alpha # Angle from flank center to involute start

        # Angle to rotate Flank A (starts locally on X-axis) to position it correctly
        # It forms the *trailing* side of the tooth space (leading side of the tooth)
        rotate_A = -self.space_angle_on_pitch / 2.0 + base_angle_shift # CCW from space edge

        # Angle to rotate mirrored Flank B (forms the *leading* side of the space)
        rotate_B = self.space_angle_on_pitch / 2.0 - base_angle_shift # CW from space edge

        # 3. Create Rotated Flank A (Root -> Tip)
        involute_A_rotated = [self._rotate_point(p, rotate_A) for p in involute_A_local]

        # 4. Create Rotated Flank B (Tip -> Root)
        involute_B_local_mirrored = [(p[0], -p[1], p[2]) for p in involute_A_local] # Mirror across local X
        involute_B_rotated = [self._rotate_point(p, rotate_B) for p in involute_B_local_mirrored]
        involute_B_final = list(reversed(involute_B_rotated)) # Reverse order to Tip -> Root

        # 5. Calculate Fillet points (using placeholder function)
        fillet_start_point_approx = involute_B_final[-1] # End of B (at root)
        fillet_end_point_approx = involute_A_rotated[0]  # Start of A (at root)

        fillet_points = self._calculate_circular_root_fillet_points(
            start_point_approx=fillet_start_point_approx,
            end_point_approx=fillet_end_point_approx,
            slot_center_angle_rad=0.0, # Fillet calculated relative to local slot center (0)
            num_points=num_fillet_points
        )
        if not fillet_points: return []

        # 6. Combine lists locally: Flank B (Tip->Root) + Fillet Points + Flank A (Root->Tip)
        combined_profile_local = []
        if involute_B_final: combined_profile_local.extend(involute_B_final)

        # Add fillet points, avoiding duplicate start point if possible (using approx comparison)
        if fillet_points:
            if combined_profile_local and \
               math.dist(fillet_points[0][:2], combined_profile_local[-1][:2]) < 1e-6:
                combined_profile_local.extend(fillet_points[1:])
            else:
                combined_profile_local.extend(fillet_points)

        # Add flank A points, avoiding duplicate start point
        if involute_A_rotated:
            if combined_profile_local and \
               math.dist(involute_A_rotated[0][:2], combined_profile_local[-1][:2]) < 1e-6:
                combined_profile_local.extend(involute_A_rotated[1:])
            else:
                combined_profile_local.extend(involute_A_rotated)

        # 7. Rotate the entire profile to the specified global slot_center_angle_rad
        final_profile = [self._rotate_point(p, slot_center_angle_rad) for p in combined_profile_local]

        return final_profile

    def generate_full_gear_vertices(
        self,
        num_profile_points: int = 20,
        num_fillet_points: int = 10
    ) -> list[tuple[float, float, float]]:
        """
        Generates vertices for the complete 2D gear profile by combining all tooth slots.

        Args:
            num_profile_points (int, optional): Resolution for the involute flank. Defaults to 20.
            num_fillet_points (int, optional): Resolution for the root fillet arc. Defaults to 10.

        Returns:
            list[tuple[float, float, float]]: A list of (x, y, z) coordinates representing the gear outline.
                                             The list forms a closed loop (ideally).
        """
        print(f"\nGenerating full gear profile vertices (N={self.num_teeth})...")
        all_vertices = []
        for i in range(self.num_teeth):
            # Calculate the center angle for this tooth space
            # We align the *space* center, assuming tooth thickness = space width at pitch
            slot_center_angle = i * self.tooth_pitch_angle

            # Generate the profile for this slot
            slot_profile = self.generate_single_tooth_slot_profile(
                slot_center_angle_rad=slot_center_angle,
                num_profile_points=num_profile_points,
                num_fillet_points=num_fillet_points
            )

            # Add points, handling potential overlap between the end of the previous
            # profile and the start of the current one.
            if not slot_profile:
                 print(f"Warning: Empty profile generated for slot {i}. Skipping.")
                 continue

            if all_vertices and math.dist(all_vertices[-1][:2], slot_profile[0][:2]) < 1e-6:
                 all_vertices.extend(slot_profile[1:]) # Skip first point if it's same as last
            else:
                 all_vertices.extend(slot_profile)

        # Optional: Check if the loop is closed (last point == first point)
        if all_vertices and math.dist(all_vertices[-1][:2], all_vertices[0][:2]) < 1e-6:
             print("  Loop closed successfully (last point matches first point).")
             # Optionally remove the last point if it's a duplicate of the first
             # all_vertices.pop()
        elif all_vertices:
             print("Warning: Loop closure check failed. Last point may not match first point exactly.")


        print(f"Full gear profile generated with {len(all_vertices)} vertices.")
        
        return all_vertices
    
    def to_mesh(self, name, module, num_teeth, pressure_angle_deg):
        gear = InvoluteSpurGear(module, num_teeth, pressure_angle_deg)
        print("-" * 20)

        # 2. Generate vertices for the full 2D profile
        num_profile_points = 20
        num_fillet_points = 10
        gear_vertices = gear.generate_full_gear_vertices(
            num_profile_points=num_profile_points,
            num_fillet_points=num_fillet_points
        )

        mesh_name = "GeneratedGearMesh"
        mesh_data = bpy.data.meshes.new(mesh_name)

        # Create Edges for the outline loop
        num_verts = len(gear_vertices)
        edges = [[i, (i + 1) % num_verts] for i in range(num_verts)]

        # Create Face (single N-gon face for the 2D profile)
        faces = [list(range(num_verts))]

        # Populate mesh data
        mesh_data.from_pydata(gear_vertices, edges, faces)
        mesh_data.update() # Ensure mesh data is valid

        # Create Object
        object_name = name
        gear_object = bpy.data.objects.new(object_name, mesh_data)

        # Link Object to the current scene's collection
        bpy.context.collection.objects.link(gear_object)

        # (Optional) Make the new object active and selected
        bpy.context.view_layer.objects.active = gear_object
        gear_object.select_set(True)

        print(f"Successfully created Blender object: '{object_name}'")
        # --- End of BPY Code ---