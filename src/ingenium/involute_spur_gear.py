import math
import bpy

class InvoluteSpurGear:
    def __init__(self, module: float, num_teeth: int, pressure_angle_deg: float):
        """
        Initializes the gear generator with basic parameters.

        Args:
            module (float): Gear module (e.g., in mm).
            num_teeth (int): Number of teeth. Must be >= 3.
            pressure_angle_deg (float): Pressure angle in degrees (e.g., 20.0).

        Raises:
            ValueError: If input parameters are invalid.
        """
        # --- Validate Inputs ---
        if num_teeth < 3:
            # Ensure a minimum number of teeth for valid geometry
            raise ValueError("Number of teeth must be 3 or greater.")
        if module <= 0:
            # Module must be a positive value
            raise ValueError("Module must be positive.")
        if pressure_angle_deg <= 0 or pressure_angle_deg >= 90:
             # Pressure angle must be within a reasonable range
             raise ValueError("Pressure angle must be between 0 and 90 degrees.")

        # --- Store Base Parameters ---
        self.module: float = module
        self.num_teeth: int = num_teeth
        self.pressure_angle_deg: float = pressure_angle_deg
        # Convert pressure angle to radians for trigonometric functions
        self.pressure_angle_rad: float = math.radians(pressure_angle_deg)

        # --- Calculate Core Geometric Parameters ---

        # Pitch circle calculations: The fundamental circle defining gear size
        self.pitch_diameter: float = self.module * self.num_teeth
        self.pitch_radius: float = self.pitch_diameter / 2.0

        # Base circle radius: The circle from which the involute curve is generated
        self.base_radius: float = self.pitch_radius * math.cos(self.pressure_angle_rad)
        # Check if base radius is valid (can become near zero with extreme parameters)
        if self.base_radius <= 1e-9:
            raise ValueError(
                f"Base radius ({self.base_radius:.4f}) is too small or zero. "
                f"Check parameters (num_teeth={num_teeth}, pressure_angle={pressure_angle_deg} deg)."
            )

        # Standard tooth dimensions (based on module)
        # These define the radial extent of the teeth
        # Note: These are standard values; non-standard gears might use different factors.
        self.addendum: float = 1.0 * self.module # Radial distance from pitch circle to tooth tip
        self.dedendum: float = 1.25 * self.module # Radial distance from pitch circle to nominal root

        # Tip circle radius: The outermost radius of the gear teeth
        self.tip_radius: float = self.pitch_radius + self.addendum
        # Root circle radius: The nominal radius at the bottom of the tooth space
        self.root_radius: float = self.pitch_radius - self.dedendum

        # Check for potential undercut condition
        if self.root_radius < self.base_radius:
            # This warning indicates that the theoretical root circle lies inside the base circle.
            # The involute profile cannot exist below the base circle.
            # The actual root fillet (trochoid) will start from the base circle.
            # This condition often implies undercut during manufacturing.
            print(f"Warning: Calculated root radius {self.root_radius:.4f} is below base radius {self.base_radius:.4f}. Potential undercut.")

        # --- Calculate Angular Parameters ---

        # Angle subtended by one tooth and one space at the gear center
        self.tooth_pitch_angle: float = (2 * math.pi) / self.num_teeth
        # Angular thickness of a tooth on the pitch circle
        # Standard assumption: tooth thickness equals space width on the pitch circle
        self.tooth_angle_on_pitch: float = self.tooth_pitch_angle / 2.0
        # Angular width of a space between teeth on the pitch circle
        self.space_angle_on_pitch: float = self.tooth_pitch_angle / 2.0

        # --- Calculate Involute Generation Angle Limit ---

        # Maximum roll angle (theta) required to generate the involute curve
        # up to the tip radius. Derived from the involute equation:
        # radius = base_radius * sqrt(1 + theta^2)
        self.theta_max: float = 0.0 # Initialize theta_max
        # Ensure base_radius is positive to avoid division by zero or math domain error
        if self.base_radius > 1e-9:
             # Calculate the squared ratio of tip radius to base radius
             tip_radius_ratio_sq = (self.tip_radius / self.base_radius) ** 2
             # Check if tip circle is outside base circle (it should be)
             if tip_radius_ratio_sq < (1.0 - 1e-9): # Use tolerance for float comparison
                 # This indicates an error in parameters or calculations
                 print(f"Error: Tip radius {self.tip_radius:.4f} seems smaller than base radius {self.base_radius:.4f}. Cannot calculate theta_max.")
                 # Optionally, raise an error instead of just printing
                 # raise ValueError("Invalid geometry: Tip radius is smaller than base radius.")
             else:
                 # Ensure the argument for sqrt is non-negative
                 sqrt_arg = max(0.0, tip_radius_ratio_sq - 1.0)
                 self.theta_max = math.sqrt(sqrt_arg)
        # else: base_radius is near zero, theta_max remains 0 (already handled by earlier check)

        # --- Debug Prints (Optional) ---
        # print(f"Debug: module={self.module}, num_teeth={self.num_teeth}, pressure_angle={self.pressure_angle_deg}")
        # print(f"Debug: pitch_radius={self.pitch_radius:.4f}, base_radius={self.base_radius:.4f}")
        # print(f"Debug: tip_radius={self.tip_radius:.4f}, root_radius={self.root_radius:.4f}")
        # print(f"Debug: theta_max={math.degrees(self.theta_max):.2f} deg")


    def _calculate_involute_point(self, theta: float) -> tuple[float, float, float]:
        """
        Calculates the (x, y, z) coordinates of a point on the involute curve
        for a given roll angle theta, relative to the base circle center.
        The curve starts at (base_radius, 0) for theta = 0.

        Args:
            theta (float): The roll angle in radians (angle unwrapped from base circle).
                           Must be non-negative.

        Returns:
            tuple[float, float, float]: The (x, y, z) coordinates of the point. z is always 0.

        Raises:
            ValueError: If theta is negative.
        """
        if theta < -1e-9: # Allow for small floating point inaccuracies near zero
            raise ValueError(f"Theta ({theta}) cannot be negative for involute calculation.")

        # Clamp theta near zero if it's very slightly negative due to float issues
        theta = max(0.0, theta)

        # Standard involute parametric equations
        x = self.base_radius * (math.cos(theta) + theta * math.sin(theta))
        y = self.base_radius * (math.sin(theta) - theta * math.cos(theta))
        return (x, y, 0.0)

    def _get_involute_profile_points(self, num_points: int) -> list[tuple[float, float, float]]:
        """
        Generates a list of points defining one side of the tooth profile
        using the involute curve, from the base circle up to the tip circle.

        Args:
            num_points (int): The number of points to generate along the curve segment
                              (excluding the start point, so num_points=10 gives 11 points total).
                              Must be >= 1.

        Returns:
            list[tuple[float, float, float]]: A list of (x, y, z) points defining the involute profile.
                                             The list starts at the base circle (theta=0) and ends
                                             at or near the tip circle (theta=theta_max).
                                             Returns an empty list if num_points < 1 or theta_max is invalid.
        """
        if num_points < 1:
            print("Warning: num_points must be >= 1 to generate involute profile. Returning empty list.")
            return []
        if self.theta_max <= 1e-9: # Check if theta_max is valid (calculated in __init__)
             print(f"Warning: theta_max ({self.theta_max:.4f}) is too small. Cannot generate involute profile. Returning empty list.")
             return []

        # Generate theta values from 0 to theta_max using a generator expression
        # We need num_points + 1 total points (including start and end)
        thetas = ((i / num_points) * self.theta_max for i in range(num_points + 1))

        # Calculate involute points for each theta using another generator expression
        # and collect them into a list
        involute_points: list[tuple[float, float, float]] = [
            self._calculate_involute_point(theta) for theta in thetas
        ]

        return involute_points

    def _rotate_point(self, point: tuple[float, float, float], angle_rad: float) -> tuple[float, float, float]:
        """Rotates a 2D point (x, y) around the origin by a given angle."""
        x, y, z = point
        cos_a = math.cos(angle_rad)
        sin_a = math.sin(angle_rad)
        new_x = x * cos_a - y * sin_a
        new_y = x * sin_a + y * cos_a
        return (new_x, new_y, z)
    
    def _generate_single_tooth_slot_profile(
        self,
        slot_center_angle_rad: float, # Center angle of the space BETWEEN teeth
        num_profile_points: int # Number of points for EACH involute flank
    ) -> list[tuple[float, float, float]]:
        """
        Generates vertices for a single tooth slot profile, consisting of two
        involute flanks. The root fillet points are NOT generated by this method
        and should be calculated and inserted separately (e.g., using trochoid).

        Points are ordered roughly from Tip of preceding flank -> Root of preceding flank,
        then Root of succeeding flank -> Tip of succeeding flank.
        The profile is generated locally around the X-axis and then rotated by slot_center_angle_rad.

        Args:
            slot_center_angle_rad (float): The global angle (in radians) to center this slot profile.
            num_profile_points (int): Number of points to generate for each involute flank segment
                                      (must be >= 1).

        Returns:
            list[tuple[float, float, float]]: A list of (x, y, z) points defining the two flanks
                                             of the tooth slot. Returns empty list on error.
        """
        if num_profile_points < 1:
             print("Error: num_profile_points must be >= 1. Cannot generate slot profile.")
             return []

        # 1. Calculate base involute profile (Flank A - local coords, starts on X-axis at base circle)
        #    Points are ordered Root -> Tip
        involute_A_local = self._get_involute_profile_points(num_points=num_profile_points)
        if not involute_A_local or len(involute_A_local) < 2:
            # Error message already printed by _get_involute_profile_points if needed
            print("Error: Failed to generate base involute profile for slot.")
            return []

        # 2. Determine rotation angles for flanks relative to the slot center (X-axis locally)
        #    Calculation based on involute properties at the pitch circle to ensure correct tooth thickness
        #    inv_alpha: Angle between radius to involute start (on base circle) and radius to point on pitch circle
        inv_pressure_angle = math.tan(self.pressure_angle_rad) - self.pressure_angle_rad
        #    half_tooth_angle_at_pitch: Half the angular thickness of the tooth on the pitch circle
        half_tooth_angle_at_pitch = self.tooth_angle_on_pitch / 2.0
        #    base_angle_shift: Angle from the center of the tooth flank (on pitch circle)
        #                      to the start of the involute curve (on base circle)
        base_angle_shift = half_tooth_angle_at_pitch + inv_pressure_angle # Corrected calculation

        # Angle to rotate Flank A (starts locally on X-axis) to position it correctly.
        # Flank A forms the *leading* side of the tooth (trailing side of the space).
        # It needs to be rotated counter-clockwise from the slot center.
        rotate_A = base_angle_shift

        # Angle to rotate mirrored Flank B.
        # Flank B forms the *trailing* side of the tooth (leading side of the space).
        # It needs to be rotated clockwise from the slot center.
        rotate_B = -base_angle_shift # Mirrored rotation

        # 3. Create Rotated Flank A (Root -> Tip)
        flank_A_rotated = [self._rotate_point(p, rotate_A) for p in involute_A_local]

        # 4. Create Rotated Flank B (Tip -> Root)
        #    Mirror the original local profile across the X-axis
        involute_B_local_mirrored = [(p[0], -p[1], p[2]) for p in involute_A_local]
        #    Rotate the mirrored profile
        flank_B_rotated = [self._rotate_point(p, rotate_B) for p in involute_B_local_mirrored]
        #    Reverse the order so it goes from Tip -> Root
        flank_B_final = list(reversed(flank_B_rotated))

        # 5. Combine Flanks (Fillet points are missing here!)
        #    The order is: Flank B (Tip -> Root) followed by Flank A (Root -> Tip)
        #    A gap exists between flank_B_final[-1] and flank_A_rotated[0] where the fillet should be.
        combined_profile_local: list[tuple[float, float, float]] = []
        combined_profile_local.extend(flank_B_final)
        # --- Placeholder for Trochoid Fillet Insertion ---
        # fillet_points = self._calculate_trochoid_fillet_points(
        #     start_flank_point=flank_B_final[-1], # End of Flank B (Root)
        #     end_flank_point=flank_A_rotated[0],  # Start of Flank A (Root)
        #     num_points=...
        # )
        # combined_profile_local.extend(fillet_points) # Add fillet points here
        # -------------------------------------------------
        combined_profile_local.extend(flank_A_rotated)

        # 6. Rotate the entire profile to the specified global slot_center_angle_rad
        final_profile = [self._rotate_point(p, slot_center_angle_rad) for p in combined_profile_local]

        return final_profile

    # ... rest of the class ...

    # def _calculate_circular_root_fillet_points(...): # This method is now obsolete/removed

    # ... (Need to add _calculate_trochoid_fillet_points later) ...

    def generate_gear_points(self, num_profile_points_per_flank: int = 20) -> list[tuple[float, float, float]]:
        """
        Generates all vertex points for the gear profile by combining multiple
        tooth slot profiles.

        Args:
            num_profile_points_per_flank (int): Number of points for each involute flank segment.

        Returns:
            list[tuple[float, float, float]]: A list of (x, y, z) points for the complete gear outline.
        """
        all_gear_points: list[tuple[float, float, float]] = []
        for i in range(self.num_teeth):
            # Calculate the center angle for each tooth *space*
            slot_center_angle = i * self.tooth_pitch_angle
            single_slot_profile = self._generate_single_tooth_slot_profile(
                slot_center_angle_rad=slot_center_angle,
                num_profile_points=num_profile_points_per_flank
            )
            if single_slot_profile:
                # Avoid duplicating the connection point between slots if possible
                if all_gear_points and single_slot_profile:
                     # Check if the start of the new profile is close to the end of the previous one
                     if math.dist(all_gear_points[-1][:2], single_slot_profile[0][:2]) < 1e-6:
                         all_gear_points.extend(single_slot_profile[1:]) # Skip duplicate point
                     else:
                         all_gear_points.extend(single_slot_profile)
                else:
                    all_gear_points.extend(single_slot_profile)
            else:
                print(f"Warning: Failed to generate profile for slot {i}. Skipping.")

        # Close the loop by connecting the last point to the first point if they are not the same
        if all_gear_points and len(all_gear_points) > 1:
             if math.dist(all_gear_points[-1][:2], all_gear_points[0][:2]) > 1e-6:
                 # This might indicate an issue if the profile isn't naturally closing
                 print("Warning: Gear profile loop doesn't seem closed. Check calculations.")
                 # Optionally, add the first point again to force closure, but it's better to fix the root cause
                 # all_gear_points.append(all_gear_points[0])
             else:
                 # Remove the last point if it's identical to the first (avoids duplicate vertex)
                 all_gear_points.pop()


        return all_gear_points

    def to_mesh(self, name: str, module: float, num_teeth: int, pressure_angle_deg: float):
        gear = InvoluteSpurGear(module, num_teeth, pressure_angle_deg)
        print("-" * 20)

        # 2. Generate vertices for the full 2D profile
        num_profile_points = 20
        num_fillet_points = 10
        gear_vertices = gear.generate_gear_points(
            num_profile_points_per_flank=num_profile_points
        )

        mesh_name = "GeneratedGearMesh"
        mesh_data = bpy.data.meshes.new(mesh_name)

        num_verts = len(gear_vertices)
        edges = [[i, (i + 1) % num_verts] for i in range(num_verts)]
        faces = [list(range(num_verts))]

        mesh_data.from_pydata(gear_vertices, edges, faces)
        mesh_data.update() # Ensure mesh data is valid

        object_name = name
        gear_object = bpy.data.objects.new(object_name, mesh_data)

        bpy.context.collection.objects.link(gear_object) # error due to .collection variable could be None
        bpy.context.view_layer.objects.active = gear_object # error due to .view_layer variable could be None

        gear_object.select_set(True)

        print(f"Successfully created Blender object: '{object_name}'")

