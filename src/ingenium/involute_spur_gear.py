import math

import bpy
import numpy

from rustipy.option import Option, Some, Nothing

class InvoluteSpurGear:
    """
    Stores input parameters and calculates derived geometric properties for an involute spur gear.
    This class does not interact with Blender.
    """
    def __init__(self, module: float, num_teeth: int, pressure_angle_deg: float):
        """
        Initializes gear parameters and calculates derived values.

        Args:
            module (float): Gear module. Must be positive.
            num_teeth (int): Number of teeth. Must be 3 or greater.
            pressure_angle_deg (float): Pressure angle in degrees. Must be between 0 (exclusive) and 90 (exclusive).

        Raises:
            ValueError: If input parameters are invalid.
        """
        # --- Validate Inputs ---
        if module <= 0:
            raise ValueError("Module must be positive.")
        if num_teeth < 3:
            # Technically possible, but practically unusual and can lead to undercut issues
            # raise ValueError("Number of teeth must be 3 or greater for practical gears.")
            print("Warning: Number of teeth is less than 3, which is highly unusual.")
        if pressure_angle_deg <= 0 or pressure_angle_deg >= 90:
             raise ValueError("Pressure angle must be between 0 (exclusive) and 90 (exclusive) degrees.")

        # --- Store Base Parameters ---
        self.module: float = module
        self.num_teeth: int = num_teeth
        self.pressure_angle_deg: float = pressure_angle_deg
        self.pressure_angle_rad: float = numpy.radians(pressure_angle_deg)

        # --- Calculate Core Geometric Parameters ---
        # Pitch diameter and radius
        self.pitch_diameter: float = self.module * self.num_teeth
        self.pitch_radius: float = self.pitch_diameter / 2.0

        # Base circle radius
        self.base_radius: float = self.pitch_radius * numpy.cos(self.pressure_angle_rad)

        # Standard tooth dimensions (based on ISO 53 or similar standards)
        # Addendum (height from pitch circle to tip circle)
        self.addendum: float = 1.0 * self.module # Standard addendum factor ha* = 1.0
        # Clearance factor (for space at the bottom of the tooth gap)
        self.clearance_factor: float = 0.25 # Standard clearance factor c* = 0.25
        # Clearance (radial distance between root circle of one gear and tip circle of mating gear)
        self.clearance: float = self.clearance_factor * self.module
        # Dedendum (depth from pitch circle to root circle)
        # hf = ha + c = (ha* + c*) * m
        self.dedendum: float = self.addendum + self.clearance # Typically 1.25 * module

        # Tip circle radius (outer radius)
        self.tip_radius: float = self.pitch_radius + self.addendum # Ra = Rp + ha

        # Root circle radius (inner radius of tooth space)
        self.root_radius: float = self.pitch_radius - self.dedendum # Rf = Rp - hf
        # Ensure root radius is not negative, although this is unlikely with standard parameters
        if self.root_radius < 0:
             print(f"Warning: Calculated root radius ({self.root_radius:.4f}) is negative. Clamping to 0.")
             self.root_radius = 0.0


        # --- Calculate Max Involute Roll Angle (theta_max) ---
        # This is the roll angle 'theta' where the involute curve intersects the tip circle.
        self.theta_max: float = 0.0
        # Check if tip radius is significantly greater than base radius before sqrt
        # Use a small tolerance for floating point comparisons
        if self.tip_radius > self.base_radius + 1e-9:
            # Formula derived from R^2 = Rb^2 * (1 + theta^2)
            # => theta = sqrt((R/Rb)^2 - 1) where R is the radius of interest (tip_radius here)
            tip_radius_ratio_sq: float = (self.tip_radius / self.base_radius) ** 2
            # Ensure the argument to sqrt is non-negative due to potential floating point inaccuracies
            sqrt_arg: float = max(0.0, tip_radius_ratio_sq - 1.0)
            self.theta_max = numpy.sqrt(sqrt_arg)
        elif abs(self.tip_radius - self.base_radius) < 1e-9:
             # If tip radius is effectively equal to base radius, the involute starts and ends at the same point.
             self.theta_max = 0.0
        else:
             # This case (tip_radius < base_radius) should not happen with standard gear parameters
             # where addendum > 0 and pressure_angle < 90.
             print(f"Warning: Tip radius ({self.tip_radius:.4f}) is less than base radius ({self.base_radius:.4f}). This indicates unusual parameters or potential undercut issues. Setting theta_max to 0.")
             self.theta_max = 0.0

        # --- Calculate other useful angles (in radians) ---
        # Pitch angle: Angle subtended by one tooth and one space at the center
        self.pitch_angle: float = 2.0 * numpy.pi / self.num_teeth
        # Tooth angle: Angle subtended by one tooth at the pitch circle
        # Assuming standard tooth thickness = half the pitch angle at the pitch circle
        self.tooth_angle_at_pitch: float = self.pitch_angle / 2.0
        # Space angle: Angle subtended by one space at the pitch circle
        self.space_angle_at_pitch: float = self.pitch_angle / 2.0

        # Involute function (inv alpha) at the pressure angle
        self.inv_alpha: float = numpy.tan(self.pressure_angle_rad) - self.pressure_angle_rad

        # Angle phi: Angle from the start of the involute (on base circle) to the point where
        # the involute intersects the pitch circle. phi = tan(alpha) - alpha = inv(alpha)
        self.phi_at_pitch: float = self.inv_alpha

        # Angle delta: Half the tooth angle at the base circle.
        # This angle positions the center of the tooth flank relative to the starting point of the involute.
        # delta = (Tooth Angle at Pitch Circle / 2) + inv(alpha)
        # delta = (pi / (2 * N)) + inv(alpha)
        self.delta_at_base: float = (numpy.pi / (2.0 * self.num_teeth)) + self.inv_alpha

        print("--- Gear Parameters Initialized ---")
        print(f"  Module: {self.module}")
        print(f"  Num Teeth: {self.num_teeth}")
        print(f"  Pressure Angle: {self.pressure_angle_deg} deg")
        print(f"  Pitch Radius: {self.pitch_radius:.4f}")
        print(f"  Base Radius: {self.base_radius:.4f}")
        print(f"  Tip Radius: {self.tip_radius:.4f}")
        print(f"  Root Radius: {self.root_radius:.4f}")
        print(f"  Addendum: {self.addendum:.4f}")
        print(f"  Dedendum: {self.dedendum:.4f}")
        print(f"  Theta Max (Tip): {self.theta_max:.4f} rad")
        print(f"  Delta at Base: {self.delta_at_base:.4f} rad")
        print("------------------------------------")

    def guideline(self) -> 'GearGuidelineGenerator':
        return GearGuidelineGenerator(self)

    def mesh(self, name: str) -> 'GearMeshGenerator':
        return GearMeshGenerator(self)

class GearGuidelineGenerator:
    def __init__(self, gear: InvoluteSpurGear):
        self.gear = gear

    def _create_circle_guideline(self, radius: float, guideline_name: str) -> Option[bpy.types.Object]:
        """
        Creates a Bezier circle curve object using bpy.ops.
        Minimal implementation without exception handling or context management.

        Args:
            radius (float): The radius of the circle guideline.
            guideline_name (str): Name for the guideline object.

        Returns:
            Option[bpy.types.Object]: Some(Object) if successful, Nothing otherwise.
        """
        # 1. Validate radius
        if radius <= 1e-9:
            print(f"Error creating '{guideline_name}': Radius ({radius:.4f}) is too small or zero.")
            return Nothing()

        # 2. Create the circle using bpy.ops
        bpy.ops.curve.primitive_bezier_circle_add( # type: ignore , bpy Type Generator Error. 
            radius=radius,
            enter_editmode=False,
            align='WORLD',
            location=(0, 0, 0),
            scale=(1, 1, 1)
        )

        # 3. Get the newly created object (assumed to be the active one)
        created_obj = bpy.context.active_object
        if created_obj is None:
            print(f"Error: Failed to get active object after creating '{guideline_name}'.")
            return Nothing() # Return Nothing if object not found

        # 4. Rename the object
        created_obj.name = guideline_name
        if created_obj.data:
            created_obj.data.name = f"{guideline_name}_Data"

        print(f"Successfully created object '{guideline_name}' using bpy.ops.")
        
        return Some(created_obj) # Return the created object wrapped in Some

    def _create_line_guideline(self, start_point: tuple[float, float, float], end_point: tuple[float, float, float], guideline_name: str) -> Option[bpy.types.Object]:
        """
        Creates a simple straight line curve object using bpy.data.
        Minimal implementation.

        Args:
            start_point (tuple[float, float, float]): The starting coordinate (x, y, z).
            end_point (tuple[float, float, float]): The ending coordinate (x, y, z).
            guideline_name (str): Name for the guideline object.

        Returns:
            Option[bpy.types.Object]: Some(Object) if successful, Nothing otherwise.
        """
        curve_data = None
        try:
            # 1. Create new Curve Data
            curve_data = bpy.data.curves.new(name=f"{guideline_name}_Data", type='CURVE')
            curve_data.dimensions = '3D'

            # 2. Create a new Spline (Poly type)
            spline = curve_data.splines.new('POLY')

            # 3. Add points to the spline
            spline.points.add(1) # type: ignore , bpy Type Generator Error. # Add 1 point (total 2 points for a line segment)
            spline.points[0].co = (*start_point, 1.0) # W coordinate is 1.0 for Poly
            spline.points[1].co = (*end_point, 1.0)

            # 4. Create the Object using the Curve Data
            guideline_object = bpy.data.objects.new(guideline_name, curve_data)

            print(f"Successfully created line object data: '{guideline_name}'")
            return Some(guideline_object)

        except Exception as e:
            print(f"An unexpected error occurred during '{guideline_name}' line creation: {e}")
            if curve_data is not None and curve_data.name in bpy.data.curves:
                 bpy.data.curves.remove(curve_data) # type: ignore , bpy Type Generator Error. 
            return Nothing()

    def draw(self) -> list[bpy.types.Object]:
        """
        Creates and links the main circular guidelines to the current scene collection.
        Handles potential errors by returning only successfully linked objects.

        Returns:
            list[bpy.types.Object]: A list of the guideline objects successfully linked to the scene.
        """
        print("\n--- Drawing Circular Guidelines ---")
        linked_objects: list[bpy.types.Object] = []
        objects_to_cleanup: list[bpy.types.Object] = [] # Store objects created but not linked

        # 1. Define guidelines to create
        guidelines_spec = {
            "Base Circle Guideline": self.gear.base_radius,
            "Pitch Circle Guideline": self.gear.pitch_radius,
            "Root Circle Guideline": self.gear.root_radius,
            "Tip Circle Guideline": self.gear.tip_radius,
        }

        # 2. Create guideline objects
        created_options: dict[str, Option[bpy.types.Object]] = {
            name: self._create_circle_guideline(radius, name)
            for name, radius in guidelines_spec.items()
        }

        # --- Immediate Check after Creation ---
        # Check if the objects were automatically linked by bpy.ops
        collection_for_check = bpy.context.collection
        if collection_for_check is not None:
            print(f"DEBUG: Immediately checking collection '{collection_for_check.name}' after object creation...")
            for name, option in created_options.items():
                if option.is_some():
                    obj_ref = option.unwrap()
                    # Check only by name, as checking by reference causes TypeError
                    # Removed the try...except ReferenceError block and is_in_collection_by_ref
                    is_in_collection_by_name = name in collection_for_check.objects
                    print(f"DEBUG: Object '{name}' ({obj_ref.name}) - In collection by name? {is_in_collection_by_name}") # Updated print
        else:
            print("DEBUG: Cannot perform immediate check, no active collection.")
        # --- End Immediate Check ---

        # 3. Check for active collection (this might seem redundant now, but keep for safety)
        collection = bpy.context.collection
        if collection is None:
            print("Error: No active collection found. Cannot link guidelines.")
            # Add all successfully created objects to cleanup list
            for option in created_options.values():
                if option.is_some():
                    objects_to_cleanup.append(option.unwrap())
        else:
            # 4. Link/Verify successfully created objects if collection exists
            print(f"Attempting to link/verify guidelines in collection: '{collection.name}'") # Modified print
            for name, option in created_options.items():
                if option.is_some():
                    obj = option.unwrap()
                    # Check if the object (by name) is already in the collection
                    existing_obj_in_collection = collection.objects.get(name) # Use get()

                    if existing_obj_in_collection is not None:
                        # An object with this name exists in the collection.
                        # Check if it's the one we just created.
                        if existing_obj_in_collection == obj:
                            # Yes, it's the same object. bpy.ops likely auto-linked it.
                            print(f"Object '{name}' already linked (likely by bpy.ops). Ensuring visibility.")
                            existing_obj_in_collection.hide_viewport = False # Ensure visible
                            linked_objects.append(existing_obj_in_collection)
                            # DO NOT add to objects_to_cleanup
                        else:
                            # An OLD object with the same name exists. The new 'obj' is redundant.
                            print(f"An older object named '{name}' exists. Cleaning up the newly created duplicate.")
                            objects_to_cleanup.append(obj) # Add the NEW object 'obj' to cleanup
                    else:
                        # Object name is not in the collection, link the new object 'obj'
                        try:
                            print(f"Linking new object '{name}' to collection '{collection.name}'.")
                            collection.objects.link(obj) # type: ignore , bpy Type Generator Error. 
                            obj.hide_viewport = False # Ensure visible
                            linked_objects.append(obj)
                            print(f"Successfully linked '{name}'.")
                        except Exception as link_e:
                            print(f"Error linking '{name}': {link_e}")
                            objects_to_cleanup.append(obj) # Add to cleanup if linking failed
                # else: object creation failed, already printed in _create_circle_guideline

        # 5. Clean up objects that were created but not linked or were duplicates of OLD objects
        if objects_to_cleanup:
            print(f"Cleaning up {len(objects_to_cleanup)} unlinked or duplicate guideline objects...")
            for obj_clean in objects_to_cleanup:
                if obj_clean.name in bpy.data.objects:
                    # Check if the object reference is still valid and matches the one in bpy.data
                    if obj_clean == bpy.data.objects.get(obj_clean.name):
                        try:
                            bpy.data.objects.remove(obj_clean, do_unlink=True) # type: ignore , bpy Type Generator Error. 
                        except Exception as remove_e:
                            print(f"Error removing '{obj_clean.name}': {remove_e}")
                    else:
                         print(f"Warning: Object reference mismatch for '{obj_clean.name}' during cleanup. Skipping removal.")
                # else: Object already removed or never properly created

        print(f"--- Finished Drawing Circular Guidelines: {len(linked_objects)} linked. ---")
        return linked_objects # Return the list of successfully linked objects


class GearMeshGenerator:
    def __init__(self, gear: InvoluteSpurGear):
        self.gear = gear

    def draw(self):
        pass