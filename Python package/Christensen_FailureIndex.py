import numpy as np

class Christensen_class:
    """
    A class for evaluating the failure index and failure number based on the Christensen failure criterion.

    This class provides methods to compute the principal stresses, transform them into polar coordinates,
    and evaluate material failure based on the Christensen criterion. The failure assessment includes
    both the invariant-based failure criterion and the fracture criterion.

    Attributes:
    -----------
    T : float
        Tensile strength of the material.
    C : float
        Compressive strength of the material.
    S11, S22, S33, S12, S13, S23 : float
        Components of the stress tensor in either three-dimensional or plane strain states.
    Principal_Stresses : array
        Eigenvalues of the stress tensor representing the principal stresses.
    rho : float
        Radial coordinate in principal stress space.
    theta, phi : float
        Polar and azimuthal angles in principal stress space.
    rho_invariant_criterion : float
        Failure surface radius computed based on the Christensen invariant criterion.
    max_ps : float
        Maximum principal stress.
    failure_index : float
        Computed failure index, indicating the relative proximity to failure.
    Fn : float
        Failure number, describing the likelihood of failure.

    Methods:
    --------
    calc_main():
        Computes the failure index based on the Christensen failure criterion.
    calc_PC():
        Calculates the principal stresses from the stress tensor.
    calc_polar():
        Converts principal stresses into polar coordinates.
    calc_rho_invariant_criterion():
        Computes the failure surface radius based on the Christensen invariant criterion.
    calc_fracture_criterion():
        Determines the maximum principal stress, used in the fracture criterion.
    calc_failure_index():
        Computes the failure index based on the Christensen failure model.
    calc_fn():
        Calculates the failure number based on the material's failure state.
    """


    def __init__(self, entries_of_stress_tensor, T=100., C=300.):
        """
        Initializes a Christensen_class object with the given stress tensor components and material strengths.

        This constructor assigns the stress tensor components based on the number of provided entries
        and sets the material properties for tensile and compressive strength.

        Parameters:
        -----------
        entries_of_stress_tensor : list or array-like
            The stress tensor components, with:
            - 6 entries for a full three-dimensional stress state.
            - 4 entries for a plane strain state (S13 and S23 are set to zero).
            - 3 entries for a plane stress state (S33, S13, and S23 are set to zero).

        T : float, optional (default=100.)
            Tensile strength of the material.

        C : float, optional (default=300.)
            Compressive strength of the material.

        Raises:
        -------
        ValueError
            If the number of stress tensor components is not 3, 4, or 6.
        """
        
        # Assign material properties
        self.T = T
        self.C = C

        # Assign stress tensor entries depending on the stress state
        if len(entries_of_stress_tensor) == 6:
            self.S11 = entries_of_stress_tensor[0]
            self.S22 = entries_of_stress_tensor[1]
            self.S33 = entries_of_stress_tensor[2]
            self.S12 = entries_of_stress_tensor[3]
            self.S13 = entries_of_stress_tensor[4]
            self.S23 = entries_of_stress_tensor[5]
        elif len(entries_of_stress_tensor) == 4:
            self.S11 = entries_of_stress_tensor[0]
            self.S22 = entries_of_stress_tensor[1]
            self.S33 = entries_of_stress_tensor[2]
            self.S12 = entries_of_stress_tensor[3]
            self.S13 = 0.
            self.S23 = 0.
        elif len(entries_of_stress_tensor) == 3:
            self.S11 = entries_of_stress_tensor[0]
            self.S22 = entries_of_stress_tensor[1]
            self.S12 = entries_of_stress_tensor[2]
            self.S33 = 0.
            self.S13 = 0.
            self.S23 = 0.
        else:
            raise ValueError("The stress tensor has to have 3, 4 or 6 entries")

    def calc_main(self):
        """
        Computes the failure index based on the Christensen failure criterion.

        This method executes the full failure evaluation process by performing the following steps:
        1. Computes the principal stresses from the given stress tensor.
        2. Converts the principal stresses into polar coordinates.
        3. Evaluates the failure surface radius using the Christensen invariant criterion.
        4. Determines the maximum principal stress as part of the fracture criterion.
        5. Calculates the failure index based on the Christensen failure model.
        6. Computes the failure number, which quantifies the failure mode.

        Returns:
        --------
        list
            [failure_index, Fn] where failure_index is the material's failure state and Fn is the failure number.
        """

        # 1. Calculation of the Principal Stresses
        self.calc_PC()

        # 2. Calculation of the polar coordinates
        self.calc_polar()

        # 3. Calculation of the two sub-criteria
        if self.rho <= 0:
            return [0, 0]

        self.calc_rho_invariant_criterion()
        self.calc_fracture_criterion()

        # 4. calculation of the single failure index of the material
        self.calc_failure_index()

        self.calc_fn()
        # Return both failure index and Fn for plugin compatibility
        return [self.failure_index, self.Fn]

    # 1. Calculation of the Principal Stresses
    def calc_PC(self):
        """
        Computes the principal stresses of the given stress tensor.

        This method constructs the stress tensor from its components and calculates its eigenvalues,
        which correspond to the principal stresses.

        Attributes Updated:
        -------------------
        Principal_Stresses : numpy array
            The computed eigenvalues representing the principal stresses.

        Raises:
        -------
        ValueError
            If the stress tensor cannot be diagonalized due to numerical issues.
        """
        try:
            Stress_tensor = np.array([
                [self.S11, self.S12, self.S13],
                [self.S12, self.S22, self.S23],
                [self.S13, self.S23, self.S33]
            ])
        
            self.Principal_Stresses, eigenvec = np.linalg.eig(Stress_tensor)

        except np.linalg.LinAlgError as e:
            raise ValueError(f"Eigenvalue computation failed: {e}")


    def calc_polar(self):
        """
        Converts the principal stresses into polar coordinates in principal stress space.

        This method calculates the radial coordinate (rho) and the two angular coordinates (theta, phi)
        based on the principal stress values. The sine and cosine values of these angles are also stored
        for further computations.

        Attributes Updated:
        -------------------
        rho : float
            Radial coordinate representing the magnitude of the principal stress vector.
        theta : float
            Polar angle in the principal stress space.
        phi : float
            Azimuthal angle in the principal stress space.
        sin_t, cos_t : float
            Sine and cosine of the polar angle theta.
        sin_p, cos_p : float
            Sine and cosine of the azimuthal angle phi.

        Raises:
        -------
        ValueError
            If the calculation fails due to invalid principal stress values or numerical instabilities.
        """
        try:
            # Compute radial coordinate (rho)
            self.rho = np.sqrt(
                self.Principal_Stresses[0]**2 +
                self.Principal_Stresses[1]**2 +
                self.Principal_Stresses[2]**2
            )

            # Ensure rho is not zero before computing angles
            if self.rho == 0:
                raise ValueError("Polar coordinate computation failed: Stresses have to be non-zero")

            # Compute polar and azimuthal angles
            self.theta = np.arctan2(
                np.sqrt(self.Principal_Stresses[0]**2 + self.Principal_Stresses[1]**2),
                self.Principal_Stresses[2]
            )
            self.phi = np.arctan2(
                self.Principal_Stresses[1],
                self.Principal_Stresses[0]
            )

            # Compute trigonometric values
            self.sin_t = np.sin(self.theta)
            self.cos_t = np.cos(self.theta)
            self.sin_p = np.sin(self.phi)
            self.cos_p = np.cos(self.phi)

        except Exception as e:
            raise ValueError(f"Polar coordinate computation failed: {e}")

    # 3.1 Calculation of the Radius rho for the failure surface given the principal stresses and material strength
    def calc_rho_invariant_criterion(self):
        """
        Computes the failure surface radius (rho) for the failure criterion based on principal stresses
        and material strength using the Christensen invariant criterion.

        This method solves a quadratic equation for rho and selects the appropriate root based on the 
        principal stress state.

        Attributes Updated:
        -------------------
        rho_invariant_criterion : float
            The computed failure surface radius, used for determining the failure index.

        Raises:
        -------
        ValueError
            If the quadratic equation has no real solutions or the calculation encounters numerical instabilities.
        """
        try:
            aux_a_1 = 1. / (2. * self.T * self.C)
            aux_a_2 = (
                (self.sin_t * self.cos_p - self.sin_t * self.sin_p)**2 +
                (self.sin_t * self.cos_p - self.cos_t)**2 +
                (self.sin_t * self.sin_p - self.cos_t)**2
            )
            self.aux_a = aux_a_1 * aux_a_2
            self.aux_b = (1. / self.T - 1. / self.C) * (
                self.sin_t * self.sin_p + self.sin_t * self.cos_p + self.cos_t
            )
            self.aux_c = -1.

            # Handle numerical instability cases
            if abs(self.aux_a) < 1e-20:
                if sum(self.Principal_Stresses) > 0:
                    if self.aux_b != 0:
                        self.rho_invariant_criterion = -self.aux_c / self.aux_b
                    else:
                        raise ValueError("Computation failed: Division by zero in linear equation.")
                else:
                    raise ValueError("Computation failed: Undefined solution for rho invariant criterion.")
            else:
                discriminant = self.aux_b**2 - 4 * self.aux_a * self.aux_c

                if discriminant < 0:
                    raise ValueError("Computation failed: No real solution for rho invariant criterion.")

                rho_solutions = [
                    np.abs((-self.aux_b + np.sqrt(discriminant)) / (2 * self.aux_a)),
                    np.abs((-self.aux_b - np.sqrt(discriminant)) / (2 * self.aux_a))
                ]

                if sum(self.Principal_Stresses) > 0:
                    self.rho_invariant_criterion = min(rho_solutions)
                else:
                    self.rho_invariant_criterion = max(rho_solutions)

            return self.rho_invariant_criterion

        except Exception as e:
            raise ValueError(f"Failure criterion computation failed: {e}")

    # 3.2
    def calc_fracture_criterion(self):
        """
        Computes the fracture criterion based on the maximum principal stress.

        According to the Christensen failure theory, the fracture criterion is governed by the
        maximum principal stress. This method identifies the largest principal stress as
        a key factor in determining material failure.

        Attributes Updated:
        -------------------
        max_ps : float
            The maximum principal stress, which serves as a measure for the fracture criterion.
        """
        self.max_ps = max(self.Principal_Stresses)


    # 4 Calculation of the Failure Index
    def calc_failure_index(self):
        """
        Computes the failure index based on the Christensen failure criterion.

        The failure index quantifies the proximity of the given stress state to material failure.
        It is determined using both the invariant-based failure criterion and the fracture criterion.
        Depending on the ratio of tensile to compressive strength, different formulations apply.

        Attributes Updated:
        -------------------
        failure_index : float
            The computed failure index, which indicates whether the material is in a failed state.
            A value greater than or equal to 1 suggests material failure.


        Notes:
        ------
        - If `rho_invariant_criterion` is undefined (None), the failure index is set to NaN.
        - If `T / C < 1/2`, the maximum of two failure measures is used.
        - Otherwise, only the ratio of `rho` to `rho_invariant_criterion` determines the failure index.
        """
        if self.rho_invariant_criterion is not None:
            if self.T / self.C < 1/2:
                self.failure_index = max(self.rho / self.rho_invariant_criterion, self.max_ps / self.T)
            else:
                self.failure_index = self.rho / self.rho_invariant_criterion
        else:
            self.failure_index = float('NaN')



    # 5 Calculation of the Failure Number
    def calc_fn(self):
        """
        Computes the failure number (Fn) based on the Christensen failure criterion.

        The failure number is a measure of the stress states damage mechanism, ranging between 0 and 1.
        It is derived from the stress state at failure and the material's tensile and compressive strengths.

        Attributes Updated:
        -------------------
        Fn : float
            The computed failure number, where:
            - Fn = 0 indicates no damage.
            - Fn = 1 indicates full failure.
            - Intermediate values represent partial damage.

        Notes:
        ------
        - The failure number is computed using the stress state at failure.
        - If Fn exceeds 1, it is capped at 1.
        - If Fn is negative, it is set to 0.
        """
        principal_stresses_at_failure = [0., 0., 0.]
        principal_stresses_at_failure[0] = (self.rho / self.failure_index) * self.sin_t * self.cos_p
        principal_stresses_at_failure[1] = (self.rho / self.failure_index) * self.sin_t * self.sin_p
        principal_stresses_at_failure[2] = (self.rho / self.failure_index) * self.cos_t

        self.Fn = 1./2. * (3*self.T/self.C - (sum(principal_stresses_at_failure) / self.C))

        # Ensure Fn is within [0,1] range
        self.Fn = max(0, min(1, self.Fn))
