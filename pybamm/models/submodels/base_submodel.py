#
# Base submodel class
#
import pybamm


class BaseSubModel(pybamm.BaseModel):
    """
    The base class for all submodels. All submodels inherit from this class and must
    only provide public methods which overwrite those in this base class. Any methods
    added to a submodel that do not overwrite those in this bass class are made
    private with the prefix '_', providing a consistent public interface for all
    submodels.

    Parameters
    ----------
    param: parameter class
        The model parameter symbols
    domain : str
        The domain of the model either 'Negative' or 'Positive'
    name: str
        A string giving the name of the submodel
    external: bool, optional
        Whether the variables defined by the submodel will be provided externally
        by the users. Default is 'False'.
    options: dict
        A dictionary of options to be passed to the model.
        See :class:`pybamm.BaseBatteryModel`
    phase : str, optional
        Phase of the particle (default is None).

    Attributes
    ----------
    param: parameter class
        The model parameter symbols
    rhs: dict
        A dictionary that maps expressions (variables) to expressions that represent
        the rhs
    algebraic: dict
        A dictionary that maps expressions (variables) to expressions that represent
        the algebraic equations. The algebraic expressions are assumed to equate
        to zero. Note that all the variables in the model must exist in the keys of
        `rhs` or `algebraic`.
    initial_conditions: dict
        A dictionary that maps expressions (variables) to expressions that represent
        the initial conditions for the state variables y. The initial conditions for
        algebraic variables are provided as initial guesses to a root finding algorithm
        that calculates consistent initial conditions.
    boundary_conditions: dict
        A dictionary that maps expressions (variables) to expressions that represent
        the boundary conditions
    variables: dict
        A dictionary that maps strings to expressions that represent
        the useful variables
    events: list
        A list of events. Each event can either cause the solver to terminate
        (e.g. concentration goes negative), or be used to inform the solver of the
        existance of a discontinuity (e.g. discontinuity in the input current)
    """

    def __init__(
        self,
        param,
        domain=None,
        name="Unnamed submodel",
        external=False,
        options=None,
        phase=None,
    ):
        super().__init__(name)
        self.domain = domain
        self.set_domain_for_broadcast()
        self.name = name

        self.external = external
        self.options = pybamm.BatteryModelOptions(options or {})

        # Save whether the submodel is a half-cell submodel
        we = self.options["working electrode"]
        self.half_cell = we != "both"

        self.param = param
        if param is None:
            self.domain_param = None
        else:
            if self.domain == "Negative":
                self.domain_param = param.n
            elif self.domain == "Positive":
                self.domain_param = param.p

            if phase == "primary":
                self.phase_param = self.domain_param.prim
            elif phase == "secondary":
                self.phase_param = self.domain_param.sec

        # Error checks for phase and domain
        self.set_phase(phase)

    def set_phase(self, phase):
        if phase is not None:
            if self.domain is None:
                raise ValueError("Phase must be None if domain is None")
            options_phase = getattr(self.options, self.domain.lower())[
                "particle phases"
            ]
            if options_phase == "1" and phase != "primary":
                raise ValueError("Phase must be 'primary' if there is only one phase")
            elif options_phase == "2" and phase not in ["primary", "secondary"]:
                raise ValueError(
                    "Phase must be either 'primary' or 'secondary' "
                    "if there are two phases"
                )

            if options_phase == "1" and phase == "primary":
                # Only one phase, no need to distinguish between
                # "primary" and "secondary"
                self.phase_name = ""
            else:
                # add a space so that we can use "" or (e.g.) "primary " interchangeably
                # when naming variables
                self.phase_name = phase + " "

        self.phase = phase

    @property
    def domain(self):
        return self._domain

    @domain.setter
    def domain(self, domain):
        if domain is not None:
            domain = domain.capitalize()
        ok_domain_list = [
            "Negative",
            "Separator",
            "Positive",
            "Negative electrode",
            "Negative electrolyte",
            "Separator electrolyte",
            "Positive electrode",
            "Positive electrolyte",
            None,
        ]
        if domain in ok_domain_list:
            self._domain = domain
        else:
            raise pybamm.DomainError(
                "Domain '{}' not recognised (must be one of {})".format(
                    domain, ok_domain_list
                )
            )

    def set_domain_for_broadcast(self):
        if hasattr(self, "_domain"):
            if self.domain in ["Negative", "Positive"]:
                self.domain_for_broadcast = self.domain.lower() + " electrode"
            elif self.domain == "Separator":
                self.domain_for_broadcast = "separator"

    def get_fundamental_variables(self):
        """
        A public method that creates and returns the variables in a submodel which can
        be created independent of other submodels. For example, the electrolyte
        concentration variables can be created independent of whether any other
        variables have been defined in the model. As a rule, if a variable can be
        created without variables from other submodels, then it should be placed in
        this method.

        Returns
        -------
        dict :
            The variables created by the submodel which are independent of variables in
            other submodels.
        """
        return {}

    def get_external_variables(self):
        """
        A public method that returns the variables in a submodel which are
        supplied by an external source.

        Returns
        -------
        list :
            A list of the external variables in the model.
        """

        external_variables = []
        list_of_vars = []

        if self.external is True:
            # look through all the variables in the submodel and get the
            # variables which are state vectors
            submodel_variables = self.get_fundamental_variables()
            for var in submodel_variables.values():
                if isinstance(var, pybamm.Variable):
                    list_of_vars += [var]

                elif isinstance(var, pybamm.Concatenation):
                    if all(
                        isinstance(child, pybamm.Variable) for child in var.children
                    ):
                        list_of_vars += [var]

            # first add only unique concatenations
            unique_ids = []
            for var in list_of_vars:
                if var.id not in unique_ids and isinstance(var, pybamm.Concatenation):
                    external_variables += [var]
                    unique_ids += [var]
                    # also add the ids of the children to unique ids
                    for child in var.children:
                        unique_ids += [child]

            # now add any unique variables that are not part of a concatentation
            for var in list_of_vars:
                if var.id not in unique_ids:
                    external_variables += [var]
                    unique_ids += [var]

        return external_variables

    def get_coupled_variables(self, variables):
        """
        A public method that creates and returns the variables in a submodel which
        require variables in other submodels to be set first. For example, the
        exchange current density requires the concentration in the electrolyte to
        be created before it can be created. If a variable can be created independent
        of other submodels then it should be created in 'get_fundamental_variables'
        instead of this method.

        Parameters
        ----------
        variables: dict
            The variables in the whole model.

        Returns
        -------
        dict :
            The variables created in this submodel which depend on variables in
            other submodels.
        """
        return {}

    def set_rhs(self, variables):
        """
        A method to set the right hand side of the differential equations which contain
        a time derivative. Note: this method modifies the state of self.rhs. Unless
        overwritten by a submodel, the default behaviour of 'pass' is used as
        implemented in :class:`pybamm.BaseSubModel`.

        Parameters
        ----------
        variables: dict
            The variables in the whole model.
        """
        pass

    def set_algebraic(self, variables):
        """
        A method to set the differential equations which do not contain a time
        derivative. Note: this method modifies the state of self.algebraic. Unless
        overwritten by a submodel, the default behaviour of 'pass' is used as
        implemented in :class:`pybamm.BaseSubModel`.

        Parameters
        ----------
        variables: dict
            The variables in the whole model.
        """
        pass

    def set_boundary_conditions(self, variables):
        """
        A method to set the boundary conditions for the submodel. Note: this method
        modifies the state of self.boundary_conditions. Unless overwritten by a
        submodel, the default behaviour of 'pass' is used a implemented in
        :class:`pybamm.BaseSubModel`.

        Parameters
        ----------
        variables: dict
            The variables in the whole model.
        """
        pass

    def set_initial_conditions(self, variables):
        """
        A method to set the initial conditions for the submodel. Note: this method
        modifies the state of self.initial_conditions. Unless overwritten by a
        submodel, the default behaviour of 'pass' is used a implemented in
        :class:`pybamm.BaseSubModel`.


        Parameters
        ----------
        variables: dict
            The variables in the whole model.
        """
        pass

    def set_events(self, variables):
        """
        A method to set events related to the state of submodel variable. Note: this
        method modifies the state of self.events. Unless overwritten by a submodel, the
        default behaviour of 'pass' is used a implemented in
        :class:`pybamm.BaseSubModel`.

        Parameters
        ----------
        variables: dict
            The variables in the whole model.
        """
        pass
